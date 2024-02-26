# S2 maps module

from import_modules import *
import set_up as setup

# ________________________________________________________________________________________________________________
# Functions
# ________________________________________________________________________________________________________________

def set_map_specs(global_vars, x_min = None, x_max = None, y_min = None, y_max = None, bin_width_in_mm = None):
    # If parameters are not provided, ask the user for input
    if x_min is None:
        x_min = input("Specify the minimum value of X in [mm] for the maps: ")
    if x_max is None:
        x_max = input("Specify the maximum value of X in [mm] for the maps: ")
    if y_min is None:
        y_min = input("Specify the minimum value of Y in [mm] for the maps: ")
    if y_max is None:
        y_max = input("Specify the maximum value of Y in [mm] for the maps: ")
    if bin_width_in_mm is None:
        bin_width_in_mm = input("Specify the bin width in [mm] for the maps: ")
    # Add more parameters as needed

    bin_width = bin_width_in_mm # [mm]

    x_lims = np.array([x_min, x_max])
    y_lims = np.array([y_min, y_max])

    # Set global variables
    vars_names = ('x_lims', 'y_lims', 'bin_width')
    vars_values = (x_lims, y_lims, bin_width)

    setup.create_or_update_global_variable(global_vars, vars_names, vars_values, verbose = False)
    setup.create_or_update_global_variable(globals(), vars_names, vars_values, verbose = False)

    # Set more global variables as needed

    print("Global parameters set successfully :)")

def sens_selection(sns_response, n_events, selected_sens):
    """
    Select a sensor and get it's charge distributio for all events.

    Parameters:
    - selected_sens (int): id of selected sensor.

    Returns:
    - sens_response (DataFrame): DataFrame with total charge detected in that sensor for each event.
    """

    # select the sensor
    sens_response = sns_response.loc[sns_response.sensor_id == selected_sens] # response of the selected sensor

    # fill in the empty events
    ## keep track of the events with no counts
    no_counts_events = set(range(n_events)) - set(sns_response.event_id.unique())
    ## print(f'Events {no_counts_events} did not have any counts')

    ## Create a DataFrame with all event IDs
    all_event_ids = pd.DataFrame({'event_id': range(n_events)})

    ## Merge the two DataFrames to align charge values
    sens_response = all_event_ids.merge(sens_response, on='event_id', how='left').fillna(0)

    return sens_response

def event_coordinates(particles):
    ev_x0 = particles.query('primary == 1').groupby("event_id").initial_x.first()
    ev_y0 = particles.query('primary == 1').groupby("event_id").initial_y.first()
    ev_z0 = particles.query('primary == 1').groupby("event_id").initial_z.first()

    return ev_x0, ev_y0, ev_z0

def map_binnin(x, y, bin_width):

    x_nbins = int((x.max() - x.min())/bin_width)
    y_nbins = int((y.max() - y.min())/bin_width)

    xedges = np.linspace(x.min(), x.max(), x_nbins + 1)
    yedges = np.linspace(y.min(), y.max(), y_nbins + 1)

    return x_nbins, y_nbins, xedges, yedges


def create_response_maps(list_of_ie_file_paths, selected_sens):


    x_nbins, y_nbins, xedges, yedges = map_binnin(x_lims, y_lims, bin_width)

    values_in_bin_dict = {}

    # Initialize lists to store statistics for each bin
    mean_per_bin = np.zeros((x_nbins, y_nbins))
    std_per_bin = np.zeros((x_nbins, y_nbins))
    entries_per_bin = np.zeros((x_nbins, y_nbins))
    mean_per_bin_err = np.zeros((x_nbins, y_nbins))

    for file_path in list_of_ie_file_paths:

        # ________________________________________________________________________________________________________________
        # Data reading
        # ________________________________________________________________________________________________________________

        particles = pd.read_hdf(file_path, "/MC/particles")
        sns_positions, sns_response = setup.read_fiber_sens(file_path)

        # ________________________________________________________________________________________________________________


        n_events = particles.event_id.max() + 1 # save number of events simulated

        ev_x0, ev_y0, ev_z0 = event_coordinates(particles)

        sens_response = sens_selection(sns_response, n_events, selected_sens)

        event_charge = sens_response.groupby("event_id").charge.sum() # total charge detected on that sensor for each event


        # Iterate over each bin
        for i in range(x_nbins):
            if i not in values_in_bin_dict.keys():
                values_in_bin_dict[i] = {}

            for j in range(y_nbins):

                if j not in values_in_bin_dict[i].keys():
                    values_in_bin_dict[i][j] = np.array([])

                # Indices of data points in the current bin
                mask = ((ev_x0 >= xedges[i]) & (ev_x0 < xedges[i + 1]) &
                        (ev_y0 >= yedges[j]) & (ev_y0 < yedges[j + 1]))

                # Extract values and weights in the current bin
                values_in_bin_dict[i][j] = np.concatenate((values_in_bin_dict[i][j], event_charge[mask]))


    # Iterate over each bin
    for i in range(x_nbins):
        for j in range(y_nbins):
            # Calculate weighted mean and standard deviation
            mean_value = np.mean(values_in_bin_dict[i][j])
            std_value = np.std(values_in_bin_dict[i][j])
            entries_value = len(values_in_bin_dict[i][j])

            # Append to lists
            mean_per_bin[i][j] = mean_value
            std_per_bin[i][j] = std_value
            entries_per_bin[i][j] = entries_value
            mean_per_bin_err[i][j] = np.where(entries_value*mean_value> 0.,
                                              std_value*100/(np.sqrt(entries_value)*mean_value),
                                              0.);

    mean_per_bin = np.nan_to_num(mean_per_bin)
    std_per_bin = np.nan_to_num(std_per_bin)
    mean_per_bin_err = np.nan_to_num(mean_per_bin_err)
    entries_per_bin = np.nan_to_num(entries_per_bin)

    print("""
            Maps created! Access the different maps following:
            maps[0] = mean_per_bin
            maps[1] = std_per_bin
            maps[2] = mean_per_bin_err
            maps[3] = entries_per_bin
          """
    )

    return selected_sens, (xedges, yedges), (mean_per_bin, std_per_bin, mean_per_bin_err, entries_per_bin)


def print_sensor_map(particles, sns_positions, selected_sens, bins, mean_per_bin):


    gs_kw = {'height_ratios': [1, 1.2], 'width_ratios': [1]}
    fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize=(15,13), gridspec_kw=gs_kw, constrained_layout=True)

    font_size = 22
    offset = 0.

    xedges, yedges = bins

    # sector plot
    ax[0].plot(particles.initial_x, particles.initial_y, 'o')

    # sensors positions plot
    ax[0].plot(sns_positions.x, sns_positions.y, 'o', markersize = 0.5*font_size)

    # selected sensor position
    this_sensor = sns_positions.loc[(sns_positions.sensor_id == selected_sens)]
    ax[0].plot(this_sensor.x, this_sensor.y, 'o', markersize = 0.5*font_size, label = f"sens_{selected_sens}")

    ax[0].legend(loc = (.25, .75), fontsize = 0.7*font_size)
    ax[0].set_aspect("equal")

    # ****************************************************************************************
    mean_map = ax[1].pcolormesh(xedges, yedges, mean_per_bin.T, cmap='inferno');
    fig.colorbar(mean_map, ax = ax[1])
    ax[1].set_title(r'$\mu$ of the detected charge [pes]', fontsize = font_size);

    ax[1].set_xlim(xedges.min() - offset, xedges.max() + offset)
    ax[1].set_ylim([yedges.min() - offset, yedges.max() + offset])

    # Annotate each cell with its numeric value
    for i in range(len(xedges) - 1):
        for j in range(len(yedges) - 1):
            value = mean_per_bin[i, j]
            if value > 0:
                ax[1].annotate(f'{value:.2f}', ((xedges[i] + xedges[i + 1]) / 2, (yedges[j] + yedges[j + 1]) / 2),
                            color='white', ha='center', va='center', fontsize = .3*font_size)


    # ****************************************************************************************
    for axx in ax:
        axx.set_xlabel('X-coordinate [mm]', fontsize = font_size)
        axx.set_ylabel('Y-coordinate [mm]', fontsize = font_size)
        axx.tick_params(axis='both', labelsize = font_size*2/3)

    #     axx.set_aspect("equal")



def print_response_maps(selected_sens, bin_edges, maps, yield_):

    fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize=(10,5), constrained_layout=True)
    font_size = 11
    offset = 0.

    xedges, yedges = bin_edges
    mean_per_bin, std_per_bin, mean_per_bin_err, entries_per_bin = maps


    # mean
    mean_map = ax[0][0].pcolormesh(xedges, yedges, mean_per_bin.T, cmap='inferno');
    fig.colorbar(mean_map, ax = ax[0][0])
    ax[0][0].set_title(r'$\mu$ of the detected photons', fontsize = font_size);

    # std
    std_map = ax[0][1].pcolormesh(xedges, yedges, std_per_bin.T, cmap='inferno');
    fig.colorbar(std_map, ax = ax[0][1])
    ax[0][1].set_title(r'$\sigma$ of the detected photons', fontsize = font_size);

    # mean_error
    err_map = ax[1][0].pcolormesh(xedges, yedges, mean_per_bin_err.T, cmap='inferno');
    fig.colorbar(err_map, ax = ax[1][0])
    ax[1][0].set_title(r'$\Delta\mu$ of the detected photons [%]', fontsize = font_size);

    # Entries
    entries_map = ax[1][1].pcolormesh(xedges, yedges, entries_per_bin.T*yield_, cmap='inferno');
    fig.colorbar(entries_map, ax = ax[1][1])
    ax[1][1].set_title(r'$N_{\gamma}$ generated ($EL_{yield}$ = %s photons/e⁻)'%(yield_), fontsize = font_size);

    for row in ax:
        for axis in row:
            axis.set_xlim(xedges.min() - offset, xedges.max() + offset)
            axis.set_ylim([yedges.min() - offset, yedges.max() + offset])
            axis.set_xlabel('X-coordinate [mm]', fontsize = font_size)
            axis.set_ylabel('Y-coordinate [mm]', fontsize = font_size)
            axis.set_aspect("equal")


    fig.suptitle(f'Maps for sens_{selected_sens}', fontsize = 1.5*font_size);



def create_s2_table(list_of_ie_file_paths, s2_table_id):

    x_nbins, y_nbins, xedges, yedges = map_binnin(x_lims, y_lims, bin_width)
    s2_dict = {}

    values_in_bin_dict = {}

    # Initialize lists to store statistics for each bin
    mean_per_bin = np.zeros((x_nbins, y_nbins))
    std_per_bin = np.zeros((x_nbins, y_nbins))
    entries_per_bin = np.zeros((x_nbins, y_nbins))
    mean_per_bin_err = np.zeros((x_nbins, y_nbins))

    sns_positions, _ = setup.read_fiber_sens(list_of_ie_file_paths[0])
    n_sens = len(sns_positions.sensor_id)

    # ****************************************************************************************
    for ff, file_path in enumerate(list_of_ie_file_paths):

        print(f'Processing file {(ff+1):.0f}/{len(list_of_ie_file_paths)}')


        # ________________________________________________________________________________________________________________
        # Data reading
        # ________________________________________________________________________________________________________________

        particles = pd.read_hdf(file_path, "/MC/particles")
        _, sns_response = setup.read_fiber_sens(file_path)

        # ________________________________________________________________________________________________________________

        n_events = particles.event_id.max() + 1 # save number of events simulated

        ev_x0, ev_y0, _ = event_coordinates(particles)

        for ii, selected_sens in enumerate(np.sort(sns_positions.sensor_id)):

            if selected_sens not in values_in_bin_dict.keys():
                values_in_bin_dict[selected_sens] = {}


            sens_response = sens_selection(sns_response, n_events, selected_sens)

            event_charge = sens_response.groupby("event_id").charge.sum() # total charge detected on that sensor for each event

            # Iterate over each bin
            for i in range(x_nbins):
                if i not in values_in_bin_dict[selected_sens].keys():
                    values_in_bin_dict[selected_sens][i] = {}

                for j in range(y_nbins):

                    if j not in values_in_bin_dict[selected_sens][i].keys():
                        values_in_bin_dict[selected_sens][i][j] = np.array([])

                    # Indices of data points in the current bin
                    mask = ((ev_x0 >= xedges[i]) & (ev_x0 < xedges[i + 1]) &
                            (ev_y0 >= yedges[j]) & (ev_y0 < yedges[j + 1]))

                    # Extract values and weights in the current bin
                    values_in_bin_dict[selected_sens][i][j] = np.concatenate((values_in_bin_dict[selected_sens][i][j], event_charge[mask]))
    # ****************************************************************************************

    print(f'Creating table ...')

    for ii, selected_sens in enumerate(np.sort(sns_positions.sensor_id)):
        print(f'Sensor {(ii+1):.0f}/{n_sens}')

        # Initialize lists to store statistics for each bin
        mean_per_bin = np.zeros((x_nbins, y_nbins))

        # Iterate over each bin
        for i in range(x_nbins):
            for j in range(y_nbins):

                # Extract values and weights in the current bin
                values_in_bin = values_in_bin_dict[selected_sens][i][j]

                # Calculate weighted mean and standard deviation
                mean_value = np.mean(values_in_bin)

                # Append to lists
                mean_per_bin[i][j] = mean_value

        mean_per_bin = np.nan_to_num(mean_per_bin)
    # ****************************************************************************************
        table_id = f'sens_{selected_sens}'
        table_data = pd.DataFrame(columns=[])

        for i in range(len(xedges) - 1):
            for j in range(len(yedges) - 1):
                bin_x0 = xedges[i]
                bin_xf = xedges[i+1]
                bin_y0 = yedges[j]
                bin_yf = yedges[j+1]
                s2 = mean_per_bin[i, j]

                new_row = {'bin_x0':bin_x0, 'bin_xf':bin_xf,
                           'bin_y0':bin_y0, 'bin_yf':bin_yf,
                           's2':s2
                          }

                table_data = table_data.append(new_row, ignore_index=True)
    #             print(table_data)

        s2_dict[table_id] = table_data


    # Save the 3D dictionary using HDF5 format
    with h5py.File(f"{s2_table_id}_s2_table.h5", 'w') as file:
        for table_id, table_data in s2_dict.items():
            file.create_dataset(table_id, data=table_data.values)

    # ****************************************************************************************

    print('Done! s2 table created :)')



# def create_s2_table(file_path, bin_width, s2_table_id):
#
#     # ________________________________________________________________________________________________________________
#     # Data reading
#     # ________________________________________________________________________________________________________________
#
#     particles = pd.read_hdf(file_path, "/MC/particles")
#     sns_positions, sns_response = setup.read_fiber_sens(file_path)
#
#     # ________________________________________________________________________________________________________________
#
#     n_sens = len(sns_positions.sensor_id)
#     n_events = particles.event_id.max() + 1 # save number of events simulated
#
#     ev_x0, ev_y0, ev_z0 = event_coordinates(particles)
#     x_nbins, y_nbins, xedges, yedges = map_binnin(ev_x0, ev_y0, bin_width)
#
#     s2_dict = {}
#
#     for ii, selected_sens in enumerate(np.sort(sns_positions.sensor_id)):
#
#         print(f'{ii:.0f}/{n_sens}')
#
#         sens_response = sens_selection(sns_response, n_events, selected_sens)
#
#         event_charge = sens_response.groupby("event_id").charge.sum() # total charge detected on that sensor for each event
#
#     # ****************************************************************************************
#         # Initialize lists to store statistics for each bin
#         mean_per_bin = np.zeros((x_nbins, y_nbins))
#
#         # Iterate over each bin
#         for i in range(x_nbins):
#             for j in range(y_nbins):
#                 # Indices of data points in the current bin
#                 mask = ((ev_x0 >= xedges[i]) & (ev_x0 < xedges[i + 1]) &
#                         (ev_y0 >= yedges[j]) & (ev_y0 < yedges[j + 1]))
#
#                 # Extract values and weights in the current bin
#                 values_in_bin = event_charge[mask]
#
#                 # Calculate weighted mean and standard deviation
#                 mean_value = np.mean(values_in_bin)
#
#                 # Append to lists
#                 mean_per_bin[i][j] = mean_value
#
#         mean_per_bin = np.nan_to_num(mean_per_bin)
#     # ****************************************************************************************
#         table_id = f'sens_{selected_sens}'
#         table_data = pd.DataFrame(columns=[])
#
#         for i in range(len(xedges) - 1):
#             for j in range(len(yedges) - 1):
#                 bin_x0 = xedges[i]
#                 bin_xf = xedges[i+1]
#                 bin_y0 = yedges[j]
#                 bin_yf = yedges[j+1]
#                 s2 = mean_per_bin[i, j]
#
#                 new_row = {'bin_x0':bin_x0, 'bin_xf':bin_xf,
#                            'bin_y0':bin_y0, 'bin_yf':bin_yf,
#                            's2':s2
#                           }
#
#                 table_data = table_data.append(new_row, ignore_index=True)
#     #             print(table_data)
#
#         s2_dict[table_id] = table_data
#
#
#     # Save the 3D dictionary using HDF5 format
#     with h5py.File(f"{s2_table_id}_s2_table.h5", 'w') as file:
#         for table_id, table_data in s2_dict.items():
#             file.create_dataset(table_id, data=table_data.values)
#
#     # ****************************************************************************************
#
#     print('Done! s2 table created :)')
