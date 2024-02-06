# S2 maps creation

import set_up


# Functions

def sens_selection(selected_sens):
    """
    Select a sensor and get it's charge distributio for all events.

    Parameters:
    - selected_sens (int): id of selected sensor.

    Returns:
    - sens_response (DataFrame): DataFrame with total charge detected in that sensor for each event.
    """

    # select the sensor
    sens_response = dst.loc[sns_response.sensor_id == selected_sens] # response of the selected sensor

    # fill in the empty events
    ## keep track of the events with no counts
    no_counts_events = set(range(n_events)) - set(sns_response.event_id.unique())
    ## print(f'Events {no_counts_events} did not have any counts')

    ## Create a DataFrame with all event IDs
    all_event_ids = pd.DataFrame({'event_id': range(n_events)})

    ## Merge the two DataFrames to align charge values
    sens_response = all_event_ids.merge(sens_response, on='event_id', how='left').fillna(0)

    return sens_response

def event_coordinates():
    ev_x0 = particles.query('primary == 1').groupby("event_id").initial_x.first()
    ev_y0 = particles.query('primary == 1').groupby("event_id").initial_y.first()
    ev_z0 = particles.query('primary == 1').groupby("event_id").initial_z.first()

    return ev_x0, ev_y0, ev_z0


def create_s2_table(table_name):

    s2_dict = {}

    for ii, selected_sens in enumerate(np.sort(sns_response.sensor_id)):

        print(f'{ii:.0f}/{len(sns_response)}')

        sens_response = sens_selection(selected_sens)

        event_charge = sens_response.groupby("event_id").charge.sum() # total charge detected on that sensor for each event

    # ****************************************************************************************

        # Create a 2D histogram
        hist, xedges, yedges = np.histogram2d(ev_x0, ev_y0,
                                              bins=bins,
                                              weights = event_charge,
                                              density=False);

        hist_counts, xedges, yedges = np.histogram2d(ev_x0, ev_y0,
                                                     bins=bins,
                                                     density=False);

        # Calculate the mean values in each bin (normalized histogram)
        hist_norm = np.where(hist_counts > 0., hist / hist_counts, 0.);

    # ****************************************************************************************
        # Initialize lists to store statistics for each bin
        mean_per_bin = np.zeros((x_nbins, y_nbins))

        # Iterate over each bin
        for i in range(x_nbins):
            for j in range(y_nbins):
                # Indices of data points in the current bin
                mask = ((ev_x0 >= xedges[i]) & (ev_x0 < xedges[i + 1]) &
                        (ev_y0 >= yedges[j]) & (ev_y0 < yedges[j + 1]))

                # Extract values and weights in the current bin
                values_in_bin = event_charge[mask]

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
    with h5py.File('s2_table.h5', 'w') as file:
        for table_id, table_data in s2_dict.items():
            file.create_dataset(table_id, data=table_data.values)

    # ****************************************************************************************




# Data reading

path = '/home/investigator/mariandbt/python/data'
path = path + '/20231025_NEXT100_full_mapping'
file_path = os.path.join(path, "Next100_full_mapp_s2_inicioEL_100Kev.next.h5")

particles = pd.read_hdf(file_path, "/MC/particles")
sns_positions, sns_response = read_sens_h5(file_path)

# Global params

n_sensors = 90
t_binning = 0.1 # [ns] Conversion constant from bin enumerations to nanoseconds (binning used in the simulation)
n_events = particles.event_id.max() + 1 # save number of events simulated
n_sens = len(sns_positions.sensor_id)

bin_width = 10 # [mm]

# Analisis

ev_x0, ev_y0, ev_z0 = event_coordinates()

x_nbins = int((ev_x0.max() - ev_x0.min())/bin_width)
y_nbins = int((ev_y0.max() - ev_y0.min())/bin_width)
bins = (x_nbins, y_nbins)
