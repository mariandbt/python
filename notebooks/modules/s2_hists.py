# S2 histograms module

from import_modules import *

import set_up as setup
import s2_signal as s2sig

def get_globals():
    return globals()

def set_global_parameters(global_vars, t_binnin_in_ns=None, fiducial_radio_in_mm=None):
    # If parameters are not provided, ask the user for input
    if t_binnin_in_ns is None:
        t_binnin_in_ns = int(input("Specify the time binning used in the simulation: "))
    if fiducial_radio_in_mm is None:
        fiducial_radio_in_mm = int(input("Specify the fiducial radio in mm cut to use: "))


    # Set global variables
    vars_names = ('t_binnin', 'fiducial_radio')
    vars_values = (t_binnin_in_ns, fiducial_radio_in_mm)

    setup.create_or_update_global_variable(global_vars, vars_names, vars_values, verbose = True)
    setup.create_or_update_global_variable(globals(), vars_names, vars_values, verbose = False)

    # Set more global variables as needed

    print("Global parameters set successfully :)")


def print_online_s2waveform(sns_response, event, sensor,
                            bin_width_in_us = 1,
                            t_window_min_in_us = 400, t_window_max_in_us = 1000,
                            new_figure = True, comment = ''):

    font_size = 22

    if new_figure:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True) # Create a new figure

    else:
        # Check if there's an existing figure and create it if there's none
        if plt.gcf().get_axes():
            ax = plt.gcf().get_axes()[0]
        else:
            fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True)

    if sensor == all:
        online_signal = sns_response.loc[(sns_response.event_id == event)].copy()
        sensor = 'all sensors'


    else:
        online_signal = sns_response.loc[(sns_response.event_id == event) &
                                         (sns_response.sensor_id == sensor)].copy()

    online_signal.time_bin = online_signal.time_bin*t_binnin # [ns]

    tt = online_signal.time_bin*1e-3 # [us]
    online_s2 = online_signal.charge # [pes]

    t_window_min = t_window_min_in_us # [us]
    t_window_max = t_window_max_in_us # [us]

    t_window = (t_window_min < tt) & (tt < t_window_max)

    on_t = tt[t_window]
    online_s2 = online_s2[t_window]


    bin_width = bin_width_in_us # time units ([us])

    on_binin = np.arange(on_t.min() - bin_width, on_t.max() + 2*bin_width, bin_width)

    events, bins, bars = ax.hist(on_t, on_binin,
                                 weights = online_s2,
                                 density=False,
                                 histtype='step',
                                 label = f'Online s2 (simulation readout) {comment}'
                                )

    ax.set_title(f's2 waveform of event {event} for sensor {sensor}', fontsize = font_size);
    ax.set_xlabel('Time [us]', fontsize = font_size);
    ax.set_ylabel('Signal [pes]', fontsize = font_size);

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax

def reallyOLD_print_offline_s2waveform(offline_s2_file_path, event, sensor, bin_width_in_us = 1, new_figure = True, comment = ''):

    if new_figure:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True) # Create a new figure

    else:
        # Check if there's an existing figure and create it if there's none
        if plt.gcf().get_axes():
            ax = plt.gcf().get_axes()[0]
        else:
            fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True)

    font_size = 22
    ev = f'{event}'
    sens = f'sens_{sensor}'

    columns = {0:'time',
               1:'s2',
               2:'prim_e_r'
              }


    # Open the HDF5 file in read mode
    with h5py.File(offline_s2_file_path, 'r') as file:

        # Get the group corresponding to the current key
        group = file[ev]

        if sensor == all:

            # Get a list of all keys (sensor names) in the group
            sensor_keys = list(group.keys())

            # Use list comprehension to get all datasets (signals) for all sensors in the group
            all_signals = [group[sensor_key][()] for sensor_key in sensor_keys]

            # Convert the list of signals to a DataFrame
            signal = pd.DataFrame(np.concatenate(all_signals))

            # Assuming 'columns' is defined elsewhere in your code
            signal.rename(columns=columns, inplace=True)

            sens = 'all sensors'


        else:

            # Get and print the value corresponding to the current subkey
            signal = group[sens][()]
            signal = pd.DataFrame(signal)
            signal.rename(columns = columns, inplace=True)


    t = signal.time*1e-3 # [us]
    s2 = signal.s2 # [pes]

    bin_width = bin_width_in_us # time units ([us])

    binin = np.arange(t.min() - bin_width, t.max() + 2*bin_width, bin_width)

    events, bins, bars = ax.hist(t, binin,
                                 weights = s2,
                                 density=False,
#                                  histtype='step',
                                 histtype='stepfilled',
                                 alpha = 0.5,
                                 label = f'Offline s2 (using maps) {comment}'
                                )

    ax.set_title(f's2 of event {ev} in {sens}', fontsize = font_size);
    ax.set_xlabel('Time [us]', fontsize = font_size);
    ax.set_ylabel('Signal [pes]', fontsize = font_size);

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax

def OLD_print_offline_s2waveform(offline_s2_file_path, event, sensor, t0_in_us = 0, bin_width_in_us = 1, new_figure = True, comment = ''):

    if new_figure:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True) # Create a new figure

    else:
        # Check if there's an existing figure and create it if there's none
        if plt.gcf().get_axes():
            ax = plt.gcf().get_axes()[0]
        else:
            fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True)

    font_size = 22
    ev = f'{event}'
    sens = f'sens_{sensor}'

    # Open the HDF5 file in read mode
    with h5py.File(offline_s2_file_path, 'r') as file:

        # Get the group corresponding to the current key
        event_group = file[ev]

        if sensor == all:
            # Get a list of all keys (sensor names) in the group
            sensor_keys = list(event_group.keys())
            s2 = 0

            # Use list comprehension to get all datasets (signals) for all sensors in the group
            for sens_key in sensor_keys:
                signal = event_group[sens_key]
                s2 = s2 + np.array(signal['s2_in_pes']) # [pes]

            samplin_rate = np.array(signal['bin_width_in_ns'])*1e-3 # [us]
            # samplin_rate = np.array(signal['samplin_rate_in_ns'])*1e-3 # [us]
            t = t0_in_us + np.arange(0, len(s2)*samplin_rate, samplin_rate)
            # t = np.array(signal['time_in_ns'])*1e-3 # [us]
            sens = 'all sensors'


        else:

            # Get and print the value corresponding to the current subkey
            signal = event_group[sens]
            s2 = np.array(signal['s2_in_pes']) # [pes]
            samplin_rate = np.array(signal['bin_width_in_ns'])*1e-3 # [us]
            # samplin_rate = np.array(signal['samplin_rate_in_ns'])*1e-3 # [us]
            t = t0_in_us + np.arange(0, len(s2)*samplin_rate, samplin_rate)
            # t = np.array(signal['time_in_ns'])*1e-3 # [us]

    # bin_width = bin_width_in_us # time units ([us])
    bin_width = bin_width_in_us # time units ([us])

    binin = np.arange(t.min() - bin_width, t.max() + 2*bin_width, bin_width)

    events, bins, bars = ax.hist(t, binin,
                                 weights = s2,
                                 density=False,
#                                  histtype='step',
                                 histtype='stepfilled',
                                 alpha = 0.5,
                                 label = f'Offline s2 (using maps) {comment}'
                                )

    ax.set_title(f's2 of event {ev} in {sens}', fontsize = font_size);
    ax.set_xlabel('Time [us]', fontsize = font_size);
    ax.set_ylabel('Signal [pes]', fontsize = font_size);

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax

def print_offline_s2waveform(offline_s2_file_path, event, sensor, t0_in_us = 0, samplin_rate_in_us = 1, new_figure = True, comment = ''):

    if new_figure:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True) # Create a new figure

    else:
        # Check if there's an existing figure and create it if there's none
        if plt.gcf().get_axes():
            ax = plt.gcf().get_axes()[0]
        else:
            fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True)

    font_size = 22
    ev = f'{event}'
    sens = f'sens_{sensor}'


    # Open the HDF5 file in read mode
    with h5py.File(offline_s2_file_path, 'r') as file:

        # Get the group corresponding to the current key
        event_group = file[ev]

        if sensor == all:
            # Get a list of all keys (sensor names) in the group
            sensor_keys = list(event_group.keys())
            s2 = 0

            # Use list comprehension to get all datasets (signals) for all sensors in the group
            for sens_key in sensor_keys:
                signal = event_group[sens_key]
                s2 = s2 + np.array(signal['s2_in_pes']) # [pes]

            original_samplin_rate_in_us = np.array(signal['samplin_rate_in_ns'])*1e-3 # [us]
            sens = 'all sensors'


        else:

            # Get and print the value corresponding to the current subkey
            signal = event_group[sens]
            s2 = np.array(signal['s2_in_pes']) # [pes]
            original_samplin_rate_in_us = np.array(signal['samplin_rate_in_ns'])*1e-3 # [us]

    
    samplin_step = int(samplin_rate_in_us//original_samplin_rate_in_us)

    s2_shaped_sampled  = s2[::samplin_step]
    t_in_us = t0_in_us + np.arange(0, len(s2_shaped_sampled)*samplin_rate_in_us, samplin_rate_in_us)

    ax.plot(t_in_us, s2_shaped_sampled, label = f'Sampling rate of [{samplin_rate_in_us}us] {comment}')


    ax.set_title(f's2 of event {ev} in {sens}', fontsize = font_size);
    ax.set_xlabel('Time [us]', fontsize = font_size);
    ax.set_ylabel('Signal [pes]', fontsize = font_size);

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return s2_shaped_sampled, t_in_us, ax

def print_offline_s2waveform_v2(offline_s2_file_path, event, sensor, bin_width_in_us = 1, new_figure = True, comment = ''):

    if new_figure:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True) # Create a new figure

    else:
        # Check if there's an existing figure and create it if there's none
        if plt.gcf().get_axes():
            ax = plt.gcf().get_axes()[0]
        else:
            fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True)

    font_size = 22
    ev = f'{event}'
    sens = f'sens_{sensor}'


    # Open the HDF5 file in read mode
    with h5py.File(offline_s2_file_path, 'r') as file:

        # Get the group corresponding to the current key
        event_group = file[ev]

        if sensor == all:
            # Get a list of all keys (sensor names) in the group
            sensor_keys = list(event_group.keys())
            s2 = 0

            # Use list comprehension to get all datasets (signals) for all sensors in the group
            for sens_key in sensor_keys:
                signal = event_group[sens_key]
                s2 = s2 + np.array(signal['s2_in_pes']) # [pes]

            sens = 'all sensors'


        else:

            # Get and print the value corresponding to the current subkey
            signal  = event_group[sens]
            s2      = np.array(signal['s2_in_pes']) # [pes]

        t   = np.array(signal['time_in_ns'])*1e-3 # [us]
    
    bin_width = bin_width_in_us # time units ([us])

    binin = np.arange(t.min() - bin_width, t.max() + 2*bin_width, bin_width)

    events, bins, bars = ax.hist(t, binin,
                                 weights = s2,
                                 density=False,
#                                  histtype='step',
                                 histtype='stepfilled',
                                 alpha = 0.5,
                                 label = f'Offline s2 (using maps) {comment}'
                                )

    ax.set_title(f's2 of event {ev} in {sens}', fontsize = font_size);
    ax.set_xlabel('Time [us]', fontsize = font_size);
    ax.set_ylabel('Signal [pes]', fontsize = font_size);

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax

# def OLD_build_offline_s2_max_dict(offline_s2_file_path, bin_width_in_us = 1):

#     # Max value of the s2 signals dictionary building

#     columns = {0:'time',
#                1:'s2',
#                2:'r'
#               }

#     bin_width = bin_width_in_us*1000 # [ns] = 1 [us]
#     s2_max_dict = {} # max s2 peak per event
#     prim_e_r_dict = {} # radial coordinate of each event

#     # Open the HDF5 file in read mode
#     with h5py.File(offline_s2_file_path, 'r') as file:
#         # Iterate through the top-level keys (groups) in the HDF5 file
#         for event in file.keys():
#             # Get the group corresponding to the current event
#             group = file[event]
#             s2_max = []

#             print(f'Event {event} processed')

#             # Iterate through the sensors (datasets) in the current group
#             for sensor in group.keys():

#                 # Get and print the value corresponding to the current sensor
#                 signal = group[sensor][()]
#                 signal = pd.DataFrame(signal)
#                 signal.rename(columns = columns, inplace=True)

#                 # print(signal.r[0])

#                 t = signal.time # [ns]
#                 s2 = signal.s2 # [pes]
#                 prim_e_r = signal.r[0]

#                 if prim_e_r > fiducial_radio:
#                     continue

#                 binin = np.arange(t.min() - bin_width, t.max() + 2*bin_width, bin_width)

#                 # Create a histogram
#                 hist_values, bin_edges = np.histogram(t, bins=binin,
#                                                       weights = s2)

#                 # Shaping
#                 tt = (binin[:-1] + binin[1:])/2 # [ns]

#                 generic_response = s2sig.sipm_response(1, tt, tt.mean())
#                 convolution_response_wvf = np.convolve(hist_values, generic_response, mode='same')

#                 # s2_max.append(hist_values.max()) # peak of s2 signal per sensor
#                 s2_max.append(convolution_response_wvf.max()) # peak of s2 signal per sensor


#             if prim_e_r > fiducial_radio:
#                 print('Discarded event by fiducial cut')
#                 continue
#             s2_max_dict[event] = max(s2_max) # max s2 peak from all sensors
#             prim_e_r_dict[event] = prim_e_r # radial coordinate of each event

#         n_sensors = len(group.keys()) # all events have all sensors, just get the last one

#     setup.create_or_update_global_variable(globals(), 'n_sensors', n_sensors, verbose = True)

#     return s2_max_dict, prim_e_r_dict

def build_offline_s2_max_dict(list_of_offline_s2_file_paths, samplin_rate_in_us = 1):

    # Max value of the s2 signals dictionary building
    s2_max_dict     = {} # max s2 peak per event
    prim_e_r_dict   = {} # radial coordinate of each event
    n_event         = 0
    n_file          = 0
    n_files         = len(list_of_offline_s2_file_paths)

    for offline_s2_file_path in list_of_offline_s2_file_paths:
        n_file   = n_file + 1

        s2_max_dict_this_file       = {}
        prim_e_r_dict_this_file     = {}
        # Open the HDF5 file in read mode
        with h5py.File(offline_s2_file_path, 'r') as file:
            # Iterate through the top-level keys (groups) in the HDF5 file
            n_events    = len(file.keys())
            for event in file.keys():
                # Get the group corresponding to the current event
                event_group = file[event]

                s2_max = []

                print(f'Event {n_event + 1}/{n_events} processed in file {n_file}/{n_files}')

                # Iterate through the sensors (datasets) in the current group
                for sensor in event_group.keys():

                    # Get and print the value corresponding to the current sensor
                    signal = event_group[sensor]

                    prim_e_r = np.array(signal['prim_e_r_in_mm']) # [mm]
                    s2 = np.array(signal['s2_in_pes']) # [pes]
                    original_samplin_rate_in_us = np.array(signal['samplin_rate_in_ns'])*1e-3 # [us]

                    if prim_e_r > fiducial_radio:
                        continue

                    samplin_step = int(samplin_rate_in_us//original_samplin_rate_in_us)

                    s2_shaped_sampled  = s2[::samplin_step]
                    s2_max.append(s2_shaped_sampled.max()) # peak of s2 signal per sensor


                if prim_e_r > fiducial_radio:
                    print('Discarded event by fiducial cut')
                    continue

                s2_max_dict_this_file[n_event] = max(s2_max) # max s2 peak from all sensors
                prim_e_r_dict_this_file[n_event] = prim_e_r # radial coordinate of each event

                n_event = n_event + 1

            s2_max_dict = {**s2_max_dict, **s2_max_dict_this_file}
            prim_e_r_dict = {**prim_e_r_dict, **prim_e_r_dict_this_file}

            n_sensors = len(event_group.keys()) # all events have all sensors, just get the last one

    setup.create_or_update_global_variable(globals(), 'n_sensors', n_sensors, verbose = True)

    return s2_max_dict, prim_e_r_dict


def print_dyn_range_hist(s2_max_dict, bin_width_in_pes = 250):

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7,7), constrained_layout=True)

    font_size = 20

    s2 = np.array(list(s2_max_dict.values()))
    n_events = len(s2)
    binin = np.arange(s2.min() - bin_width_in_pes, s2.max() + 2*bin_width_in_pes, bin_width_in_pes)

    events, bins, bars = ax.hist(s2, binin,
                                 density=False,
                                 label='s2 max value in each event distribution',
                                 histtype='step')


    ax.text(0.6, .85, 'max value =%.2f'%(s2.max()),
    transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

    ax.text(0.6, .8, '$\mu$=%.2f'%(s2.mean()),
    transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

    ax.text(0.6, .75, '$\sigma$=%.2f'%(s2.std()),
    transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

    ax.text(0.6, .7, '$N_{entries}$ = %s'%(int(events.sum())),
    transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))

    ax.text(0.6, .65, 'Fiducial radio cut = %.2f [mm]'%(fiducial_radio),
    transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))


    ax.set_title(f'Max s2 signal of all {n_sensors} sensors for {n_events} events', fontsize = font_size);
    ax.set_xlabel('s2 signal max [pes]', fontsize = font_size);
    ax.set_ylabel('Counts', fontsize = font_size);
    ax.set_yscale('log')

    ax.legend(fontsize=0.7*font_size, loc='best')

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax

