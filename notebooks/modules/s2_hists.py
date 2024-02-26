# S2 histograms module

from import_modules import *

import set_up as setup

def set_global_parameters(global_vars, t_binnin_in_ns=None):
    # If parameters are not provided, ask the user for input
    if t_binnin_in_ns is None:
        t_binnin_in_ns = int(input("Specify the time binning used in the simulation: "))

    # Set global variables
    vars_names = ('t_binnin')
    vars_values = (t_binnin_in_ns)

    setup.create_or_update_global_variable(global_vars, vars_names, vars_values, verbose = False)
    setup.create_or_update_global_variable(globals(), vars_names, vars_values, verbose = False)

    # Set more global variables as needed

    print("Global parameters set successfully :)")


def build_offline_s2_max_dict(offline_s2_file_path, bin_width_in_us = 1):

    # Max value of the s2 signals dictionary building

    columns = {0:'time',
               1:'s2'
              }

    bin_width = bin_width_in_us*1000 # [ns] = 1 [us]
    s2_max_dict = {} # max s2 peak per event

    # Open the HDF5 file in read mode
    with h5py.File(offline_s2_file_path, 'r') as file:
        # Iterate through the top-level keys (groups) in the HDF5 file
        for event in file.keys():
            # Get the group corresponding to the current event
            group = file[event]
            s2_max = []

            print(f'Event {event} processed')

            # Iterate through the sensors (datasets) in the current group
            for sensor in group.keys():

                # Get and print the value corresponding to the current sensor
                signal = group[sensor][()]
                signal = pd.DataFrame(signal)
                signal.rename(columns = columns, inplace=True)

                t = signal.time
                s2 = signal.s2
                binin = np.arange(t.min() - bin_width, t.max() + 2*bin_width, bin_width)

                # Create a histogram
                hist_values, bin_edges = np.histogram(t, bins=binin,
                                                      weights = s2)

                s2_max.append(hist_values.max()) # peak of s2 signal per sensor

            s2_max_dict[event] = max(s2_max) # max s2 peak from all sensors

        n_sensors = len(group.keys()) # all events have all sensors, just get the last one

    setup.create_or_update_global_variable(globals(), 'n_sensors', n_sensors, verbose = True)

    return s2_max_dict


def print_online_s2waveform(sns_response, event, sensor, bin_width_in_us = 1, new_figure = True, comment = ''):

    font_size = 15

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
    online_s2 = online_signal.charge # [e]

    t_window_min = 400 # [us]
    t_window_max = 1000 # [us]

    t_window = (t_window_min < tt) & (tt < t_window_max)

    on_t = tt[t_window]
    online_s2 = online_s2[t_window]


    bin_width = bin_width_in_us # time units ([us])

    on_binin = np.arange(on_t.min() - bin_width, on_t.max() + 2*bin_width, bin_width)

    events, bins, bars = ax.hist(on_t, on_binin,
                                 weights = online_s2,
                                 density=False,
                                 histtype='step',
                                 label = f'Online s2 of event {event} in {sensor} (simulation readout) {comment}'
                                )

    ax.set_title(f's2 waveform for sensor {sensor}', fontsize = font_size);
    ax.set_xlabel('Time [us]', fontsize = font_size);
    ax.set_ylabel('Signal [e]', fontsize = font_size);

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax

def print_offline_s2waveform(offline_s2_file_path, event, sensor, bin_width_in_us = 1, new_figure = True, comment = ''):

    if new_figure:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True) # Create a new figure

    else:
        # Check if there's an existing figure and create it if there's none
        if plt.gcf().get_axes():
            ax = plt.gcf().get_axes()[0]
        else:
            fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True)

    font_size = 15
    ev = f'{event}'
    sens = f'sens_{sensor}'

    columns = {0:'time',
               1:'s2'
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
    s2 = signal.s2 # [e]

    bin_width = bin_width_in_us # time units ([us])

    binin = np.arange(t.min() - bin_width, t.max() + 2*bin_width, bin_width)

    events, bins, bars = ax.hist(t, binin,
                                 weights = s2,
                                 density=False,
#                                  histtype='step',
                                 histtype='stepfilled',
                                 alpha = 0.5,
                                 label = f'Offline s2 of event {ev} in {sens} (using maps) {comment}'
                                )

    ax.set_title(f's2 of event {ev} in {sens}', fontsize = font_size);
    ax.set_xlabel('Time [us]', fontsize = font_size);
    ax.set_ylabel('Signal [e]', fontsize = font_size);

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax


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

    ax.text(0.6, .7, '$N_{events}$ = %s'%(int(events.sum())),
    transform=ax.transAxes, fontsize=0.5*font_size, bbox=dict(facecolor='1.', edgecolor='none', pad=3.0))


    ax.set_title(f'Max s2 signal of all {n_sensors} sensors for {n_events} events', fontsize = font_size);
    ax.set_xlabel('s2 signal max [pes]', fontsize = font_size);
    ax.set_ylabel('Counts', fontsize = font_size);
    ax.set_yscale('log')

    ax.legend(fontsize=0.7*font_size, loc='best')

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax
