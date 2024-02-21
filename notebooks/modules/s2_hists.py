# S2 histograms module

from import_modules import *

import set_up as setup

# ________________________________________________________________________________________________________________
# Functions
# ________________________________________________________________________________________________________________

def online_s2_waveform(filename, event, sensor, bin_width_in_us = 1, new_figure = True, comment = ''):

    font_size = 15
    ev = event
    sens = sensor

    if new_figure:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True) # Create a new figure

    else:
        # Check if there's an existing figure and create it if there's none
        if plt.gcf().get_axes():
            ax = plt.gcf().get_axes()[0]
        else:
            fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), constrained_layout=True)

    if sensor == all:
        online_signal = sns_response.loc[(sns_response.event_id == ev)].copy()
        sens = 'all sensors'


    else:
        online_signal = sns_response.loc[(sns_response.event_id == ev) &
                                         (sns_response.sensor_id == sens)].copy()

    online_signal.time_bin = online_signal.time_bin*t_binning # [ns]

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
                                 label = f'Online s2 of event {0} in {sens} (simulation readout) {comment}'
                                )

    ax.set_title(f's2 waveform for sensor {sens}', fontsize = font_size);
    ax.set_xlabel('Time [us]', fontsize = font_size);
    ax.set_ylabel('Signal [e]', fontsize = font_size);

    ax.tick_params(axis='both', labelsize = font_size*2/3)

    return events, bins, ax





def offline_s2_waveform(filename, event, sensor, bin_width_in_us = 1, new_figure = True, comment = ''):

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
    with h5py.File(filename, 'r') as file:

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
