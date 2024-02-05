# S2 maps creation

import numpy              as np

import data_prep


# Functions

def fill_empty_events(selected_sens):

    sens_response = dst.loc[sns_response.sensor_id == selected_sens] # response of the selected sensor

    # keep track of the events with no counts
    no_counts_events = set(range(n_events)) - set(sns_response.event_id.unique())
    # print(f'Events {no_counts_events} did not have any counts')

    # Create a DataFrame with all event IDs
    all_event_ids = pd.DataFrame({'event_id': range(n_events)})

    # Merge the two DataFrames to align charge values
    sens_response = all_event_ids.merge(sens_response, on='event_id', how='left').fillna(0)


# Data reading

particles = pd.read_hdf(file_path, "/MC/particles")
sns_positions, sns_response = read_sens_h5(file_path)

# Global params

n_sensors = 90
t_binning = 0.1 # [ns] Conversion constant from bin enumerations to nanoseconds (binning used in the simulation)
n_events = particles.event_id.max() + 1 # save number of events simulated

# Data prep
