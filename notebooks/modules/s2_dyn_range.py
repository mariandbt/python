# S2 maximums for dynamic range

def s2_max_dict(file_path):
    # Max value of the s2 signals
    # Load the 3D dictionary from the HDF5 file
    columns = {0:'time',
               1:'s2'
              }

    bin_width = 1000 # [ns] = 1 [us]
    s2_max_dict = {} # s2 peak per sensor
    full_s2_max_dict = {} # s2 peak per event

    # Open the HDF5 file in read mode
    with h5py.File(file_path, 'r') as file:
        # Iterate through the top-level keys (groups) in the HDF5 file
        for key in file.keys():
            # Get the group corresponding to the current key
            group = file[key]
            full_s2 = []

            # Print the top-level key
            print(f'Top-level key: {key}')

            # Iterate through the subkeys (datasets) in the current group
            for subkey in group.keys():

                # Check if subkey is already in the dictionary
                if subkey not in s2_max_dict:
                    s2_max_dict[subkey] = []

                # Get and print the value corresponding to the current subkey
                signal = group[subkey][()]
                signal = pd.DataFrame(signal)
                signal.rename(columns = columns, inplace=True)

                t = signal.time
                s2 = signal.s2
                binin = np.arange(t.min() - bin_width, t.max() + 2*bin_width, bin_width)

                # Create a histogram
                hist_values, bin_edges = np.histogram(t, bins=binin,
                                                      weights = s2)

    #             print(f'Subkey: {subkey}, Value: {value}')

                full_s2.append(hist_values.max()) # peak of s2 signal
                s2_max_dict[subkey].append(hist_values.max()) # peak of s2 signal

            full_s2_max_dict[key] = max(full_s2)

    return full_s2_max_dict, s2_max_dict
