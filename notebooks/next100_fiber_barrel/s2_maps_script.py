# S2 maps creation script

# Script to create a s2 table from an .h5 file
import sys
sys.path.append('/scratch/marian.dbt/s2simulation/modules')
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')

from import_modules import *

import set_up as setup
import s2_table as s2tab


path = '/scratch/marian.dbt/s2simulation/data/ie'
# path = '/home/investigator/mariandbt/python/data/20240313_s2simulation/ie'

list_of_ie_file_paths = []
for i in range(500):
    file_path = os.path.join(path, f'20240405_Next100_ie_s2_{i+1}.next.h5')
    list_of_ie_file_paths.append(file_path)

# file_path = os.path.join(path, f'20240228_Next100_ie_s2_2.next.h5')
# list_of_ie_file_paths.append(file_path)

sns_positions, _ = setup.read_fiber_sens(list_of_ie_file_paths[0])

# Global params
# ________________________________________________________________________________________________________________
x_min = 0
x_max = 490
y_min = -84
y_max = 84

bin_width = 10 # [mm]
yield_ = 1050 # ph/e⁻

s2tab.set_map_specs(globals(), x_min = x_min, x_max = x_max,
                    y_min = y_min, y_max = y_max, bin_width_in_mm = bin_width)
# create and print an example map to check

entries_map_dict = {}

for i, sens_id in enumerate(sns_positions.sensor_id):

    _, bins, maps = s2tab.create_response_maps(list_of_ie_file_paths, sens_id)

    entries_map_dict[sens_id] = {}
    entries_map_dict[sens_id]['bins_x'] = bins[0]
    entries_map_dict[sens_id]['bins_y'] = bins[1]
    entries_map_dict[sens_id]['map'] = maps[3]

    print(f'Sensor {i+1}/{len(sns_positions.sensor_id)}')

# print(entries_map_dict[240])


s2_entries_map_id = 20240412
s2_entries_map_name = f'{s2_entries_map_id}_s2_entries_map.h5'

# Create an HDF5 file
with h5py.File(s2_entries_map_name, 'w') as f:
    # Traverse the dictionary and save each matrix as a dataset
    for key, sub_dict in entries_map_dict.items():
        key = str(key)
        group = f.create_group(key)
        for sub_key, data in sub_dict.items():
            sub_key = str(sub_key)
            group.create_dataset(sub_key, data=data)



print('Done! s2 uncertainty table created :)')
