# S2 maps creation script

# Script to create a s2 table from an .h5 file
import sys
sys.path.append('/scratch/marian.dbt/s2simulation/modules')

from import_modules import *

import s2_table as s2tab


path = '/home/investigator/mariandbt/python/data'
path = path + '/20231025_NEXT100_full_mapping'

list_of_ie_file_paths = []
for i in range(100):
    file_path = os.path.join(path, f'20240228_Next100_ie_s2_{i+1}.next.h5')
    list_of_ie_file_paths.append(file_path)

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

selected_sensors = [204, 240]

for selected_sens in selected_sensors:
    _, bins, maps = s2tab.create_response_maps(list_of_ie_file_paths, selected_sens)

    # Specify the PDF file name
    png_filename = f'sens_{selected_sens}_maps.png'

    s2tab.print_response_maps(selected_sens, bins, maps, yield_)
    plt.savefig(png_filename)
