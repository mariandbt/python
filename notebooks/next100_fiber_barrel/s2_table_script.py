# S2 maps creation script

# Script to create a s2 table from an .h5 file
import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')

from import_modules import *

import s2_table as s2tab


path = '/home/investigator/mariandbt/python/data'
path = path + '/20231025_NEXT100_full_mapping'
file_path = os.path.join(path, "Next100_full_mapp_s2_inicioEL_100Kev.next.h5")

list_of_ie_file_paths = [file_path]


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

s2tab.create_s2_table(list_of_ie_file_paths, 20240206) # we use the current date as s2_table_id
