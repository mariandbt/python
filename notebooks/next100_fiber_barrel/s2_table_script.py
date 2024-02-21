# S2 maps creation script

# Script to create a s2 table from an .h5 file
import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')

from import_modules import *

import s2_table as s2tab


print('START')

path = '/home/investigator/mariandbt/python/data'
path = path + '/20231025_NEXT100_full_mapping'
file_path = os.path.join(path, "Next100_full_mapp_s2_inicioEL_100Kev.next.h5")


bin_width = 10 # [mm]

s2tab.create_s2_table(file_path, bin_width, 20240206) # we use the current date as s2_table_id


print('DONE! :)')
