# S2 signal script
import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')

from import_modules import *

import set_up as setup
import s2_signal as s2sig


# ________________________________________________________________________________________________________________
# Data reading
# ________________________________________________________________________________________________________________

path = '/home/investigator/mariandbt/python/data'
path = path + '/20231025_NEXT100_full_mapping'


# s2 table
s2_table_path = os.path.join(path, "s2_table.h5")
s2_table = setup.read_s2_table(s2_table_path)

# bb
# bb_file_path = os.path.join(path, "next100_fibers/20240122_Next100_bb_1.next.h5") # 1 full bb w s2
bb_file_path = os.path.join(path, "next100_fibers/20240111_Next100_bb_3.next.h5") # 200

list_of_bb_file_paths = [bb_file_path]

# ________________________________________________________________________________________________________________
# Global params
# ________________________________________________________________________________________________________________

n_panels = 18
n_sensors = 90

n_bb_files = 1
n_bb_events_per_file = 200

output_file_path = os.path.join(path, '20240222_bb0nu_200ev_s2_signal.h5')

s2sig.set_global_parameters(n_bb_files, n_bb_events_per_file, n_panels, n_sensors)
s2sig.create_s2_signal(s2_table, list_of_bb_file_paths, output_file_path)
