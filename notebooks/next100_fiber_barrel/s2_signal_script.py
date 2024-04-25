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
path = path + '/20240313_s2simulation'


# s2 table
# s2_table_path = os.path.join(path, "s2_table.h5")
s2_table_path = os.path.join(path, "s2tab/20240405_s2_table.h5")
s2_table = setup.read_s2_table(s2_table_path)

# sensors' info
sns_path = os.path.join(path, "ie/20240405_Next100_ie_s2_1.next.h5")

# bb events
# bb_file_path = os.path.join(path, "bb/20240228_Next100_10ev_ELon_bb_1.next.h5") # 10 full bb w s2
# list_of_bb_file_paths = [bb_file_path]
# n_bb_files = len(list_of_bb_file_paths)

list_of_bb_file_paths = []
n_bb_files = 10
for i in range(n_bb_files):
    bb_file_path = os.path.join(path, f'20240306_Next100_200ev_ELoff_bb_{i+1}.next.h5')
    list_of_bb_file_paths.append(bb_file_path)

# ________________________________________________________________________________________________________________
# Global params
# ________________________________________________________________________________________________________________

n_panels = 18
n_sensors = 108

v_drift_EL = 2.5e-3 # [mm]/[ns] = 2.5 [mm]/[us]

n_bb_events_per_file = 200
# n_bb_events_per_file = 3

output_file_path = os.path.join(path, 's2signals/20240423_TEST3shaped_bb0nu_10ev_ELon_s2_signal_20240405s2table.h5')

# ________________________________________________________________________________________________________________

s2sig.set_global_parameters(globals(),
                            n_bb_files = n_bb_files,
                            n_bb_events_per_file = n_bb_events_per_file,
                            n_panels = n_panels,
                            n_sensors = n_sensors,
                            v_drift_EL = v_drift_EL
                           )

s2sig.create_s2_signal(s2_table, sns_path, list_of_bb_file_paths, output_file_path)
