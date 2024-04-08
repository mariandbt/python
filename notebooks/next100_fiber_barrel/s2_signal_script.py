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
# s2_table_path = os.path.join(path, "s2_table.h5")
s2_table_path = os.path.join(path, "20240226_s2_table.h5")
s2_table = setup.read_s2_table(s2_table_path)

# bb events
# bb_file_path = os.path.join(path, "next100_fibers/20240122_Next100_bb_1.next.h5") # 1 full bb w s2
# bb_file_path = os.path.join(path, "next100_fibers/20240111_Next100_bb_3.next.h5") # 200

# list_of_bb_file_paths = [bb_file_path]
list_of_bb_file_paths = []
n_bb_files = 10
for i in range(n_bb_files):
    bb_file_path = os.path.join(path, f'20240306_Next100_200ev_ELoff_bb_{i+1}.next.h5')
    list_of_bb_file_paths.append(bb_file_path)

# ________________________________________________________________________________________________________________
# Global params
# ________________________________________________________________________________________________________________

n_panels = 18
n_sensors = 90

v_drift_EL = 2.5e-3 # [mm]/[ns] = 2.5 [mm]/[us]

n_bb_events_per_file = 200

output_file_path = os.path.join(path, '20240226_bb0nu_200ev_s2_signal.h5')

# ________________________________________________________________________________________________________________

s2sig.set_global_parameters(globals(),
                            n_bb_files = n_bb_files,
                            n_bb_events_per_file = n_bb_events_per_file,
                            n_panels = n_panels,
                            n_sensors = n_sensors,
                            v_drift_EL = v_drift_EL
                           )

s2sig.create_s2_signal(s2_table, list_of_bb_file_paths, output_file_path)
