# S2 signal script
import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')
sys.path.append('/scratch/marian/s2simulation/modules')

from import_modules import *

import s2_simulation as s2sim
from s2_simulation import unit

# ________________________________________________________________________________________________________________
# Data reading
# ________________________________________________________________________________________________________________

path    = '/home/investigator/mariandbt/python/data'
path    = path + '/20240313_s2simulation'
# path    = '/scratch/marian/s2simulation/data'

# sensors' info
sns_path    = os.path.join(path, "ie/20240405_Next100_ie_s2_1.next.h5")

# s2 table
s2_table_path   = os.path.join(path, "s2tab/20240405_s2_table.h5")
# s2_table_path   = os.path.join(path, "s2tables/20240405_s2_table.h5")

# bb events
event_path = os.path.join(path, "bb/20240228_Next100_10ev_ELon_bb_1.next.h5") # 10 full bb w s2
# list_of_bb_file_paths   = []
# n_bb_files              = 10
# for i in range(n_bb_files):
#     bb_file_path    = os.path.join(path, f'bb/20240614_Next100_bb_300ev_ELoff_{i+1}.next.h5')
#     list_of_bb_file_paths.append(bb_file_path)

Events_Dict             = {}
# Events_Dict['bb0nu']    = list_of_bb_file_paths
Events_Dict['Kr']     = [event_path]

# Create TPC
TPC     = s2sim.HPGXeTPC()
TPC.SetDefaults(sns_path)

# Create s2 table object
s2table = s2sim.s2Table(s2_table_path, light_yield = 1050, n_ie = 50000)


output_file_path  = os.path.join(path, 's2signals/TEST0.h5')
# output_file_path    = os.path.join(path, 's2signals/20240625_bb0nu_3Kev_ELoff_s2_signal_20240405s2table.h5')
s2sim.CreateSignalHDF5(output_file_path, 
                       s2table, 
                       TPC, 
                       Events_Dict,
                       fiducial_radio     = 490 *unit.mm,
                       shapin_tau         = 155 *unit.ns,
                       samplin_rate       = 25*unit.ns,
                       t_binin            = 1*unit.ns,
                       units              = 'pes/ns',
                       impedance_in_ohm   = 50 # only used if units = mV
                       )
