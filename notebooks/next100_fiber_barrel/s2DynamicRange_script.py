
# S2 Dynamic Range script
import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')

from import_modules import *

import s2_simulation as s2sim
from s2_simulation import unit

# ________________________________________________________________________________________________________________
# Data reading
# ________________________________________________________________________________________________________________

path = '/home/investigator/mariandbt/python/data'
path = path + '/20240313_s2simulation'

# sensors' info
sns_path = os.path.join(path, "ie/20240405_Next100_ie_s2_1.next.h5")

# s2 table
s2_table_path = os.path.join(path, "s2tab/20240405_s2_table.h5")

# bb events
event_path = os.path.join(path, "bb/20240228_Next100_10ev_ELon_bb_1.next.h5") # 10 full bb w s2
# event_path2 = os.path.join(path, "bb/20240306_Next100_200ev_ELoff_bb_1.next.h5") # 200 evs w/o s2

# List_bb_Paths = [event_path, event_path2]
List_bb_Paths = [event_path]

# Create TPC
TPC     = s2sim.HPGXeTPC()
TPC.SetDefaults(sns_path)

# Create s2 table object
s2table = s2sim.s2Table(s2_table_path)

DynamicRange = s2sim.DynamicRange(s2table, 
                                  TPC, 
                                  List_bb_Paths = List_bb_Paths,
                                  fiducial_radio = 490 *unit.mm,
                                  shapin_tau = 155 *unit.ns,
                                  samplin_rate = 25*unit.ns,
                                  t_binin = 1*unit.ns,
                                  units = 'mV',
                                  impedance_in_ohm = 50 # only used if units = mV
                                 )


hdf5output_path     = os.path.join(path, "DynRange/20240617_Next100_bb_10ev_DynRange.h5")
DynamicRange.CreateHDF5(hdf5output_path)