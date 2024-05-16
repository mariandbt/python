# S2 signal script FULL EVENT
import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')
sys.path.append('/scratch/marian.dbt/s2simulation/modules')

from import_modules import *

import set_up as setup
import s2_bb_events as s2bb
from s2_bb_events import unit


# ________________________________________________________________________________________________________________
# Data reading
# ________________________________________________________________________________________________________________

path    = '/home/investigator/mariandbt/python/data'
path    = path + '/20240313_s2simulation'

# Sensors' info
sns_path = os.path.join(path, "ie/20240405_Next100_ie_s2_1.next.h5")

# Create S2 table
s2_table_path = os.path.join(path, "s2tab/20240405_s2_table.h5")
s2_table = s2bb.s2Table(s2_table_path)

# Create TPC
TPC     = s2bb.FiberBarrelTPC()
TPC.SetActiveDriftVelocity()
TPC.SetRecombinationFactor()
TPC.SetElectronLifetime()
TPC.SetActiveLongDiffusion()
TPC.SetActiveTransDiffusion()
TPC.SetEL()
TPC.SetSensors(sns_path)

# Create nexus event
nexus_event_path    = os.path.join(path, "bb/20240228_Next100_10ev_ELon_bb_1.next.h5")

nexus_event_dict    = {}
n_events            = 10
for ev in range(n_events):
    nexus_event_dict[ev] = s2bb.nexusEvent(nexus_event_path, ev)

event = 5
nexus_event = nexus_event_dict[event]
nexus_event.AddDriftAndDiffusion(TPC)