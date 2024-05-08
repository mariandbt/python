# S2 signal shpain and samplin script
import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')
sys.path.append('/scratch/marian.dbt/s2simulation/modules')

from import_modules import *

import set_up as setup
import s2_signal as s2sig


# ________________________________________________________________________________________________________________
# Data reading
# ________________________________________________________________________________________________________________

path = '/home/investigator/mariandbt/python/data'
path = path + '/20240313_s2simulation'
# path = '/scratch/marian.dbt/s2simulation/data'


# s2 signal NOT shaped NOR sampled
signal_not_shaped_path = os.path.join(path, 's2signals/20240503_bb0nu_2Kev_ELoff_s2_signal_20240405s2table.h5.h5')

# ________________________________________________________________________________________________________________
# Global params
# ________________________________________________________________________________________________________________

shapin_tau_in_ns   = 155 # [ns]
# shapin_tau_in_ns   = 20  # [ns]
samplin_rate_in_ns = 25  # [ns]
t_binin_in_ns      = 1   # [ns]

# ________________________________________________________________________________________________________________

s2sig.shapin_and_samplin(signal_not_shaped_path, shapin_tau_in_ns, samplin_rate_in_ns, t_binin_in_ns)
