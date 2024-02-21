# S2 histograms script

import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')

from import_modules import *

# ________________________________________________________________________________________________________________
# Data reading
# ________________________________________________________________________________________________________________


path = '/home/investigator/mariandbt/python/data'
path = path + '/20231025_NEXT100_full_mapping'


# online s2 simulation
online_s2_filename = os.path.join(path, "next100_fibers/20240122_Next100_bb_1.next.h5") # 1 full bb w s2

# offline s2 reconstruction
off_s2_filename = os.path.join(path, "20240213_bb0nu_1fullev_s2_signal.h5") # has fluctuations
off_s2_filename_no_fluct = os.path.join(path, "20240206_bb0nu_1fullev_s2_signal.h5") # no fluctuations


# ________________________________________________________________________________________________________________
# Global params
# ________________________________________________________________________________________________________________



# ________________________________________________________________________________________________________________
# Analisis
# ________________________________________________________________________________________________________________

sens = 240

ev = '0'
sensor = f'sens_{sens}'
