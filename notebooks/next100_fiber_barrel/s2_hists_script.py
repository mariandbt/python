# S2 dynamic range script

import sys
sys.path.append('/home/investigator/mariandbt/python/notebooks/modules')

from import_modules import *

import set_up as setup
import s2_hists as s2hist


# ________________________________________________________________________________________________________________
# Data reading
# ________________________________________________________________________________________________________________

path = '/home/investigator/mariandbt/python/data'
path = path + '/20231025_NEXT100_full_mapping'


# off_s2_file_path = os.path.join(path, "20240219_bb0nu_200ev_s2_signal.h5") # has fluctuations
off_s2_file_path = os.path.join(path, "20240226_bb0nu_200ev_s2_signal.h5") # has fluctuations

online_s2_filename = os.path.join(path, "next100_fibers/20240122_Next100_bb_1.next.h5") # 1 full bb w s2
sns_positions, sns_response = setup.read_fiber_sens(online_s2_filename)

"""
NOTE: if you want to compare online and offline waveforms
it should be done with the offline file and the original online file from which this was created
"""

# ________________________________________________________________________________________________________________
# Global params
# ________________________________________________________________________________________________________________

t_binnin = 0.1 # [ns] Conversion constant from bin enumerations to nanoseconds (binning used in the simulation)

s2hist.set_global_parameters(globals(), t_binnin)

# ________________________________________________________________________________________________________________
# Analisis
# ________________________________________________________________________________________________________________

s2_max_dict = s2hist.build_offline_s2_max_dict(off_s2_file_path)

pdf_filename = f'dynamic_range_hist.pdf'

# Create a PdfPages object to save the figures in the PDF
with PdfPages(pdf_filename) as pdf:

    events, _, ax = s2hist.print_dyn_range_hist(s2_max_dict, bin_width_in_pes = 250)

    pdf.savefig()
    plt.close()
