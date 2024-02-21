# S2 simulation

from import_modules import *

import set_up as setup


# ________________________________________________________________________________________________________________
# Functions
# ________________________________________________________________________________________________________________



# ________________________________________________________________________________________________________________
# Data reading
# ________________________________________________________________________________________________________________

path = '/home/investigator/mariandbt/python/data'
path = path + '/20231025_NEXT100_full_mapping'


# s2 table
s2_table_path = os.path.join(path, "s2_table.h5")
s2_table = read_s2_table(s2_table_path)

# bb
bb_filename = os.path.join(path, "next100_fibers/20240122_Next100_bb_1.next.h5") # 1 full bb w s2
bb_sns_pos, bb_sns_res = setup.read_fiber_sens(bb_filename)


# ________________________________________________________________________________________________________________
# Global params
# ________________________________________________________________________________________________________________

n_panels = 18
n_sens = 90
dtheta = 2*np.pi/n_panels # rad

dpos = int(n_sens/n_panels) # number of sensors in 1 panel

chunksize = int(2e5) # aprox length of an event to read the tables


# s2 map specs

s2_tab = s2_table[f'sens_{200}'] # all maps have the same specs, so we get whichever

x_nbins, y_nbins = len(s2_tab.bin_initial_x.unique()), len(s2_tab.bin_initial_y.unique())

x_bin_width = (s2_tab.bin_final_x - s2_tab.bin_initial_x)[0]
y_bin_width = (s2_tab.bin_final_y - s2_tab.bin_initial_y)[0]

x_min = s2_tab.bin_initial_x.min()
y_min = s2_tab.bin_initial_y.min()

x_max = s2_tab.bin_final_x.max()
y_max = s2_tab.bin_final_y.max()

fiducial_radio = x_max - 20 # [mm] fiducial cut of 2cm away from the panels


# ________________________________________________________________________________________________________________
# Analisis
# ________________________________________________________________________________________________________________
