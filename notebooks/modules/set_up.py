
# Preparing the data

from import_modules import *


def check_h5(file_path):
    """
    Check an HDF5 file is not empty and see it's tables using tb.open_file().

    Parameters:
    - file_path (str): The path to the HDF5 file.

    Returns:
    - DataFrame: The loaded DataFrame.
    """

    with tb.open_file(file_path) as file:
        print(file)



def read_fiber_sens(file_path):
    """
    Read an HDF5 file using pd.read_hdf() with an optional condition.

    Parameters:
    - file_path (str): The path to the HDF5 file.

    Returns:
    - DataFrame: The loaded DataFrame.
    """

    sns_positions = pd.read_hdf(file_path, "/MC/sns_positions", where='sensor_name == F_SENSOR')
    sns_response = pd.read_hdf(file_path, "/MC/sns_response")


    sns_response = sns_response.loc[sns_response.sensor_id.isin(sns_positions.sensor_id)] # get the positions of said sensors


    return sns_positions, sns_response

def read_s2_table(s2_table_path):
    # Load the 3D dictionary from the HDF5 file
    s2_table = {}

    columns = {0:'bin_initial_x',
               1:'bin_final_x',
               2:'bin_initial_y',
               3:'bin_final_y',
               4:'s2'
              }

    with h5py.File(s2_table_path, 'r') as file:
        for table_id in file.keys():
            # Get the column names from the HDF5 attributes
            s2_table[table_id] = pd.DataFrame(file[table_id][:])
            s2_table[table_id].rename(columns = columns, inplace=True)

    return s2_table

def print_sens_geometry(file_path, selected_sens = None):

    particles = pd.read_hdf(file_path, "/MC/particles")
    sns_positions, _ = read_fiber_sens(file_path)

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(11,11), constrained_layout=True)

    ax.plot(particles.initial_x, particles.initial_y, 'o')


    font_size = 22
    ax.plot(sns_positions.x, sns_positions.y, 'o', markersize = font_size)

    labels_fontsize = font_size

    ax.set_xlabel('X-coordinate [mm]', fontsize = labels_fontsize)
    ax.set_ylabel('Y-coordinate [mm]', fontsize = labels_fontsize)
    ax.tick_params(axis='both', labelsize = labels_fontsize*2/3)

    for sens_id in sns_positions.sensor_id:
        xx = float(sns_positions.loc[sns_positions.sensor_id == sens_id].x)
        yy = float(sns_positions.loc[sns_positions.sensor_id == sens_id].y)

        ax.annotate(f'{sens_id:.0f}', (xx, yy),
                    color='black', ha='center', va='center', fontsize = 0.5*font_size)

    if selected_sens in sns_positions.sensor_id.values:
        ax.plot(sns_positions.loc[sns_positions.sensor_id == selected_sens].x,
                sns_positions.loc[sns_positions.sensor_id == selected_sens].y, 'om',
                markersize = font_size, label = 'Selected sensor')

        ax.legend(fontsize = labels_fontsize, loc = (.4, .6))

    else:
        ax.annotate('Not (valid) sensor ID selected', (0, 100),
                    color='black', ha='center', va='center', fontsize = font_size)


    return ax
