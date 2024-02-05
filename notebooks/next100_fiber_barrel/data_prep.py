
# Preparing the data

import tables             as tb
import pandas             as pd

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


def read_sens_h5(file_path):
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
