# Basic functions

from import_modules import *


def sum_dirac_deltas(x, x0):
    if np.isscalar(x0):  # Check if x0 is a scalar (single number)
        x0 = np.array([x0])  # Convert x0 to a 1-element array
    x0 = x0[:, np.newaxis]  # Convert x0 to a column vector
    all_deltas = np.where(x == x0, 1, 0)
    return all_deltas.sum(axis = 0)