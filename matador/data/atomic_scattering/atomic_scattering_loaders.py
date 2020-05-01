# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule provides data loaders for the tabulated atomic
scattering data stored within the module.

"""


def load_scattering_file(filename, a_inds, b_inds, c_ind):
    """ Load a scattering file with a, b and c coefficients
    in the columns provided.

    """

    import numpy as np

    with open(filename, "r") as f:
        flines = [line for line in f.readlines() if not line.startswith("#")]

    atomic_scattering_coeffs = {}

    for line in flines:
        x = line.split()
        label = x[0]
        a = np.array([float(x[ind]) for ind in a_inds])
        b = np.array([float(x[ind]) for ind in b_inds])
        c = float(x[c_ind])
        atomic_scattering_coeffs[label] = [a, b, c]

    return atomic_scattering_coeffs
