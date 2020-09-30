# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule provides data loaders for the tabulated electric
quadrupole moment data stored within the module.

"""


def _load_quadrupole_moments_from_file():
    """ Load the electric quadrupole moment file into the constant dictionary
    ELECTRIC_QUADRUPOLE_MOMENTS.

    """

    import numpy as np
    from periodictable import elements
    from pathlib import Path

    filename = Path(__file__).parent.joinpath("electric_quadrupole_moments.csv")
    moments = np.loadtxt(filename, comments="#", delimiter=",")
    moments_dict = {}
    for el, col in zip(elements, moments):
        moments_dict[el.symbol] = moments[el.number-1]

    return moments_dict


ELECTRIC_QUADRUPOLE_MOMENTS = _load_quadrupole_moments_from_file()
