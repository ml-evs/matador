# coding: utf-8
# Distributed under the terms of the MIT License.

from pathlib import Path

from .atomic_scattering_loaders import load_scattering_file

__all__ = ["RASPA_ATOMIC_SCATTERING_COEFFS"]


def load_atomic_scattering_factors_raspa():
    """Load atomic scattering factors used by RASPA2, with unknown
    reference.

    https://github.com/numat/RASPA2/blob/master/src/scattering_factors.c

    """
    filename = Path(__file__).parent.joinpath("atomic_scattering_factors.dat")
    return load_scattering_file(
        filename, a_inds=(1, 3, 5, 7), b_inds=(2, 4, 6, 8), c_ind=9
    )


RASPA_ATOMIC_SCATTERING_COEFFS = load_atomic_scattering_factors_raspa()
