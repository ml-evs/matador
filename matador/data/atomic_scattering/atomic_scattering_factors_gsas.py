# coding: utf-8
# Distributed under the terms of the MIT License.

from .atomic_scattering_loaders import load_scattering_file

__all__ = ["GSAS_ATOMIC_SCATTERING_COEFFS"]


def load_atomic_scattering_factors_gsas():
    """Load atomic scattering factors used by GSAS-II, from the paper:

    "New Analytical Scattering Factor Functions for Free Atoms and Ions
    for Free Atoms and Ions", D. Waasmaier & A. Kirfel, Acta Cryst.
    (1995). A51, 416-413.

    """
    from pathlib import Path

    filename = Path(__file__).parent.joinpath("atomic_scattering_factors_gsas.dat")
    return load_scattering_file(
        filename, a_inds=(1, 2, 3, 4, 5), b_inds=(7, 8, 9, 10, 11), c_ind=6
    )


GSAS_ATOMIC_SCATTERING_COEFFS = load_atomic_scattering_factors_gsas()
