# coding: utf-8
# Distributed under the terms of the MIT License.

""" The plotting module provides functions to make publication quality
plots for bandstructures, phase diagrams, densities of states, voltages
and volume curves.

Plot style settings can be configured by specifying a matplotlib style
(either built-in or a custom file). The default "matador" style can be
found in the matador.config submodule. User default styles can be set
with the plotting.default_style option in user's matadorrc. Minor tweaks
can be made by passing style dictionaries on a per plot level.

"""


__all__ = [
    "plot_spectral",
    "plot_voltage_curve",
    "plot_volume_curve",
    "plot_2d_hull",
    "plot_ternary_hull",
    "plot_ir_spectrum",
    "get_linear_cmap",
    "plot_free_energy",
    "plot_cutoff_kpt_grid",
    "plot_ensemble_hull",
    "plot_temperature_hull",
    "plot_pdf",
    "plot_projected_pdf",
    "plot_diff_overlap",
    "plot_magres",
    "plot_projected_diff_overlap",
    "plot_pxrd",
    "set_style",
]
__author__ = "Matthew Evans"
__maintainer__ = "Matthew Evans"

from .plotting import get_linear_cmap, set_style
from .spectral_plotting import plot_spectral
from .convergence_plotting import plot_cutoff_kpt_grid
from .battery_plotting import plot_voltage_curve, plot_volume_curve
from .hull_plotting import (
    plot_2d_hull,
    plot_ternary_hull,
    plot_temperature_hull,
    plot_ensemble_hull,
)
from .temperature_plotting import plot_free_energy
from .pdf_plotting import (
    plot_pdf,
    plot_projected_pdf,
    plot_diff_overlap,
    plot_projected_diff_overlap,
)
from .pxrd_plotting import plot_pxrd
from .ir_plotting import plot_ir_spectrum
from .magres_plotting import plot_magres
