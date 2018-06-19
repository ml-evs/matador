# coding: utf-8
# Distributed under the terms of the MIT License.

""" The plotting module provides functions to make publication quality
plots for bandstructures, phase diagrams, densities of states, voltages
and volume curves.

Plot style settings can be configured by specifying a matplotlib style
(either built-in or a custom file). The default "matador" style can be
found in the matador.config submodule. User default styles can be set
with the plotting.default_style option in user's matadorrc.

"""


__all__ = ['plot_spectral', 'plot_voltage_curve', 'plot_thermo_curves', 'plot_volume_curve',
           'plot_2d_hull', 'plot_ternary_hull', 'get_linear_cmap', 'match_bands']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.plotting.plotting import get_linear_cmap, match_bands
from matador.plotting.plotting import plot_spectral, plot_voltage_curve, plot_volume_curve
from matador.plotting.plotting import plot_thermo_curves, plot_2d_hull, plot_ternary_hull
