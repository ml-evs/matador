# coding: utf-8
# Distributed under the terms of the MIT License.

""" The plotting module provides functions to make publication quality
plots for bandstructures, phase diagrams, densities of states, voltages
and volume curves.

Plot style can be configured partially via settings in matplotlibrc.
Recommended settings:

    font.family      : sans-serif
    font.sans-serif  :  helvetica
    mathtext.rm      : sans
    mathtext.fontset : custom

"""


__all__ = ['plot_spectral', 'plot_voltage_curve', 'plot_thermo_curves', 'plot_volume_curve',
           'plot_2d_hull', 'plot_ternary_hull', 'get_linear_cmap', 'set_seaborn_style']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.plotting.plotting import get_linear_cmap, set_seaborn_style
from matador.plotting.plotting import plot_spectral, plot_voltage_curve, plot_volume_curve
from matador.plotting.plotting import plot_thermo_curves, plot_2d_hull, plot_ternary_hull
