# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot battery-relevant curves,
such as voltages and volume expansions.

"""

from typing import List, Optional, Dict, Union

import numpy as np
from matador.utils.chem_utils import get_formula_from_stoich
from matador.plotting.plotting import plotting_function, SAVE_EXTS
from matador.plotting.hull_plotting import _get_hull_labels
from matador.battery import VoltageProfile


__all__ = ['plot_voltage_curve', 'plot_volume_curve']


@plotting_function
def plot_voltage_curve(
    profiles: Union[List[VoltageProfile], VoltageProfile],
    ax=None,
    show: bool = False,
    savefig: Optional[str] = None,
    curve_labels: Optional[Union[str, List[str]]] = None,
    line_kwargs: Optional[Union[Dict, List[Dict]]] = None,
    expt: Optional[str] = None,
    expt_label: Optional[str] = None
):
    """ Plot voltage curve calculated for phase diagram.

    Parameters:
        profiles (list/VoltageProfile): list of/single voltage profile(s).

    Keyword arguments:
        ax (matplotlib.axes.Axes): an existing axis on which to plot.
        show (bool): whether to show plot in an X window.
        savefig (str): filename to use to save the plot.
        curve_labels (list): optional list of labels for the curves in
            the profiles list.
        line_kwargs (list or dict): parameters to pass to the curve plotter,
            if a list then the line kwargs will be passed to each line individually.
        expt (str): string to a filename of a csv Q, V to add to the plot.
        expt_label (str): label for any experimental profile passed to the plot.

    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax_volt = fig.add_subplot(111)
    else:
        ax_volt = ax

    if not isinstance(profiles, list):
        profiles = [profiles]
    if curve_labels is not None and not isinstance(curve_labels, list):
        curve_labels = [curve_labels]
    if line_kwargs is not None and not isinstance(line_kwargs, list):
        line_kwargs = [line_kwargs]

    if curve_labels is not None and len(curve_labels) != len(profiles):
        raise RuntimeError(
            "Wrong number of labels passed for number of profiles: {} vs {}"
            .format(len(curve_labels), len(profiles))
        )

    if line_kwargs is not None and len(line_kwargs) != len(profiles):
        raise RuntimeError(
            "Wrong number of line kwargs passed for number of profiles: {} vs {}"
            .format(len(line_kwargs), len(profiles))
        )

    dft_label = None
    if expt is not None:
        expt_data = np.loadtxt(expt, delimiter=',')
        if expt_label:
            ax_volt.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label=expt_label)
        else:
            ax_volt.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label='Experiment')

        if len(profiles) == 1:
            dft_label = 'DFT (this work)'

    for ind, profile in enumerate(profiles):
        if dft_label is None and curve_labels is None:
            stoich_label = get_formula_from_stoich(profile.starting_stoichiometry, tex=True)
        else:
            stoich_label = None

        label = stoich_label if dft_label is None else dft_label
        if curve_labels is not None and len(curve_labels) > ind:
            label = curve_labels[ind]

        _line_kwargs = {'c': list(plt.rcParams['axes.prop_cycle'].by_key()['color'])[ind+2]}
        if line_kwargs is not None:
            _line_kwargs.update(line_kwargs[ind])

        _add_voltage_curve(profile.capacities, profile.voltages, ax_volt, label=label, **_line_kwargs)

    if expt or len(profiles) > 1:
        ax_volt.legend()

    ax_volt.set_ylabel('Voltage (V) vs {ion}$^+/${ion}'.format(ion=profile.active_ion))
    ax_volt.set_xlabel('Gravimetric cap. (mAh/g)')

    _, end = ax_volt.get_ylim()
    from matplotlib.ticker import MultipleLocator
    ax_volt.yaxis.set_major_locator(MultipleLocator(0.2))
    ax_volt.set_ylim(0, 1.1 * end)
    _, end = ax_volt.get_xlim()
    ax_volt.set_xlim(0, 1.1 * end)
    ax_volt.grid(False)
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)

    if savefig:
        plt.savefig(savefig)
        print('Wrote {}'.format(savefig))

    elif show:
        plt.show()

    return ax_volt


def _add_voltage_curve(capacities, voltages, ax_volt, label=None, **kwargs):
    """ Add the voltage curves stored under hull['voltage_data'] to the plot.

    Parameters:
        capacities (list): list or numpy array of capacities.
        voltages (list): list or numpy array of voltages.
        ax_volt (matplotlib.axes.Axes): an existing axis object on which to plot.

    Keyword arguments:
        **kwargs (dict): to pass to matplotlib, using abbreviated names (e.g. 'c' not 'color')
        label (str): if present, add this label to the first segment of the voltage curve
            so it can be added to the legend.
        alpha (float): transparency of line

    """
    line_kwargs = {
        'lw': 2,
        'alpha': 1,
    }
    line_kwargs.update(kwargs)

    for i in range(1, len(voltages) - 1):
        if i == 1 and label is not None:
            ax_volt.plot([capacities[i - 1], capacities[i]], [voltages[i], voltages[i]], label=label, **line_kwargs)
        else:
            ax_volt.plot([capacities[i - 1], capacities[i]], [voltages[i], voltages[i]], **line_kwargs)
        if i != len(voltages) - 2:
            ax_volt.plot([capacities[i], capacities[i]], [voltages[i], voltages[i + 1]], **line_kwargs)


@plotting_function
def plot_volume_curve(hull, ax=None, show=True, legend=False, **kwargs):
    """ Plot volume curve calculated for phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        show (bool): whether or not to display plot in X-window.

    """
    import matplotlib.pyplot as plt

    if ax is None:
        if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
            fig = plt.figure(facecolor=None, figsize=(8, 6))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
    else:
        ax = ax

    for j in range(len(hull.volume_data['electrode_volume'])):
        c = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])[j+2]
        stable_hull_dist = hull.volume_data['hull_distances'][j]
        if len(stable_hull_dist) != len(hull.volume_data['Q'][j]):
            raise RuntimeError("This plot does not support --hull_cutoff.")

        ax.plot(
            [q for ind, q in enumerate(hull.volume_data['Q'][j][:-1]) if stable_hull_dist[ind] == 0],
            [v for ind, v in enumerate(hull.volume_data['volume_ratio_with_bulk'][j]) if stable_hull_dist[ind] == 0],
            marker='o', markeredgewidth=1.5, markeredgecolor='k', c=c, zorder=1000, lw=0,
        )

        ax.plot(
            [q for ind, q in enumerate(hull.volume_data['Q'][j][:-1]) if stable_hull_dist[ind] == 0],
            [v for ind, v in enumerate(hull.volume_data['volume_ratio_with_bulk'][j]) if stable_hull_dist[ind] == 0],
            lw=2, c=c,
            label=("Volume expansion from bulk {}"
                   .format(get_formula_from_stoich(hull.volume_data['endstoichs'][j], tex=True)))
        )

    ax.set_xlabel("Gravimetric capacity (mAh/g)")
    ax.set_ylabel('Volume ratio with starting electrode')
    if legend or len(hull.volume_data['Q']) > 1:
        ax.legend()
    fname = '{}_volume'.format(''.join(hull.elements))

    if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                plt.savefig('{}.{}'.format(fname, ext), transparent=True)
                print('Wrote {}.{}'.format(fname, ext))

    if show:
        plt.show()
