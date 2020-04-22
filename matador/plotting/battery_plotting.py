# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot battery-relevant curves,
such as voltages and volume expansions.

"""

import numpy as np
from matador.utils.chem_utils import get_formula_from_stoich
from matador.plotting.plotting import plotting_function, SAVE_EXTS
from matador.plotting.hull_plotting import _get_hull_labels

__all__ = ['plot_voltage_curve', 'plot_volume_curve']


@plotting_function
def plot_voltage_curve(hull, ax=None, show=False, curve_label=None, line_kwargs=None, **kwargs):
    """ Plot voltage curve calculated for phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        ax (matplotlib.axes.Axes): an existing axis on which to plot.
        show (bool): whether to show plot in an X window.

    """
    import matplotlib.pyplot as plt
    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    if ax is None:
        if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
            fig = plt.figure(facecolor=None, figsize=(4, 3.5))
        else:
            fig = plt.figure(facecolor=None)

        ax_volt = fig.add_subplot(111)
    else:
        ax_volt = ax

    dft_label = None
    if hull.args.get('expt') is not None:
        expt_data = np.loadtxt(hull.args.get('expt'), delimiter=',')
        if hull.args.get('expt_label'):
            ax_volt.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label=hull.args.get('expt_label'))
        else:
            ax_volt.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label='Experiment')

        dft_label = 'DFT (this work)'

    if curve_label is not None:
        dft_label = curve_label

    for ind, (capacities, voltages) in enumerate(zip(hull.voltage_data['Q'], hull.voltage_data['voltages'])):
        if dft_label is None and len(hull.voltage_data['voltages']) > 1:
            stoich_label = get_formula_from_stoich(hull.voltage_data['endstoichs'][ind], tex=True)
        else:
            stoich_label = None
        label = stoich_label if dft_label is None else dft_label
        if line_kwargs is None:
            _line_kwargs = {'c': list(plt.rcParams['axes.prop_cycle'].by_key()['color'])[ind+1]}
        else:
            _line_kwargs = line_kwargs
        _add_voltage_curve(capacities, voltages, ax_volt, label=label, **_line_kwargs)

    if hull.args.get('labels') or hull.args.get('label_cutoff') is not None:
        label_cursor = _get_hull_labels(hull, num_species=2)
        # for testing purposes only
        if 'label_cursor' in kwargs:
            kwargs['label_cursor'].extend(label_cursor)
        for i, doc in enumerate(label_cursor):
            ax_volt.annotate(get_formula_from_stoich(doc['stoichiometry'],
                                                     elements=hull.elements, tex=True),
                             xy=(hull.voltage_data['Q'][0][i+1]+0.02*max(hull.voltage_data['Q'][0]),
                                 hull.voltage_data['voltages'][0][i+1]+0.02*max(hull.voltage_data['voltages'][0])),
                             textcoords='data',
                             ha='center',
                             zorder=9999)

    if hull.args.get('expt') or len(hull.voltage_data['voltages']) != 1:
        ax_volt.legend(loc=1)
    ax_volt.set_ylabel('Voltage (V) vs {}$^+/${}'.format(hull.elements[0], hull.elements[0]))
    ax_volt.set_xlabel('Gravimetric cap. (mAh/g)')
    _, end = ax_volt.get_ylim()
    from matplotlib.ticker import MultipleLocator
    ax_volt.yaxis.set_major_locator(MultipleLocator(0.2))
    ax_volt.set_ylim(0, 1.1 * end)
    _, end = ax_volt.get_xlim()
    ax_volt.set_xlim(0, 1.1 * end)
    ax_volt.grid(False)
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)

    if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
        fname = ''.join(hull.elements) + '_voltage'
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                plt.savefig('{}.{}'.format(fname, ext), transparent=True)
                print('Wrote {}.{}'.format(fname, ext))
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
    line_kwargs = {'lw': 2,
                   'alpha': 1}
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

    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    if ax is None:
        if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
            fig = plt.figure(facecolor=None, figsize=(4, 3.5))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
    else:
        ax = ax

    for j in range(len(hull.volume_data['electrode_volume'])):
        stable_hull_dist = hull.volume_data['hull_distances'][j]
        if len(stable_hull_dist) != len(hull.volume_data['Q'][j]):
            raise RuntimeError("This plot does not support --hull_cutoff.")

        ax.plot(
            [q for ind, q in enumerate(hull.volume_data['Q'][j][:-1]) if stable_hull_dist[ind] == 0],
            [v for ind, v in enumerate(hull.volume_data['volume_ratio_with_bulk'][j]) if stable_hull_dist[ind] == 0],
            marker='o', markeredgewidth=1.5, markeredgecolor='w',
            zorder=100,
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
