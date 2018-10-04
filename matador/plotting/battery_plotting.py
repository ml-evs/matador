# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot battery-relevant curves,
such as voltages and volume expansions.

"""

import numpy as np
from matador.utils.chem_utils import get_formula_from_stoich
from matador.plotting.plotting import plotting_function
from matador.plotting.hull_plotting import get_hull_labels


@plotting_function
def plot_voltage_curve(hull, ax=None, show=False):
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
        if hull.savefig:
            if len(hull.voltage_data['voltages']) != 1:
                fig = plt.figure(facecolor=None, figsize=(4, 3.5))
            else:
                fig = plt.figure(facecolor=None, figsize=(4, 3.5))
        else:
            fig = plt.figure(facecolor=None)
        axQ = fig.add_subplot(111)
    else:
        axQ = ax
    if hull.args.get('expt') is not None:
        expt_data = np.loadtxt(hull.args.get('expt'), delimiter=',')
        if hull.args.get('expt_label'):
            axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label=hull.args.get('expt_label'))
        else:
            axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label='Experiment')
    for ind, voltage in enumerate(hull.voltage_data['voltages']):
        for i in range(1, len(voltage) - 1):
            if i == 1 and hull.args.get('expt'):
                axQ.plot([hull.voltage_data['Q'][ind][i - 1], hull.voltage_data['Q'][ind][i]], [voltage[i], voltage[i]],
                         marker='*',
                         lw=2,
                         c=hull.colours[ind],
                         label='DFT (this work)')
            elif i == 1 and len(hull.voltage_data['voltages']) != 1:
                axQ.plot([hull.voltage_data['Q'][ind][i - 1], hull.voltage_data['Q'][ind][i]], [voltage[i], voltage[i]],
                         # marker='o',
                         markersize=5,
                         lw=2,
                         c=hull.colours[ind],
                         label=get_formula_from_stoich(hull.endstoichs[ind], tex=True))
            else:
                axQ.plot([hull.voltage_data['Q'][ind][i - 1], hull.voltage_data['Q'][ind][i]], [voltage[i], voltage[i]],
                         # marker='o',
                         markersize=5,
                         lw=2,
                         c=hull.colours[ind])
            if i != len(voltage) - 2:
                axQ.plot([hull.voltage_data['Q'][ind][i], hull.voltage_data['Q'][ind][i]], [voltage[i], voltage[i + 1]],
                         # marker='o',
                         markersize=5,
                         lw=2,
                         c=hull.colours[ind])
    if hull.args.get('labels') or hull.args.get('label_cutoff') is not None:
        label_cursor = get_hull_labels(hull, num_species=2)
        for i in range(len(label_cursor)):
            axQ.annotate(get_formula_from_stoich(label_cursor[i]['stoichiometry'],
                                                 elements=hull.elements, tex=True),
                         xy=(hull.voltage_data['Q'][0][i+1]+0.02*max(hull.voltage_data['Q'][0]),
                             hull.voltage_data['voltages'][0][i+1]+0.02*max(hull.voltage_data['voltages'][0])),
                         textcoords='data',
                         ha='center',
                         zorder=9999)
    if hull.args.get('expt') or len(hull.voltage_data['voltages']) != 1:
        axQ.legend(loc=1)
    axQ.set_ylabel('Voltage (V) vs {}$^+/${}'.format(hull.elements[0], hull.elements[0]))
    axQ.set_xlabel('Gravimetric cap. (mAh/g)')
    _, end = axQ.get_ylim()
    from matplotlib.ticker import MultipleLocator
    axQ.yaxis.set_major_locator(MultipleLocator(0.2))
    axQ.set_ylim(0, 1.1 * end)
    _, end = axQ.get_xlim()
    axQ.set_xlim(0, 1.1 * end)
    axQ.grid(False)
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)

    if hull.savefig:
        if hull.args.get('pdf'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_voltage.pdf', dpi=500, transparent=True)
        if hull.args.get('svg'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_voltage.svg', dpi=500, transparent=True)
        if hull.args.get('png'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_voltage.png', dpi=500, transparent=True)
    elif show:
        plt.show()

    return axQ


@plotting_function
def plot_volume_curve(hull, show=False):
    """ Plot volume curve calculated for phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        show (bool): whether or not to display plot in X-window.

    """
    import matplotlib.pyplot as plt
    from matador.utils.cursor_utils import get_array_from_cursor
    from matador.utils.chem_utils import get_generic_grav_capacity
    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    if hull.savefig:
        fig = plt.figure(facecolor=None, figsize=(4, 3.5))
    else:
        fig = plt.figure(facecolor=None)
    ax = fig.add_subplot(111)
    stable_hull_dist = get_array_from_cursor(hull.hull_cursor, 'hull_distance')

    hull_vols = []
    hull_comps = []
    for i in range(len(hull.volume_data['vol_per_y'])):
        if stable_hull_dist[i] <= 0 + 1e-16:
            hull_vols.append(hull.volume_data['volume_ratio_with_bulk'][i])
            hull_comps.append(hull.volume_data['x'][i])
            s = 40
            zorder = 1000
            markeredgewidth = 1.5
            c = hull.colours[1]
            alpha = 1
        else:
            s = 30
            zorder = 900
            alpha = 0.3
            markeredgewidth = 0
            c = 'grey'

        ax.scatter(hull.volume_data['x'][i] / (1 + hull.volume_data['x'][i]),
                   hull.volume_data['volume_ratio_with_bulk'][i],
                   marker='o', s=s, edgecolor='k',
                   lw=markeredgewidth, c=c, zorder=zorder, alpha=alpha)

    hull_comps, hull_vols = np.asarray(hull_comps), np.asarray(hull_vols)
    ax.plot(hull_comps / (1 + hull_comps), hull_vols, marker='o', lw=4, c=hull.colours[0], zorder=100)

    ax.set_xlabel(r'$x$ in ' + hull.elements[0] + '$_x$' + hull.elements[1] + '$_{1-x}$')
    ax.set_ylabel('Volume ratio with bulk {}'.format(hull.volume_data['bulk_species']))
    ax.set_ylim(0, np.max(hull.volume_data['volume_ratio_with_bulk']))
    ax.set_xlim(-0.05, 1.05)
    ax.yaxis.set_label_position('left')
    ax.grid(False)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    tick_locs = [0, 0.2, 0.4, 0.6, 0.8, 1]
    ax2.set_xticks(tick_locs)
    new_tick_labels = [int(get_generic_grav_capacity([loc, 1-loc], [hull.elements[0], hull.elements[1]]))
                       for loc in tick_locs[:-1]]
    new_tick_labels[0] = 0
    new_tick_labels.append(r'$\infty$')
    ax2.set_xlabel('Gravimetric capacity (mAh/g)')
    ax2.set_xticklabels(new_tick_labels)
    ax2.grid(False)
    dark_grey = '#262626'
    for spine in ['left', 'top', 'right', 'bottom']:
        ax.spines[spine].set_color(dark_grey)
        ax2.spines[spine].set_color(dark_grey)
        ax.spines[spine].set_linewidth(0.5)
        ax2.spines[spine].set_linewidth(0.5)
    # ax.yaxis.set_ticks(range(0, int(end)+1, 5))
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
    if hull.savefig:
        if hull.args.get('pdf'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_volume.pdf', dpi=300)
        if hull.args.get('svg'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_volume.svg', dpi=300)
        if hull.args.get('png'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_volume.png', dpi=300, bbox_inches='tight')
    elif show:
        plt.show()