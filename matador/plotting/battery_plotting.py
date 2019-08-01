# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot battery-relevant curves,
such as voltages and volume expansions.

"""

import numpy as np
from matador.utils.chem_utils import get_formula_from_stoich
from matador.plotting.plotting import plotting_function, SAVE_EXTS
from matador.plotting.hull_plotting import get_hull_labels


@plotting_function
def plot_voltage_curve(hull, ax=None, show=False, **kwargs):
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

    for ind, (capacities, voltages) in enumerate(zip(hull.voltage_data['Q'], hull.voltage_data['voltages'])):
        if dft_label is None and len(hull.voltage_data['voltages']) > 1:
            stoich_label = get_formula_from_stoich(hull.voltage_data['endstoichs'][ind], tex=True)
        else:
            stoich_label = None
        label = stoich_label if dft_label is None else dft_label
        add_voltage_curve(capacities, voltages, ax_volt, c=hull.colours[ind], label=label)

    if hull.args.get('labels') or hull.args.get('label_cutoff') is not None:
        label_cursor = get_hull_labels(hull, num_species=2)
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


def add_voltage_curve(capacities, voltages, ax_volt, label=None, **kwargs):
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
    import matplotlib.pyplot as plt
    line_kwargs = {'lw': 2,
                   'c': list(plt.rcParams['axes.prop_cycle'].by_key()['color'])[0],
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
def plot_volume_curve(hull, ax=None, show=False, **kwargs):
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

    if ax is None:
        if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
            fig = plt.figure(facecolor=None, figsize=(4, 3.5))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
    else:
        ax = ax

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
    fname = '{}_volume'.format(''.join(hull.elements))

    if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                plt.savefig('{}.{}'.format(fname, ext), transparent=True)
                print('Wrote {}.{}'.format(fname, ext))

    if show:
        plt.show()

def plot_beef_voltage(hull, **kwargs):
    raise NotImplementedError

    beef_cursor = deepcopy(hull.cursor)
    n_beef = len(beef_cursor[0]['_beef']['thetas'])
    if kwargs.get('n_beef') is not None:
        n_beef = min([n_beef, kwargs.get('n_beef')])

    # parameters for voltage heat map
    max_voltage = 2
    min_voltage = 0
    max_q = 1.5*max(hull.voltage_data['Q'][0])
    min_q = 0
    grid_scale = kwargs.get('grid_scale', 100)
    voltage_heatmap = np.zeros((grid_scale, grid_scale), dtype=np.int)
    q_grid, v_grid = np.meshgrid(np.linspace(min_q, max_q, num=grid_scale), np.linspace(min_voltage, max_voltage, num=grid_scale))

    print('Calculating voltage heat map from {} to {} V'.format(min_voltage, max_voltage))

    # collect minimum formation energy for ylimits
    min_ef = 0

    for beef_ind in range(n_beef):
        previous_qindex = 0
        V = beef_hull.voltage_data['voltages'][0]
        Q = beef_hull.voltage_data['Q'][0]
        for step, (qs, vs) in enumerate(zip(Q, V)):
            if step == len(V) - 1 or np.isnan(qs) or np.isnan(vs):
                continue
            v_index = max([0, min([int((vs-min_voltage)/(max_voltage-min_voltage) * grid_scale), grid_scale-1])])
            q_index = max([0, min([int((qs-min_q)/(max_q-min_q)*grid_scale), grid_scale-1])])
            voltage_heatmap[v_index, previous_qindex:q_index+1] += 1
            if kwargs.get('hist_vertical'):
                if not np.isnan(V[step+1]):
                    next_vindex = max([0, min([int(V[step+1]/(max_voltage-min_voltage) * grid_scale), grid_scale-1])])
                    voltage_heatmap[next_vindex+1:v_index, q_index] += 1
            previous_qindex = q_index

        for doc in beef_hull.hull_cursor[1:-1]:
            hull_concs[doc['_formula']] += 1

    ax.set_ylim(min_ef)

    ax_volt.pcolor(q_grid, v_grid, voltage_heatmap, cmap='viridis')#, norm=mpl.colors.LogNorm(vmin=1, vmax=np.max(voltage_heatmap)))

    voltage_curve = calculate_average_voltage_from_heatmap(voltage_heatmap, v_grid, q_grid)
    ax_volt.plot(q_grid[0, :], voltage_curve, c='white', zorder=1e20, lw=3, ls='--')
    c = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])[1]
    add_voltage_curve(hull.voltage_data['Q'][0], hull.voltage_data['voltages'][0], ax_volt, lw=5, ls='--', c=c)
