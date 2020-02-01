# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot convex hulls and phase
diagrams generally.

"""


import numpy as np
from matador.utils.chem_utils import get_stoich_from_formula, get_formula_from_stoich
from matador.plotting.plotting import plotting_function, get_linear_cmap, SAVE_EXTS

EPS = 1e-12


def get_hull_labels(hull, label_cutoff=None, num_species=None, exclude_edges=True):
    """ Return list of structures to labels on phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): phase diagram to plot.

    Keyword arguments:
        num_species (int): structures containing this number of species
            will be labelled.
        exclude_edges (bool): ignore any structure that has zero of any
            chemical potential (i.e the edge of the hull).

    Returns:
        label_cursor (list(dict)): list of matador documents to label.

    """
    eps = 1e-9
    if label_cutoff is None:
        label_cutoff = 0.0
    if isinstance(label_cutoff, list) and len(label_cutoff) == 2:
        label_cutoff = sorted(label_cutoff)
        # first, only apply upper limit as we need to filter by stoich aftewards
        label_cursor = [doc for doc in hull.cursor if doc['hull_distance'] <= label_cutoff[1]]
    else:
        if isinstance(label_cutoff, list):
            assert len(label_cutoff) == 1, 'Incorrect number of label_cutoff values passed, should be 1 or 2.'
            label_cutoff = label_cutoff[0]
        label_cursor = [doc for doc in hull.cursor if doc['hull_distance'] <= label_cutoff + eps]

    num_labels = len({get_formula_from_stoich(doc['stoichiometry']) for doc in label_cursor})
    if num_labels < len(label_cursor):
        tmp_cursor = []
        for doc in label_cursor:
            if doc['stoichiometry'] not in [_doc['stoichiometry'] for _doc in tmp_cursor]:
                tmp_cursor.append(doc)
            else:
                label_cursor = tmp_cursor
    if isinstance(label_cutoff, list) and len(label_cutoff) == 2:
        # now apply lower limit
        label_cursor = [doc for doc in label_cursor if label_cutoff[0] <= doc['hull_distance'] <= label_cutoff[1]]
    # remove chemical potentials and unwanted e.g. binaries
    if num_species is not None:
        label_cursor = [doc for doc in label_cursor if len(doc['stoichiometry']) == num_species]
    if exclude_edges:
        label_cursor = [doc for doc in label_cursor if (all(doc['concentration']) > 0 and
                                                        sum(doc['concentration']) <= 1 - 1e-6)]

    label_cursor = sorted(label_cursor, key=lambda doc: doc['concentration'])

    return label_cursor


@plotting_function
def plot_2d_hull(hull, ax=None, show=True, plot_points=True, plot_tie_line=True,
                 plot_hull_points=True, labels=None, label_cutoff=None, colour_by_source=False,
                 sources=None, hull_label=None, source_labels=None, title=True, plot_fname=None, show_cbar=True,
                 label_offset=(1.15, 0.05), eform_limits=None,
                 **kwargs):
    """ Plot calculated hull, returning ax and fig objects for further editing.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        ax (matplotlib.axes.Axes): an existing axis on which to plot,
        show (bool): whether or not to display the plot in an X window,
        plot_points (bool): whether or not to display off-hull structures,
        plot_hull_points (bool): whether or not to display on-hull structures,
        labels (bool): whether to label formulae of hull structures, also read from
            hull.args.
        label_cutoff (float/:obj:`tuple` of :obj:`float`): draw labels less than or
            between these distances form the hull, also read from hull.args.
        colour_by_source (bool): plot and label points by their sources
        alpha (float): alpha value of points when colour_by_source is True
        sources (list): list of possible provenances to colour when colour_by_source
            is True (others will be grey)
        title (str/bool): whether to include a plot title.
        png/pdf/svg (bool): whether or not to write the plot to a file.
        plot_fname (str): filename to write plot to, without file extension.

    Returns:
        matplotlib.axes.Axes: matplotlib axis with plot.

    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as colours

    if ax is None:
        fig = plt.figure(facecolor=None, figsize=(8, 6))
        ax = fig.add_subplot(111, facecolor=None)

    if not hasattr(hull, 'colours'):
        hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    hull.default_cmap_list = get_linear_cmap(hull.colours[1:4], list_only=True)
    hull.default_cmap = get_linear_cmap(hull.colours[1:4], list_only=False)

    if labels is None:
        labels = hull.args.get('labels', False)
    if label_cutoff is None:
        label_cutoff = hull.args.get('label_cutoff')

    scale = 1
    scatter = []
    chempot_labels = [get_formula_from_stoich(get_stoich_from_formula(species, sort=False), tex=True)
                      for species in hull.species]
    tie_line = hull.structure_slice[hull.convex_hull.vertices]

    # plot hull structures
    if plot_hull_points:
        ax.scatter(tie_line[:, 0], tie_line[:, 1],
                   c=hull.colours[1],
                   marker='o', zorder=99999, edgecolor='k',
                   s=scale*40, lw=1.5)
        if plot_tie_line:
            ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1],
                    c=hull.colours[0], zorder=1, label=hull_label,
                    marker='o', markerfacecolor=hull.colours[0],
                    markeredgecolor='k', markeredgewidth=1.5, markersize=np.sqrt(scale*40))
    if plot_tie_line:
        ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1],
                c=hull.colours[0], zorder=1, label=hull_label, markersize=0)

    if hull.hull_cutoff > 0:
        ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1] + hull.hull_cutoff,
                '--', c=hull.colours[1], alpha=0.5, zorder=1, label='')

    # annotate hull structures
    if labels or label_cutoff is not None:
        label_cursor = get_hull_labels(hull, num_species=2, label_cutoff=label_cutoff)
        already_labelled = []
        for ind, doc in enumerate(label_cursor):
            formula = get_formula_from_stoich(doc['stoichiometry'], sort=True)
            if formula not in already_labelled:
                arrowprops = dict(arrowstyle="-|>", lw=2, alpha=1, zorder=1, shrinkA=2, shrinkB=4)
                min_comp = tie_line[np.argmin(tie_line[:, 1]), 0]
                e_f = label_cursor[ind]['formation_' + str(hull.energy_key)]
                conc = label_cursor[ind]['concentration'][0]
                if conc < min_comp:
                    position = (0.8 * conc, label_offset[0] * (e_f - label_offset[1]))
                elif label_cursor[ind]['concentration'][0] == min_comp:
                    position = (conc, label_offset[0] * (e_f - label_offset[1]))
                else:
                    position = (min(1.1 * conc + 0.15, 0.95), label_offset[0] * (e_f - label_offset[1]))
                ax.annotate(get_formula_from_stoich(doc['stoichiometry'],
                                                    latex_sub_style=r'\mathregular',
                                                    tex=True,
                                                    elements=hull.species,
                                                    sort=False),
                            xy=(conc, e_f),
                            xytext=position,
                            textcoords='data',
                            ha='right',
                            va='bottom',
                            arrowprops=arrowprops,
                            zorder=1)
                already_labelled.append(formula)

    # points for off hull structures; we either colour by source or by energy
    if plot_points and not colour_by_source:

        if hull.hull_cutoff == 0:
            # if no specified hull cutoff, ignore labels and colour by hull distance
            cmap = hull.default_cmap
            if plot_points:
                scatter = ax.scatter(hull.structures[np.argsort(hull.hull_dist), 0][::-1],
                                     hull.structures[np.argsort(hull.hull_dist), -1][::-1],
                                     s=scale*40,
                                     c=np.sort(hull.hull_dist)[::-1],
                                     zorder=10000,
                                     cmap=cmap, norm=colours.LogNorm(0.02, 2))

                if show_cbar:
                    cbar = plt.colorbar(scatter, aspect=30, pad=0.02,
                                        ticks=[0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                    cbar.ax.tick_params(length=0)
                    cbar.ax.set_yticklabels([0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                    cbar.ax.yaxis.set_ticks_position('right')
                    cbar.ax.set_frame_on(False)
                    cbar.outline.set_visible(False)
                    cbar.set_label('Distance from hull (eV/atom)')

        elif hull.hull_cutoff != 0:
            # if specified hull cutoff colour those below
            c = hull.colours[1]
            for ind in range(len(hull.structures)):
                if hull.hull_dist[ind] <= hull.hull_cutoff or hull.hull_cutoff == 0:
                    if plot_points:
                        scatter.append(ax.scatter(hull.structures[ind, 0], hull.structures[ind, 1],
                                                  s=scale*40,
                                                  alpha=0.9, c=c,
                                                  zorder=300))
            if plot_points:
                ax.scatter(hull.structures[1:-1, 0], hull.structures[1:-1, 1],
                           s=scale*30, lw=0,
                           alpha=0.3, c=hull.colours[-2],
                           edgecolor='k', zorder=10)

    elif colour_by_source:
        _scatter_plot_by_source(
            hull, ax, scale, kwargs,
            sources=None, source_labels=None, plot_hull_points=plot_hull_points)

    if eform_limits is None:
        eform_limits = (np.min(hull.structures[:, 1]), np.max(hull.structures[:, 1]))
        lims = (-0.1 if eform_limits[0] >= 0 else 1.4*eform_limits[0],
                eform_limits[1] if eform_limits[0] >= 0 else 0.1)
    else:
        lims = sorted(eform_limits)
    ax.set_ylim(lims)

    if isinstance(title, bool) and title:
        if hull._non_elemental:
            ax.set_title(r'({d[0]})$_\mathrm{{x}}$({d[1]})$_\mathrm{{1-x}}$'.format(d=chempot_labels))
        else:
            ax.set_title(r'{d[0]}$_\mathrm{{x}}${d[1]}$_\mathrm{{1-x}}$'.format(d=chempot_labels))
    elif isinstance(title, str) and title != '':
        ax.set_title(title)

    plt.locator_params(nbins=3)
    if hull._non_elemental:
        ax.set_xlabel(r'x in ({d[0]})$_\mathrm{{x}}$({d[1]})$_\mathrm{{1-x}}$'.format(d=chempot_labels))
    else:
        ax.set_xlabel(r'x in {d[0]}$_\mathrm{{x}}${d[1]}$_\mathrm{{1-x}}$'.format(d=chempot_labels))

    ax.grid(False)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(ax.get_xticks())
    ax.set_ylabel('Formation energy (eV/atom)')

    if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
        import os
        fname = plot_fname or ''.join(hull.species) + '_hull'
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                fname_tmp = fname
                ind = 0
                while os.path.isfile('{}.{}'.format(fname_tmp, ext)):
                    ind += 1
                    fname_tmp = fname + str(ind)

                fname = fname_tmp
                plt.savefig('{}.{}'.format(fname, ext),
                            bbox_inches='tight', transparent=True)
                print('Wrote {}.{}'.format(fname, ext))

    if show:
        plt.show()

    return ax


@plotting_function
def plot_ensemble_hull(hull, data_key,
                       ax=None,
                       formation_energy_key='formation_enthalpy_per_atom',
                       plot_points=False,
                       plot_hull_points=True,
                       alpha_scale=0.25,
                       plot_hulls=True,
                       voltages=False,
                       show=True,
                       plot_fname=None,
                       **kwargs):
    """ Plot and generate an ensemble of hulls. If axis not requested,
    a histogram of frequency of a particular concentration appearing on
    the convex hull is also generated as a second axis.

    Parameters:
        hull (QueryConvexHull): hull object created with a cursor that
            contains the data key for all entries in hull.cursor.
        data_key (str): the key under which ensemble data is stored.

    Keyword arguments:
        ax (matplotlib.axes.Axes): matplotlib axis object on which to plot.
        formation_energy_key (str): the key under which formation energies have been stored.
        alpha_scale (float): value by which to scale transparency of hulls.
        plot_points (bool): whether to plot the hull points for each hull in the ensemble.
        plot_hulls (bool): whether to plot the hull tie-lines for each hull in the ensemble.
        voltages (bool): compute average voltage and heatmaps.

    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

    n_hulls = len(hull.phase_diagrams)
    plot_2d_hull(hull, ax=ax, plot_points=False, plot_hull_points=True, show=False, **kwargs)
    min_ef = 0
    colours_list = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    alpha = min([1, max([1/(alpha_scale*n_hulls), 0.01])])
    for ind, _ in enumerate(hull.phase_diagrams):
        hull_cursor = [doc for doc in hull.cursor if doc[data_key]['hull_distance'][ind] <= 0.0 + EPS]
        min_ef = np.min([doc[data_key][formation_energy_key][ind] for doc in hull_cursor] + [min_ef])
        if plot_hulls:
            ax.plot([doc['concentration'][0] for doc in hull_cursor],
                    [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
                    alpha=alpha, c='k', lw=0.5, zorder=0)
        if plot_hull_points:
            ax.scatter([doc['concentration'][0] for doc in hull_cursor],
                       [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
                       alpha=alpha, marker='o', c=colours_list[1], lw=0, zorder=0)
        if plot_points:
            ax.scatter([doc['concentration'][0] for doc in hull.cursor],
                       [doc[data_key][formation_energy_key][ind] for doc in hull.cursor],
                       alpha=alpha, marker='o', c='k', s=5, lw=0, zorder=0)

    ax.set_ylim(1.1*min_ef)

    if hull.savefig or any(kwargs.get(ext) for ext in SAVE_EXTS):
        if plot_fname is not None:
            fname = plot_fname
        else:
            fname = ''.join(hull.species) + data_key + '_hull'
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                plt.savefig('{}.{}'.format(fname, ext),
                            bbox_inches='tight', transparent=True)
                print('Wrote {}.{}'.format(fname, ext))


@plotting_function
def plot_temperature_hull(
    hull,
    cmap='plasma',
    cmap_limits=(0.2, 0.8),
    ax=None,
    formation_energy_key='formation_free_energy_per_atom',
    plot_points=False,
    plot_hull_points=True,
    alpha_scale=0.25,
    lw_scale=1,
    plot_hulls=True,
    voltages=False,
    show=True,
    plot_fname=None,
    **kwargs
):
    """ Plot and generate an ensemble of hulls. If axis not requested,
    a histogram of frequency of a particular concentration appearing on
    the convex hull is also generated as a second axis.

    Parameters:
        hull (QueryConvexHull): hull object created with a cursor that
            contains the data key for all entries in hull.cursor.

    Keyword arguments:
        ax (matplotlib.axes.Axes): matplotlib axis object on which to plot.
        formation_energy_key (str): the key under which formation energies have been stored.
        alpha_scale (float): value by which to scale transparency of hulls.
        plot_points (bool): whether to plot the hull points for each hull in the ensemble.
        plot_hulls (bool): whether to plot the hull tie-lines for each hull in the ensemble.
        voltages (bool): compute average voltage and heatmaps.

    """
    import matplotlib.pyplot as plt
    data_key = hull.data_key

    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

    n_hulls = len(hull.phase_diagrams)
    colours = plt.cm.get_cmap(cmap)(np.linspace(*cmap_limits, n_hulls))
    min_ef = 0
    alpha = alpha_scale

    # hack the energy key so that labels work
    _cached_key = hull.energy_key
    hull.energy_key = hull.chempot_energy_key
    # set up initial plot without plotting tie line
    ax = plot_2d_hull(
        hull, ax=ax,
        plot_tie_line=False,
        plot_points=False,
        plot_hull_points=False,
        show=False,
        **kwargs
    )
    hull.energy_key = _cached_key

    ax.plot([doc['concentration'][0] for doc in hull.hull_cursor],
            [doc['formation_' + hull.chempot_energy_key] for doc in hull.hull_cursor],
            marker='o',
            alpha=1, c='k', lw=2*lw_scale, ls='--', label='Static', zorder=1e5)

    ind = 0
    hull_cursor = [doc for doc in hull.cursor if doc[data_key]['hull_distance'][ind] <= 0.0 + EPS]
    ax.plot([doc['concentration'][0] for doc in hull_cursor],
            [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
            marker='o', alpha=1, c=colours[ind], lw=2*lw_scale,
            markeredgewidth=1.5, markeredgecolor='k',
            ls='--', zorder=1e5, label='Static + ZPE')

    # plot remaining temperatures
    for ind, _ in enumerate(hull.phase_diagrams[1:]):
        hull_cursor = [doc for doc in hull.cursor if doc[data_key]['hull_distance'][ind] <= 0.0 + EPS]
        min_ef = np.min([doc[data_key][formation_energy_key][ind] for doc in hull_cursor] + [min_ef])
        if plot_hulls:
            ax.plot([doc['concentration'][0] for doc in hull_cursor],
                    [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
                    alpha=alpha, c=colours[ind], lw=1*lw_scale, zorder=0)
        if plot_hull_points:
            ax.scatter([doc['concentration'][0] for doc in hull_cursor],
                       [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
                       alpha=alpha, marker='o', c=colours[ind], lw=0.5,
                       edgecolor='k',
                       zorder=1e4)
        if plot_points:
            ax.scatter([doc['concentration'][0] for doc in hull.cursor],
                       [doc[data_key][formation_energy_key][ind] for doc in hull.cursor],
                       alpha=alpha, marker='o', c=colours[ind], lw=0, zorder=1e-3)

    ax.set_ylim(1.1*min_ef)

    if hull.savefig or any(kwargs.get(ext) for ext in SAVE_EXTS):
        if plot_fname is not None:
            fname = plot_fname
        else:
            fname = ''.join(hull.species) + data_key + '_hull'
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                plt.savefig('{}.{}'.format(fname, ext),
                            bbox_inches='tight', transparent=True)
                print('Wrote {}.{}'.format(fname, ext))

    ax.legend()

    return ax


@plotting_function
def plot_ternary_hull(hull, axis=None, show=True, plot_points=True, hull_cutoff=None,
                      label_cutoff=None, label_corners=True, expecting_cbar=True, labels=None, plot_fname=None, **kwargs):
    """ Plot calculated ternary hull as a 2D projection.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        axis (matplotlib.axes.Axes): matplotlib axis object on which to plot.
        show (bool): whether or not to show plot in X window.
        plot_points (bool): whether or not to plot each structure as a point.
        label_cutoff (float/:obj:`tuple` of :obj:`float`): draw labels less than or
            between these distances form the hull, also read from hull.args.
        expecting_cbar (bool): whether or not to space out the plot to preserve
            aspect ratio if a colourbar is present.
        labels (bool): whether or not to label on-hull structures
        label_corners (bool): whether or not to put axis labels on corners or edges.
        png/pdf/svg (bool): whether or not to write the plot to a file.
        plot_fname (str): filename to write plot to.

    Returns:
        matplotlib.axes.Axes: matplotlib axis with plot.

    """
    import ternary
    import matplotlib.pyplot as plt
    import matplotlib.colors as colours
    from matador.utils.chem_utils import get_generic_grav_capacity

    plt.rcParams['axes.linewidth'] = 0
    plt.rcParams['xtick.major.size'] = 0
    plt.rcParams['ytick.major.size'] = 0
    plt.rcParams['xtick.minor.size'] = 0
    plt.rcParams['ytick.minor.size'] = 0

    if labels is None:
        labels = hull.args.get('labels')
    if label_cutoff is None:
        label_cutoff = hull.args.get('label_cutoff')
        if label_cutoff is None:
            label_cutoff = 0
    else:
        labels = True

    if hull_cutoff is None and hull.hull_cutoff is None:
        hull_cutoff = 0
    else:
        hull_cutoff = hull.hull_cutoff

    print('Plotting ternary hull...')
    if hull.args.get('capmap') or hull.args.get('efmap'):
        scale = 100
    elif hull.args.get('sampmap'):
        scale = 20
    else:
        scale = 1
    fontsize = plt.rcParams['font.size']

    if axis is not None:
        fig, ax = ternary.figure(scale=scale, ax=axis)
    else:
        fig, ax = ternary.figure(scale=scale)
    if hull.args.get('capmap') or hull.args.get('efmap') or hull.args.get('sampmap'):
        fig.set_size_inches(8, 5)
    elif not expecting_cbar:
        fig.set_size_inches(5, 5)
    else:
        fig.set_size_inches(6.67, 5)

    ax.boundary(linewidth=2.0, zorder=99)
    ax.clear_matplotlib_ticks()

    chempot_labels = [get_formula_from_stoich(get_stoich_from_formula(species, sort=False), sort=False, tex=True) for species in hull.species]

    ax.gridlines(color='black', multiple=scale * 0.1, linewidth=0.5)
    ticks = [float(val) for val in np.linspace(0, 1, 6)]
    if label_corners:
        # remove 0 and 1 ticks when labelling corners
        ticks = ticks[1:-1]
        ax.left_corner_label(chempot_labels[2], fontsize=fontsize + 4)
        ax.right_corner_label(chempot_labels[0], fontsize=fontsize + 4)
        ax.top_corner_label(chempot_labels[1], fontsize=fontsize + 4, offset=0.16)
    else:
        ax.left_axis_label(chempot_labels[2], fontsize=fontsize + 2, offset=0.12)
        ax.right_axis_label(chempot_labels[1], fontsize=fontsize + 2, offset=0.12)
        ax.bottom_axis_label(chempot_labels[0], fontsize=fontsize + 2, offset=0.08)
        ax.set_title('-'.join(['{}'.format(label) for label in chempot_labels]), fontsize=fontsize + 2, y=1.02)

    ax.ticks(axis='lbr', linewidth=1, offset=0.02, fontsize=fontsize-2,
             locations=(scale * np.asarray(ticks)).tolist(),
             ticks=ticks, tick_formats='%.1f')

    concs = np.zeros((len(hull.structures), 3))
    concs[:, :-1] = hull.structures[:, :-1]
    for i in range(len(concs)):
        # set third triangular coordinate
        concs[i, -1] = 1 - concs[i, 0] - concs[i, 1]

    stable = concs[np.where(hull.hull_dist <= 0 + EPS)]

    # sort by hull distances so things are plotting the right order
    concs = concs[np.argsort(hull.hull_dist)].tolist()
    hull_dist = np.sort(hull.hull_dist)

    filtered_concs = []
    filtered_hull_dists = []
    for ind, conc in enumerate(concs):
        if conc not in filtered_concs:
            if hull_dist[ind] <= hull.hull_cutoff or (hull.hull_cutoff == 0 and hull_dist[ind] < 0.1):
                filtered_concs.append(conc)
                filtered_hull_dists.append(hull_dist[ind])
    if hull.args.get('debug'):
        print('Trying to plot {} points...'.format(len(filtered_concs)))

    concs = np.asarray(filtered_concs)
    hull_dist = np.asarray(filtered_hull_dists)

    min_cut = 0.0
    max_cut = 0.2

    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    hull.default_cmap_list = get_linear_cmap(hull.colours[1:4], list_only=True)
    hull.default_cmap = get_linear_cmap(hull.colours[1:4], list_only=False)
    n_colours = len(hull.default_cmap_list)
    colours_hull = hull.default_cmap_list

    cmap = hull.default_cmap
    cmap_full = plt.cm.get_cmap('Pastel2')
    pastel_cmap = colours.LinearSegmentedColormap.from_list('Pastel2', cmap_full.colors)

    for plane in hull.convex_hull.planes:
        plane.append(plane[0])
        plane = np.asarray(plane)
        ax.plot(scale * plane, c=hull.colours[0], lw=1.5, alpha=1, zorder=98)

    if hull.args.get('pathways'):
        for phase in stable:
            if phase[0] == 0 and phase[1] != 0 and phase[2] != 0:
                ax.plot([scale * phase, [scale, 0, 0]], c='r', alpha=0.2, lw=6, zorder=99)

    # add points
    if plot_points:
        colours_list = []
        colour_metric = hull_dist
        for i, _ in enumerate(colour_metric):
            if colour_metric[i] >= max_cut:
                colours_list.append(n_colours - 1)
            elif colour_metric[i] <= min_cut:
                colours_list.append(0)
            else:
                colours_list.append(int((n_colours - 1) * (colour_metric[i] / max_cut)))
        colours_list = np.asarray(colours_list)
        ax.scatter(scale*stable, marker='o', color=hull.colours[1], edgecolors='black', zorder=9999999,
                   s=150, lw=1.5)
        ax.scatter(scale*concs, colormap=cmap, colorbar=True, cbarlabel='Distance from hull (eV/atom)',
                   c=colour_metric, vmax=max_cut, vmin=min_cut, zorder=1000, s=40, alpha=0)
        for i, _ in enumerate(concs):
            ax.scatter(
                scale * concs[i].reshape(1, 3),
                color=colours_hull[colours_list[i]],
                marker='o',
                zorder=10000 - colours_list[i],
                s=70 * (1 - float(colours_list[i]) / n_colours) + 15,
                lw=1,
                edgecolors='black'
            )

    # add colourmaps
    if hull.args.get('capmap'):
        capacities = dict()
        from ternary.helpers import simplex_iterator
        for (i, j, k) in simplex_iterator(scale):
            capacities[(i, j, k)] = get_generic_grav_capacity([
                float(i) / scale, float(j) / scale, float(scale - i - j) / scale
            ], hull.species)
        ax.heatmap(capacities, style="hexagonal", cbarlabel='Gravimetric capacity (mAh/g)',
                   vmin=0, vmax=3000, cmap=pastel_cmap)
    elif hull.args.get('efmap'):
        energies = dict()
        fake_structures = []
        from ternary.helpers import simplex_iterator
        for (i, j, k) in simplex_iterator(scale):
            fake_structures.append([float(i) / scale, float(j) / scale, 0.0])
        fake_structures = np.asarray(fake_structures)
        plane_energies = hull.get_hull_distances(fake_structures, precompute=False)
        ind = 0
        for (i, j, k) in simplex_iterator(scale):
            energies[(i, j, k)] = -1 * plane_energies[ind]
            ind += 1
        if isinstance(hull.args.get('efmap'), str):
            efmap = hull.args.get('efmap')
        else:
            efmap = 'BuPu_r'
        ax.heatmap(energies, style="hexagonal", cbarlabel='Formation energy (eV/atom)', vmax=0, cmap=efmap)
    elif hull.args.get('sampmap'):
        sampling = dict()
        from ternary.helpers import simplex_iterator
        eps = 1.0 / float(scale)
        for (i, j, k) in simplex_iterator(scale):
            sampling[(i, j, k)] = np.size(np.where((concs[:, 0] <= float(i)/scale + eps) *
                                                   (concs[:, 0] >= float(i)/scale - eps) *
                                                   (concs[:, 1] <= float(j)/scale + eps) *
                                                   (concs[:, 1] >= float(j)/scale - eps) *
                                                   (concs[:, 2] <= float(k)/scale + eps) *
                                                   (concs[:, 2] >= float(k)/scale - eps)))
        ax.heatmap(sampling, style="hexagonal", cbarlabel='Number of structures', cmap='afmhot')

    # add labels
    if labels:
        label_cursor = get_hull_labels(hull, label_cutoff=label_cutoff)
        if len(label_cursor) == 1:
            label_coords = [[0.25, 0.5]]
        else:
            label_coords = [[0.1+(val-0.5)*0.3, val] for val in np.linspace(0.5, 0.8, int(round(len(label_cursor)/2.)+1))]
            label_coords += [[0.9-(val-0.5)*0.3, val+0.2] for val in np.linspace(0.5, 0.8, int(round(len(label_cursor)/2.)))]
        from matador.utils.hull_utils import barycentric2cart
        for ind, doc in enumerate(label_cursor):
            conc = np.asarray(doc['concentration'] + [1 - sum(doc['concentration'])])
            formula = get_formula_from_stoich(doc['stoichiometry'], sort=False, tex=True, latex_sub_style=r'\mathregular', elements=hull.species)
            arrowprops = dict(arrowstyle="-|>", color='k', lw=2, alpha=0.5, zorder=1, shrinkA=2, shrinkB=4)
            cart = barycentric2cart([doc['concentration'] + [0]])[0][:2]
            min_dist = 1e20
            closest_label = 0
            for coord_ind, coord in enumerate(label_coords):
                dist = np.sqrt((cart[0] - coord[0])**2 + (cart[1] - coord[1])**2)
                if dist < min_dist:
                    min_dist = dist
                    closest_label = coord_ind
            ax.annotate(formula, scale*conc,
                        textcoords='data',
                        xytext=[scale*val for val in label_coords[closest_label]],
                        ha='right',
                        va='bottom',
                        arrowprops=arrowprops)
            del label_coords[closest_label]

    plt.tight_layout(w_pad=0.2)

    if hull.savefig:
        fname = plot_fname or ''.join(hull.species) + '_hull'
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                plt.savefig('{}.{}'.format(fname, ext),
                            bbox_inches='tight', transparent=True)
                print('Wrote {}.{}'.format(fname, ext))
    elif show:
        print('Showing plot...')
        plt.show()

    return ax


def _scatter_plot_by_source(hull, ax, scale, kwargs,
                            sources=None, source_labels=None, plot_hull_points=True):
    """ Add scatter points to the hull depending on the guessed
    provenance of a structure.

    """
    from matador.utils.cursor_utils import get_guess_doc_provenance
    if sources is None:
        sources = ['AIRSS', 'GA', 'OQMD', 'SWAPS', 'ICSD', 'DOI', 'SM', 'MP', 'PF', 'Other']
    if source_labels is None:
        source_labels = sources
    else:
        assert len(source_labels) == len(sources)

    # hack: double length of hull colours
    hull.colours.extend(hull.colours)

    colour_choices = {source: hull.colours[ind + 1] for ind, source in enumerate(sources)}
    colours = []
    concs = []
    energies = []
    zorders = []
    for doc in hull.cursor:
        source = get_guess_doc_provenance(doc['source'])
        if source not in sources:
            # use grey for undesired sources
            colours.append(hull.colours[-2])
            source = 'Other'
            if 'Other' not in sources:
                sources.append('Other')
                source_labels.append('Other')
                colour_choices['Other'] = hull.colours[-2]
        else:
            colours.append(source)
        zorders.append(sources.index(source))
        concs.append(doc['concentration'])
        energies.append(doc['formation_{}'.format(hull.energy_key)])

    sources_present = set(colours)
    sources_present = [source for source in sources if source in sources_present]
    colour_choices = {source: hull.colours[ind + 1] for ind, source in enumerate(sources_present)}
    colours = [colour_choices[src] for src in colours]

    alpha = kwargs.get('alpha')
    if alpha is None:
        alpha = 0.2

    for ind, conc in enumerate(concs):
        if hull.cursor[ind]['hull_distance'] <= 0 + 1e-9 and not plot_hull_points:
            ax.scatter(conc, energies[ind],
                       c=colours[ind], alpha=1, s=scale*40, edgecolor='k',
                       zorder=zorders[ind]+1e5, lw=1.5)
        else:
            ax.scatter(conc, energies[ind],
                       c=colours[ind], alpha=alpha, s=scale*20,
                       zorder=zorders[ind]+100)

    for ind, source in enumerate(sources_present):
        ax.scatter(1e10, 1e10, c=colour_choices[source], label=sources_present[ind], alpha=alpha, lw=1)

    legend = ax.legend(loc='lower right', facecolor='w', frameon=True, fancybox=False, shadow=False, ncol=2)
    legend.set_zorder(1e20)

    return ax
