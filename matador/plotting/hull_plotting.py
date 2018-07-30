# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot convex hulls and phase
diagrams generally.

"""


import numpy as np
from matador.utils.chem_utils import get_formula_from_stoich
from matador.plotting.plotting import plotting_function, get_linear_cmap


def get_hull_labels(hull, label_cutoff=None, num_species=2):
    """ Return list of structures to labels on phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): phase diagram to plot.

    Keyword arguments:
        num_species (int): structures containing this number of species
            will be labelled.

    Returns:
        label_cursor (list(dict)): list of matador documents to label.

    """
    eps = 1e-9
    if label_cutoff is None:
        label_cutoff = 0.0
    if isinstance(label_cutoff, list) and len(label_cutoff) == 2:
        label_cutoff = sorted(label_cutoff)
        # first, only apply upper limit as we need to filter by stoich aftewards
        label_cursor = [doc for doc in hull.hull_cursor if doc['hull_distance'] <= label_cutoff[1]]
    else:
        if isinstance(label_cutoff, list):
            assert len(label_cutoff) == 1, 'Incorrect number of label_cutoff values passed, should be 1 or 2.'
            label_cutoff = label_cutoff[0]
        label_cursor = [doc for doc in hull.hull_cursor if doc['hull_distance'] <= label_cutoff + eps]

    num_labels = len(set([get_formula_from_stoich(doc['stoichiometry']) for doc in label_cursor]))
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
    label_cursor = [doc for doc in label_cursor if len(doc['stoichiometry']) == num_species]

    return label_cursor


@plotting_function
def plot_2d_hull(hull, ax=None, show=False, plot_points=True,
                 plot_hull_points=True, labels=None, label_cutoff=None, colour_by_source=False,
                 sources=None, source_labels=None, title=True,
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

    Returns:
        matplotlib.axes.Axes: matplotlib axis with plot.

    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as colours

    if ax is None:
        fig = plt.figure(facecolor=None, figsize=(8, 6))
        ax = fig.add_subplot(111)

    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    hull.default_cmap_list = get_linear_cmap(hull.colours[1:4], list_only=True)
    hull.default_cmap = get_linear_cmap(hull.colours[1:4], list_only=False)

    if labels is None:
        labels = hull.args.get('labels', False)
    if label_cutoff is None:
        label_cutoff = hull.args.get('label_cutoff')

    scale = 1
    scatter = []
    x_elem = [hull.elements[0]]
    one_minus_x_elem = list(hull.elements[1:])
    tie_line = hull.structure_slice[hull.hull.vertices]

    # plot hull structures
    if plot_hull_points:
        ax.scatter(tie_line[:, 0], tie_line[:, 1],
                   c=hull.colours[1],
                   marker='o', zorder=99999, edgecolor='k',
                   s=scale*40, lw=1.5)
    ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1],
            c=hull.colours[0], zorder=1)

    if hull.hull_cutoff > 0:
        ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1] + hull.hull_cutoff,
                '--', c=hull.colours[1], alpha=0.5, zorder=1, label='')

    # annotate hull structures
    if labels or label_cutoff is not None:
        label_cursor = get_hull_labels(hull, num_species=2, label_cutoff=label_cutoff)
        for ind, doc in enumerate(label_cursor):
            arrowprops = dict(arrowstyle="-|>", lw=2, alpha=1, zorder=1, shrinkA=2, shrinkB=4)
            min_comp = tie_line[np.argmin(tie_line[:, 1]), 0]
            e_f = label_cursor[ind]['formation_' + str(hull._energy_key)]
            conc = label_cursor[ind]['concentration'][0]
            if conc < min_comp:
                position = (0.8 * conc, 1.15 * (e_f - 0.05))
            elif label_cursor[ind]['concentration'][0] == min_comp:
                position = (conc, 1.15 * (e_f - 0.05))
            else:
                position = (min(1.1 * conc + 0.15, 0.95),1.15 * (e_f - 0.05))
            # if (ind + 2) < np.argmin(tie_line[:, 1]):
                # position = (0.8 * tie_line[ind + 2, 0], 1.15 * (tie_line[ind + 2, 1]) - 0.05)
            # elif (ind + 2) == np.argmin(tie_line[:, 1]):
                # position = (tie_line[ind + 2, 0], 1.15 * (tie_line[ind + 2, 1]) - 0.05)
            # else:
                # position = (min(1.1 * tie_line[ind + 2, 0] + 0.15, 0.95), 1.15 * (tie_line[ind + 2, 1]) - 0.05)
            ax.annotate(get_formula_from_stoich(doc['stoichiometry'],
                                                elements=hull.elements,
                                                latex_sub_style=r'\mathregular',
                                                tex=True),
                        xy=(conc, e_f),
                        xytext=position,
                        textcoords='data',
                        ha='right',
                        va='bottom',
                        arrowprops=arrowprops,
                        zorder=1)

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

                cbar = plt.colorbar(scatter, aspect=30, pad=0.02,
                                    ticks=[0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                cbar.ax.tick_params(length=0)
                cbar.ax.set_yticklabels([0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                cbar.ax.yaxis.set_ticks_position('right')
                cbar.ax.set_frame_on(False)
                cbar.outline.set_visible(False)
                cbar.set_label('Distance from hull (eV)')

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
        from matador.utils.cursor_utils import get_guess_doc_provenance
        if sources is None:
            sources = ['AIRSS', 'GA', 'OQMD', 'SWAPS', 'ICSD']
        if source_labels is None:
            source_labels = sources
        else:
            assert len(source_labels) == len(sources)

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
                if 'Other' not in sources:
                    sources.append('Other')
                    labels.append('Other')
                    colour_choices['Other'] = hull.colours[-2]
            else:
                colours.append(colour_choices[source])
            zorders.append(sources.index(source))
            concs.append(doc['concentration'])
            energies.append(doc['formation_enthalpy_per_atom'])

        alpha = kwargs.get('alpha')
        if alpha is None:
            alpha = 0.2

        for ind, conc in enumerate(concs):
            if hull.cursor[ind]['hull_distance'] <= 0 + 1e-9 and not plot_hull_points:
                ax.scatter(conc, energies[ind],
                           c=colours[ind], alpha=alpha, s=scale*40,
                           zorder=zorders[ind]+1e5, lw=1.5)
            else:
                ax.scatter(conc, energies[ind],
                           c=colours[ind], alpha=alpha, s=scale*20,
                           zorder=zorders[ind]+100)

        for ind, source in enumerate(sources):
            ax.scatter(1e10, 1e10, c=colour_choices[source], label=source_labels[ind], alpha=alpha, lw=1)

        legend = ax.legend(loc=9, facecolor='w', frameon=True, fancybox=False, shadow=False)
        legend.set_zorder(1e20)

    eform_limits = (np.min(hull.structures[:, 1]), np.max(hull.structures[:, 1]))
    lims = (-0.1 if eform_limits[0] >= 0 else 1.4*eform_limits[0],
            eform_limits[1] if eform_limits[0] >= 0 else 0.1)
    ax.set_ylim(lims)

    if isinstance(title, bool) and title:
        if len(one_minus_x_elem) == 1:
            ax.set_title(x_elem[0] + r'$_\mathrm{x}$' + one_minus_x_elem[0] + r'$_\mathrm{1-x}$')
        if hull._non_binary:
            ax.set_title(r'{d[0]}$_\mathrm{{x}}$({d[1]})$_\mathrm{{1-x}}$'.format(d=hull.chempot_search))
    elif title:
        ax.set_title(title)

    plt.locator_params(nbins=3)
    ax.set_xlabel(r'x in {}$_\mathrm{{x}}${}$_\mathrm{{1-x}}$'.format(x_elem[0], one_minus_x_elem[0]))
    ax.grid(False)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.25, 0.33, 0.5, 0.66, 0.75, 1])
    ax.set_xticklabels(ax.get_xticks())
    # ax.set_yticks(np.arange(, np.min(hull.structure_slice[hull.hull.vertices, 1]) - 0.15, -0.1))
    # ax.set_yticklabels(['{:.1f}'.format(val) for val in ax.get_yticks()])
    ax.set_ylabel('Formation energy (eV/atom)')

    if hull.savefig:
        fname = ''.join(hull.elements) + '_hull'
        exts = ['pdf', 'svg', 'png']
        for ext in exts:
            if hull.args.get(ext):
                plt.savefig('{}.{}'.format(fname, ext),
                            dpi=500, bbox_inches='tight', transparent=True)
                print('Wrote {}.{}'.format(fname, ext))
    elif show:
        plt.show()

    return ax


@plotting_function
def plot_ternary_hull(hull, axis=None, show=False, plot_points=True, hull_cutoff=None,
                      label_cutoff=None, expecting_cbar=True, labels=None, **kwargs):
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
    ax.gridlines(color='black', multiple=scale * 0.1, linewidth=0.5)

    ax.clear_matplotlib_ticks()
    ticks = [float(val) for val in np.linspace(0.0, 1.0, 6)]
    ax.ticks(axis='lbr', linewidth=1, multiple=scale*0.2, offset=0.02, fontsize=fontsize-2,
             ticks=ticks, tick_formats='%.1f')

    ax.set_title('-'.join(hull.elements), fontsize=fontsize + 2, y=1.02)
    ax.left_axis_label(hull.elements[2], fontsize=fontsize + 2)
    ax.right_axis_label(hull.elements[1], fontsize=fontsize + 2)
    ax.bottom_axis_label(hull.elements[0], fontsize=fontsize + 2)

    concs = np.zeros((len(hull.structures), 3))

    concs[:, :-1] = hull.structures[:, :-1]
    for i in range(len(concs)):
        # set third triangular coordinate
        concs[i, -1] = 1 - concs[i, 0] - concs[i, 1]

    stable = np.asarray([concs[ind] for ind in hull.hull.vertices])

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

    for plane in hull.hull.planes:
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
                   c=hull_dist, vmax=max_cut, vmin=min_cut, zorder=1000, s=40, alpha=0)
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
            ], hull.elements)
        ax.heatmap(capacities, style="hexagonal", cbarlabel='Gravimetric capacity (maH/g)',
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
        ax.heatmap(energies, style="hexagonal", cbarlabel='Formation energy (eV/atom)', vmax=0, cmap='bone')
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
        label_cursor = get_hull_labels(hull, label_cutoff=label_cutoff, num_species=3)
        if len(label_cursor) == 1:
            label_coords = [[0.25, 0.5]]
        else:
            label_coords = [[0.1+(val-0.5)*0.3, val] for val in np.linspace(0.5, 0.8, int(round(len(label_cursor)/2.)+1))]
            label_coords += [[0.9-(val-0.5)*0.3, val+0.2] for val in np.linspace(0.5, 0.8, int(round(len(label_cursor)/2.)))]
        from matador.utils.hull_utils import barycentric2cart
        for ind, doc in enumerate(label_cursor):
            conc = np.asarray(doc['concentration'] + [1 - sum(doc['concentration'])])
            formula = get_formula_from_stoich(doc['stoichiometry'], tex=True, elements=hull.elements, latex_sub_style=r'\mathregular')
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
        if hull.args.get('png'):
            plt.savefig(''.join(hull.elements) + '_hull.png', dpi=400, transparent=True, bbox_inches='tight')
        if hull.args.get('svg'):
            plt.savefig(''.join(hull.elements) + '_hull.svg', dpi=400, transparent=True, bbox_inches='tight')
        if hull.args.get('pdf'):
            plt.savefig(''.join(hull.elements) + '_hull.pdf', dpi=400, transparent=True, bbox_inches='tight')
    elif show:
        ax.show()

    return ax
