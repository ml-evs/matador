# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot convex hulls and phase
diagrams generally.

"""


from collections import defaultdict
import numpy as np
from matador.utils.chem_utils import get_stoich_from_formula, get_formula_from_stoich
from matador.utils.viz_utils import (
    get_element_colours,
    colour_from_ternary_concentration,
)
from matador.plotting.plotting import plotting_function, get_linear_cmap, SAVE_EXTS

EPS = 1e-12

__all__ = [
    "plot_2d_hull",
    "plot_ternary_hull",
    "plot_ensemble_hull",
    "plot_temperature_hull",
]


def _get_hull_labels(hull, label_cutoff=None, num_species=None, exclude_edges=True):
    """Return list of structures to labels on phase diagram.

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
    if label_cutoff is None:
        label_cutoff = 0.0
    if isinstance(label_cutoff, list) and len(label_cutoff) == 2:
        label_cutoff = sorted(label_cutoff)
        # first, only apply upper limit as we need to filter by stoich aftewards
        label_cursor = [
            doc for doc in hull.cursor if doc["hull_distance"] <= label_cutoff[1]
        ]
    else:
        if isinstance(label_cutoff, list):
            assert (
                len(label_cutoff) == 1
            ), "Incorrect number of label_cutoff values passed, should be 1 or 2."
            label_cutoff = label_cutoff[0]
        label_cursor = [
            doc for doc in hull.cursor if doc["hull_distance"] <= label_cutoff + EPS
        ]

    label_cursor = sorted(label_cursor, key=lambda _: _["hull_distance"], reverse=False)

    num_labels = len(
        {get_formula_from_stoich(doc["stoichiometry"]) for doc in label_cursor}
    )
    if num_labels < len(label_cursor):
        tmp_cursor = []
        for doc in label_cursor:
            if doc["stoichiometry"] not in [
                _doc["stoichiometry"] for _doc in tmp_cursor
            ]:
                tmp_cursor.append(doc)
            else:
                label_cursor = tmp_cursor
    if isinstance(label_cutoff, list) and len(label_cutoff) == 2:
        # now apply lower limit
        label_cursor = [
            doc
            for doc in label_cursor
            if label_cutoff[0] <= doc["hull_distance"] <= label_cutoff[1]
        ]
    # remove chemical potentials and unwanted e.g. binaries
    if num_species is not None:
        label_cursor = [
            doc for doc in label_cursor if len(doc["stoichiometry"]) == num_species
        ]
    if exclude_edges:
        label_cursor = [
            doc
            for doc in label_cursor
            if (all(doc["concentration"]) > 0 and sum(doc["concentration"]) <= 1 - EPS)
        ]

    label_cursor = sorted(label_cursor, key=lambda doc: doc["concentration"])

    return label_cursor


@plotting_function
def plot_2d_hull(
    hull,
    ax=None,
    show=True,
    plot_points=True,
    plot_tie_line=True,
    plot_hull_points=True,
    labels=None,
    label_cutoff=None,
    colour_by_source=False,
    sources=None,
    hull_label=None,
    source_labels=None,
    title=True,
    plot_fname=None,
    show_cbar=True,
    colour_by_composition=False,
    label_offset=(1, 0.1),
    specific_label_offset=None,
    eform_limits=None,
    legend_kwargs=None,
    hull_dist_unit="meV",
    **kwargs,
):
    """Plot calculated hull, returning ax and fig objects for further editing.

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
        specific_label_offset (dict): A dictionary keyed by chemical formula strings
            containing 2-length float tuples of label offsets. A value of `None` for
            a particular formula will exclude that label from the plot.
        colour_by_source (bool): plot and label points by their sources
        show_cbar: Whether or not to plot the hull distance colourbar.
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
    import matplotlib.colors

    # cache the specific label offset dict
    _specific_label_offset = None
    if specific_label_offset:
        _specific_label_offset = dict(specific_label_offset)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if not hasattr(hull, "colours"):
        hull.colours = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    hull.default_cmap_list = get_linear_cmap(hull.colours[1:4], list_only=True)
    hull.default_cmap = get_linear_cmap(hull.colours[1:4], list_only=False)

    if colour_by_composition:
        try:
            element_colours = [get_element_colours()[s] for s in hull.species]
        except Exception:
            raise RuntimeError(
                f"Cannot `colour_by_composition`: no colour found for species {hull.species}."
            )
        conc_cmap = get_linear_cmap(element_colours, list_only=False)

    if labels is None:
        labels = hull.args.get("labels", False)
    if label_cutoff is None:
        label_cutoff = hull.args.get("label_cutoff")

    scale = 1
    scatter = []
    chempot_labels = [
        get_formula_from_stoich(get_stoich_from_formula(species, sort=False), tex=True)
        for species in hull.species
    ]
    tie_line = hull.convex_hull.points[hull.convex_hull.vertices]

    # plot hull structures
    if plot_hull_points:
        if colour_by_composition:
            hull_point_scale = 75
            edgewidth = 1.5
            edgecolor = kwargs.get("markeredgecolor", "black")
        else:
            hull_point_scale = 40
            edgewidth = 1.5
            edgecolor = kwargs.get("markeredgecolor", "black")

        if colour_by_composition:
            point_colours = tie_line[:, 0]
            cmap = conc_cmap
        else:
            point_colours = hull.colours[1]
            cmap = None
        ax.scatter(
            tie_line[:, 0],
            tie_line[:, 1],
            c=point_colours,
            cmap=cmap,
            marker="o",
            zorder=99999,
            edgecolor=edgecolor,
            s=scale * hull_point_scale,
            lw=edgewidth,
        )
        if plot_tie_line:
            ax.plot(
                np.sort(tie_line[:, 0]),
                tie_line[np.argsort(tie_line[:, 0]), 1],
                c=hull.colours[0],
                zorder=1,
                label=hull_label,
                marker="o",
                markerfacecolor=hull.colours[0],
                markeredgecolor=edgecolor,
                markeredgewidth=edgewidth,
                markersize=np.sqrt(scale * hull_point_scale),
            )
    if plot_tie_line:
        ax.plot(
            np.sort(tie_line[:, 0]),
            tie_line[np.argsort(tie_line[:, 0]), 1],
            c=hull.colours[0],
            zorder=1,
            label=hull_label,
            markersize=0,
        )

    if hull.hull_cutoff > 0:
        ax.plot(
            np.sort(tie_line[:, 0]),
            tie_line[np.argsort(tie_line[:, 0]), 1] + hull.hull_cutoff,
            "--",
            c=hull.colours[1],
            alpha=0.5,
            zorder=1,
            label="",
        )

    # annotate hull structures
    if labels or label_cutoff is not None:
        label_cursor = _get_hull_labels(hull, num_species=2, label_cutoff=label_cutoff)
        already_labelled = []
        for ind, doc in enumerate(label_cursor):
            formula = get_formula_from_stoich(doc["stoichiometry"], sort=True)
            if formula not in already_labelled:
                arrowprops = dict(
                    arrowstyle="-|>",
                    facecolor="k",
                    edgecolor="k",
                    lw=2,
                    alpha=1,
                    zorder=1,
                    shrinkA=2,
                    shrinkB=4,
                )
                min_comp = tie_line[np.argmin(tie_line[:, 1]), 0]
                e_f = label_cursor[ind]["formation_" + str(hull.energy_key)]
                conc = label_cursor[ind]["concentration"][0]
                if conc < min_comp:
                    position = (0.8 * conc, label_offset[0] * (e_f - label_offset[1]))
                elif label_cursor[ind]["concentration"][0] == min_comp:
                    position = (conc, label_offset[0] * (e_f - label_offset[1]))
                else:
                    position = (
                        min(1.1 * conc + 0.15, 0.95),
                        label_offset[0] * (e_f - label_offset[1]),
                    )
                    position = (
                        min(1.05 * conc + 0.15, 0.95),
                        label_offset[0] * (e_f - label_offset[1]),
                    )
                if _specific_label_offset:
                    plain_formula = get_formula_from_stoich(
                        doc["stoichiometry"], tex=False
                    )
                    if plain_formula in _specific_label_offset:
                        offset = _specific_label_offset.pop(plain_formula)
                        if offset is None:
                            continue
                        position = (position[0] + offset[0], position[1] + offset[1])
                if colour_by_composition:
                    text_colour = conc_cmap(conc)
                else:
                    text_colour = "k"
                ax.annotate(
                    get_formula_from_stoich(
                        doc["stoichiometry"],
                        latex_sub_style=r"\mathregular",
                        tex=True,
                        elements=hull.species,
                        sort=False,
                    ),
                    xy=(conc, e_f),
                    xytext=position,
                    color=text_colour,
                    textcoords="data",
                    ha="right",
                    va="bottom",
                    arrowprops=arrowprops,
                    zorder=1,
                )
                already_labelled.append(formula)

    if _specific_label_offset:
        import warnings

        warnings.warn(f"Found unused requested offsets: {_specific_label_offset}")

    # points for off hull structures; we either colour by source or by energy
    if plot_points and not colour_by_source:

        if hull.hull_cutoff == 0:
            # if no specified hull cutoff, ignore labels and colour by hull distance
            if plot_points:
                concs = hull.structures[np.argsort(hull.hull_dist), 0][::-1]
                energies = (hull.structures[np.argsort(hull.hull_dist), -1][::-1],)
                if colour_by_composition:
                    point_colours = concs
                    norm = None
                    cmap = conc_cmap
                else:
                    point_colours = np.sort(hull.hull_dist)[::-1]
                    norm = matplotlib.colors.LogNorm(0.01, 1, clip=True)
                    cmap = hull.default_cmap
                scatter = ax.scatter(
                    concs,
                    energies,
                    s=scale * 40,
                    c=point_colours,
                    zorder=100,
                    cmap=cmap,
                    norm=norm,
                    rasterized=True,
                )

                if show_cbar:
                    ticks = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64]
                    cbar = plt.colorbar(
                        scatter, aspect=30, pad=0.02, ticks=ticks, extend="both"
                    )
                    if hull_dist_unit.lower() == "mev":
                        ticks = [int(1000 * t) for t in ticks]
                    cbar.ax.set_yticklabels(ticks)
                    cbar.ax.tick_params(which="minor", length=0)
                    unit = "eV/atom"
                    if hull_dist_unit.lower() == "mev":
                        unit = "m" + unit
                    cbar.set_label(f"Distance from hull ({unit})")

        elif hull.hull_cutoff != 0:
            # if specified hull cutoff colour those below
            c = hull.colours[1]
            for ind in range(len(hull.structures)):
                if hull.hull_dist[ind] <= hull.hull_cutoff or hull.hull_cutoff == 0:
                    if plot_points:
                        scatter.append(
                            ax.scatter(
                                hull.structures[ind, 0],
                                hull.structures[ind, 1],
                                s=scale * 40,
                                alpha=0.9,
                                c=c,
                                zorder=300,
                            )
                        )
            if plot_points:
                ax.scatter(
                    hull.structures[1:-1, 0],
                    hull.structures[1:-1, 1],
                    s=scale * 30,
                    lw=0,
                    alpha=0.3,
                    c=hull.colours[-2],
                    edgecolor="k",
                    zorder=10,
                    rasterized=True,
                )

    elif colour_by_source:
        _scatter_plot_by_source(
            hull,
            ax,
            scale,
            kwargs,
            sources=sources,
            source_labels=source_labels,
            plot_hull_points=plot_hull_points,
            legend_kwargs=legend_kwargs,
        )

    if eform_limits is None:
        eform_limits = (np.min(hull.structures[:, 1]), np.max(hull.structures[:, 1]))
        lims = (
            -0.1 if eform_limits[0] >= 0 else 1.25 * eform_limits[0],
            eform_limits[1] if eform_limits[0] >= 0 else 0.05,
        )
    else:
        lims = sorted(eform_limits)
    ax.set_ylim(lims)

    if isinstance(title, bool) and title:
        if hull._non_elemental:
            ax.set_title(
                r"({d[0]})$_\mathrm{{x}}$({d[1]})$_\mathrm{{1-x}}$".format(
                    d=chempot_labels
                )
            )
        else:
            ax.set_title(
                r"{d[0]}$_\mathrm{{x}}${d[1]}$_\mathrm{{1-x}}$".format(d=chempot_labels)
            )
    elif isinstance(title, str) and title != "":
        ax.set_title(title)

    if hull._non_elemental:
        ax.set_xlabel(
            r"x in ({d[0]})$_\mathrm{{x}}$({d[1]})$_\mathrm{{1-x}}$".format(
                d=chempot_labels
            )
        )
    else:
        ax.set_xlabel(
            r"x in {d[0]}$_\mathrm{{x}}${d[1]}$_\mathrm{{1-x}}$".format(
                d=chempot_labels
            )
        )

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylabel("Formation energy (eV/atom)")

    if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
        import os

        fname = plot_fname or "".join(hull.species) + "_hull"
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                fname_tmp = fname
                ind = 0
                while os.path.isfile("{}.{}".format(fname_tmp, ext)):
                    ind += 1
                    fname_tmp = fname + str(ind)

                fname = fname_tmp
                plt.savefig(
                    "{}.{}".format(fname, ext), bbox_inches="tight", transparent=True
                )
                print("Wrote {}.{}".format(fname, ext))

    return ax


@plotting_function
def plot_ensemble_hull(
    hull,
    data_key,
    ax=None,
    formation_energy_key="formation_enthalpy_per_atom",
    plot_points=False,
    plot_hull_points=True,
    alpha_scale=0.25,
    plot_hulls=True,
    **kwargs,
):
    """Plot and generate an ensemble of hulls. If axis not requested,
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
        fig = plt.figure()
        ax = fig.add_subplot(111)

    n_hulls = len(hull.phase_diagrams)
    plot_2d_hull(
        hull, ax=ax, plot_points=False, plot_hull_points=True, show=False, **kwargs
    )
    min_ef = 0
    colours_list = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    alpha = min([1, max([1 / (alpha_scale * n_hulls), 0.01])])
    for ind, _ in enumerate(hull.phase_diagrams):
        hull_cursor = [
            doc
            for doc in hull.cursor
            if doc[data_key]["hull_distance"][ind] <= 0.0 + EPS
        ]
        min_ef = np.min(
            [doc[data_key][formation_energy_key][ind] for doc in hull_cursor] + [min_ef]
        )
        if plot_hulls:
            ax.plot(
                [doc["concentration"][0] for doc in hull_cursor],
                [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
                alpha=alpha,
                c="k",
                lw=0.5,
                zorder=0,
            )
        if plot_hull_points:
            ax.scatter(
                [doc["concentration"][0] for doc in hull_cursor],
                [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
                alpha=alpha,
                marker="o",
                c=colours_list[1],
                lw=0,
                zorder=0,
            )
        if plot_points:
            ax.scatter(
                [doc["concentration"][0] for doc in hull.cursor],
                [doc[data_key][formation_energy_key][ind] for doc in hull.cursor],
                alpha=alpha,
                marker="o",
                c="k",
                s=5,
                lw=0,
                zorder=0,
            )

    ax.set_ylim(1.1 * min_ef)


@plotting_function
def plot_temperature_hull(
    hull,
    cmap="plasma",
    cmap_limits=(0.2, 0.8),
    ax=None,
    formation_energy_key="formation_free_energy_per_atom",
    plot_points=False,
    show_cbar=True,
    plot_hull_points=True,
    alpha_scale=1,
    lw_scale=1,
    plot_hulls=True,
    plot_fname=None,
    **kwargs,
):
    """Plot and generate an ensemble of hulls. If axis not requested,
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
    from matplotlib.colors import rgb2hex

    data_key = hull.data_key

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    n_hulls = len(hull.phase_diagrams)
    colours = [
        rgb2hex(col)
        for col in plt.cm.get_cmap(cmap)(np.linspace(*cmap_limits, n_hulls)).tolist()
    ]

    min_ef = 0
    alpha = alpha_scale

    # hack the energy key so that labels work
    _cached_key = hull.energy_key
    hull.energy_key = hull.chempot_energy_key
    max_temperature = max(hull.temperatures)
    # set up initial plot without plotting tie line
    ax = plot_2d_hull(
        hull,
        ax=ax,
        plot_tie_line=False,
        plot_points=False,
        plot_hull_points=False,
        show=False,
        **kwargs,
    )
    hull.energy_key = _cached_key

    ax.plot(
        [doc["concentration"][0] for doc in hull.hull_cursor],
        [doc["formation_" + hull.chempot_energy_key] for doc in hull.hull_cursor],
        marker="o",
        alpha=1,
        c="k",
        lw=2 * lw_scale,
        ls="--",
        label="Static",
        zorder=1e5,
    )

    ind = 0
    hull_cursor = [
        doc for doc in hull.cursor if doc[data_key]["hull_distance"][ind] <= 0.0 + EPS
    ]
    ax.plot(
        [doc["concentration"][0] for doc in hull_cursor],
        [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
        marker="o",
        c=colours[ind],
        lw=2 * lw_scale,
        markeredgewidth=1.5,
        markeredgecolor="k",
        ls="--",
        zorder=1e5,
        label="Static + ZPE",
    )

    # plot remaining temperatures
    for ind, _ in enumerate(hull.temperatures[1:]):
        hull_cursor = [
            doc
            for doc in hull.cursor
            if doc[data_key]["hull_distance"][ind] <= 0.0 + EPS
        ]
        min_ef = np.min(
            [doc[data_key][formation_energy_key][ind] for doc in hull_cursor] + [min_ef]
        )
        colour_ind = int(n_hulls * hull.temperatures[ind] / max_temperature)
        if plot_hulls:
            ax.plot(
                [doc["concentration"][0] for doc in hull_cursor],
                [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
                alpha=alpha,
                c=colours[colour_ind],
                lw=1 * lw_scale,
                zorder=0,
            )
        if plot_hull_points:
            ax.scatter(
                [doc["concentration"][0] for doc in hull_cursor],
                [doc[data_key][formation_energy_key][ind] for doc in hull_cursor],
                alpha=1,
                marker="o",
                c=colours[colour_ind],
                lw=0.5,
                edgecolor="k",
                label="Structures on hull" if ind == n_hulls - 2 else None,
                zorder=1e4,
            )
        if plot_points:
            ax.scatter(
                [doc["concentration"][0] for doc in hull.cursor],
                [doc[data_key][formation_energy_key][ind] for doc in hull.cursor],
                label="Structures above hull" if ind == n_hulls - 2 else None,
                alpha=1,
                marker="o",
                edgecolor="w",
                c=colours[colour_ind],
                lw=0.5,
                zorder=1e-3,
            )

    ax.set_ylim(1.1 * min_ef)

    if show_cbar:
        import matplotlib.colors

        mappable = plt.cm.ScalarMappable(
            cmap=matplotlib.colors.LinearSegmentedColormap.from_list("cut", colours),
            norm=plt.Normalize(vmin=0, vmax=np.max(hull.temperatures)),
        )
        mappable._A = hull.temperatures
        cbar = plt.colorbar(mappable, alpha=alpha)
        cbar.ax.tick_params(length=0)
        cbar.ax.yaxis.set_ticks_position("right")
        cbar.ax.set_frame_on(False)
        cbar.outline.set_visible(False)
        cbar.set_label("Temperature (K)")

    if hull.savefig or any(kwargs.get(ext) for ext in SAVE_EXTS):
        if plot_fname is not None:
            fname = plot_fname
        else:
            fname = "".join(hull.species) + data_key + "_hull"
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                plt.savefig(
                    "{}.{}".format(fname, ext), bbox_inches="tight", transparent=True
                )
                print("Wrote {}.{}".format(fname, ext))

    ax.legend()

    return ax


@plotting_function
def plot_ternary_hull(
    hull,
    axis=None,
    show=True,
    plot_points=True,
    hull_cutoff=None,
    fig_height=None,
    fig_width=None,
    label_cutoff=None,
    label_corners=True,
    expecting_cbar=True,
    labels=None,
    plot_fname=None,
    colour_points_by="hull_distance",
    hull_dist_unit="meV",
    efmap=None,
    sampmap=None,
    concmap=None,
    capmap=None,
    pathways=False,
    scale_override=None,
    label_spacing=0.05,
    label_offset=0.8,
    **kwargs,
):
    """Plot calculated ternary hull as a 2D projection.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        axis (matplotlib.axes.Axes): matplotlib axis object on which to plot.
        show (bool): whether or not to show plot in X window.
        plot_points (bool): whether or not to plot each structure as a point.
        colour_points_by (str): either `hull_distance` or `concentration`.
        label_cutoff (float/:obj:`tuple` of :obj:`float`): draw labels less than or
            between these distances form the hull, also read from hull.args.
        expecting_cbar (bool): whether or not to space out the plot to preserve
            aspect ratio if a colourbar is present.
        labels (bool): whether or not to label on-hull structures
        label_corners (bool): whether or not to put axis labels on corners or edges.
        label_offset (float): Location of the top left rcorner of the labels.
        label_spacing (float): Spacing between label list.
        hull_dist_unit (str): either "eV" or "meV",
        png/pdf/svg (bool): whether or not to write the plot to a file.
        plot_fname (str): filename to write plot to.
        efmap (bool/str): plot heatmap of formation energy,
        sampmap (bool/str): plot heatmap showing sampling density,
        concmap (bool/str): colour background by concentration and VESTA colours
        capmap (bool/str): plot heatmap showing gravimetric capacity.
        pathways (bool): plot the pathway from the starting electrode to active ion.

    Returns:
        matplotlib.axes.Axes: matplotlib axis with plot.

    """
    import ternary
    import matplotlib.pyplot as plt
    import matplotlib.colors
    from matador.utils.chem_utils import get_generic_grav_capacity

    _colour_points_by_values = ("hull_distance", "concentration")
    if colour_points_by not in _colour_points_by_values:
        raise RuntimeError(
            f"colour_points_by` must be in {_colour_points_by_values}, not {colour_points_by}"
        )

    if efmap is None:
        efmap = hull.args.get("efmap")
    if sampmap is None:
        sampmap = hull.args.get("sampmap")
    if capmap is None:
        capmap = hull.args.get("capmap")
    if pathways is None:
        pathways = hull.args.get("pathways")

    if labels is None:
        labels = hull.args.get("labels")
    if label_cutoff is None:
        label_cutoff = hull.args.get("label_cutoff")
        if label_cutoff is None:
            label_cutoff = 0
    else:
        labels = True

    _hull_cutoff = 0.1
    if hull_cutoff is not None:
        _hull_cutoff = hull_cutoff
    elif hull.hull_cutoff is not None:
        _hull_cutoff = hull.hull_cutoff
    hull_cutoff = _hull_cutoff

    print("Plotting ternary hull...")
    if capmap or efmap or concmap:
        scale = 100
    elif sampmap:
        scale = 25
    else:
        scale = 1

    if scale_override:
        scale = scale_override

    if axis is not None:
        fig, ax = ternary.figure(scale=scale, ax=axis)
    else:
        fig, ax = ternary.figure(scale=scale)

    # maintain aspect ratio of triangle
    if fig_height is None:
        _user_height = plt.rcParams.get("figure.figsize", (8, 6))[0]
    else:
        _user_height = fig_height
    if fig_width is None:
        if capmap or efmap or sampmap:
            _user_width = 8 / 5 * _user_height
        elif not expecting_cbar:
            _user_width = _user_height
        else:
            _user_width = 6.67 / 5 * _user_height
    else:
        _user_width = fig_width

    fig.set_size_inches(_user_width, _user_height)

    ax.boundary(linewidth=2.0, zorder=99)
    ax.get_axes().axis("off")
    ax.clear_matplotlib_ticks()

    chempot_labels = [
        get_formula_from_stoich(
            get_stoich_from_formula(species, sort=False), sort=False, tex=True
        )
        for species in hull.species
    ]

    ax.gridlines(color="black", multiple=scale * 0.1, linewidth=0.5)
    ticks = [float(val) for val in np.linspace(0, 1, 6)]

    element_colours = [get_element_colours()[s] for s in hull.species]
    if label_corners:
        # remove 0 and 1 ticks when labelling corners
        ticks = ticks[1:-1]
        ax.left_corner_label(
            chempot_labels[2], fontsize="x-large", color=element_colours[2]
        )
        ax.right_corner_label(
            chempot_labels[0], fontsize="x-large", color=element_colours[0]
        )
        ax.top_corner_label(
            chempot_labels[1], fontsize="x-large", color=element_colours[1], offset=0.16
        )
    else:
        ax.left_axis_label(chempot_labels[2], fontsize="large", offset=0.12)
        ax.right_axis_label(chempot_labels[1], fontsize="large", offset=0.12)
        ax.bottom_axis_label(chempot_labels[0], fontsize="large", offset=0.08)
        ax.set_title(
            "-".join(["{}".format(label) for label in chempot_labels]),
            fontsize="large",
            y=1.02,
        )

    ax.ticks(
        axis="lbr",
        linewidth=1,
        offset=0.025,
        fontsize="medium",
        locations=(scale * np.asarray(ticks)).tolist(),
        ticks=ticks,
        tick_formats="%.1f",
    )

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
            if hull_dist[ind] <= hull_cutoff:
                filtered_concs.append(conc)
                filtered_hull_dists.append(hull_dist[ind])

    unfiltered_concs = None
    if sampmap:
        unfiltered_concs = np.asarray(concs)

    concs = np.asarray(filtered_concs)
    hull_dist = np.asarray(filtered_hull_dists)

    min_cut = 0.0
    max_cut = 0.1

    if hull_dist_unit.lower() == "mev":
        hull_dist *= 1000
        min_cut *= 1000
        max_cut *= 1000

    if not hasattr(hull, "colours"):
        hull.colours = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    hull.default_cmap_list = get_linear_cmap(hull.colours[1:4], list_only=True)
    hull.default_cmap = get_linear_cmap(hull.colours[1:4], list_only=False)
    n_colours = len(hull.default_cmap_list)
    colours_hull = hull.default_cmap_list
    cmap = hull.default_cmap

    for plane in hull.convex_hull.planes:
        plane.append(plane[0])
        plane = np.asarray(plane)
        ax.plot(
            scale * plane,
            c=kwargs.get("tielinecolor", hull.colours[0]),
            lw=1.5,
            alpha=1,
            zorder=98,
        )

    if pathways:
        for phase in stable:
            if phase[0] == 0 and phase[1] != 0 and phase[2] != 0:
                ax.plot(
                    [scale * phase, [scale, 0, 0]], c="r", alpha=0.2, lw=6, zorder=99
                )

    # add points
    if plot_points:
        colours_list = []
        colour_metric = hull_dist
        conc_colours = []
        for i, _ in enumerate(colour_metric):
            if colour_points_by == "concentration":
                conc_colours.append(
                    colour_from_ternary_concentration(concs[i], hull.species)
                )
            if hull_dist[i] >= max_cut:
                colours_list.append(n_colours - 1)
            elif hull_dist[i] <= min_cut:
                colours_list.append(0)
            else:
                colours_list.append(int((n_colours - 1) * (hull_dist[i] / max_cut)))
        colours_list = np.asarray(colours_list)
        if colour_points_by == "concentration":
            colours = [
                colour_from_ternary_concentration(stable[i], hull.species)
                for i in range(len(stable))
            ]
        else:
            colours = hull.colours[1]

        # Plot stable points
        ax.scatter(
            scale * stable,
            marker="o",
            color=colours,
            edgecolors=kwargs.get("markeredgecolor", "black"),
            zorder=9999999,
            s=150,
            lw=1.5,
        )

        # Plot unstable points
        for i, _ in enumerate(concs):
            if colour_points_by == "concentration":
                c = conc_colours[i]
            else:
                c = colours_hull[colours_list[i]]

            ax.scatter(
                scale * concs[i].reshape(1, 3),
                color=c,
                marker="o",
                zorder=10000 - colours_list[i],
                s=70 * (1 - float(colours_list[i]) / n_colours) + 15,
                lw=1,
                edgecolors=kwargs.get("markeredgecolor", "black"),
            )

        # Test points for colourbar normalization
        if colour_points_by == "hull_distance":
            ax.scatter(
                scale * concs,
                colormap=cmap,
                colorbar=(colour_points_by == "hull_distance"),
                cbarlabel="Distance from hull ({}eV/atom)".format(
                    "m" if hull_dist_unit.lower() == "mev" else ""
                ),
                c=colour_metric,
                vmax=max_cut,
                vmin=min_cut,
                zorder=1000,
                s=40,
                alpha=0,
            )

        # Set up legend in case of colour by concentration
        elif colour_points_by == "concentration":
            labelled_hull_dists = [0, 25, 100]
            for dist in labelled_hull_dists:
                if dist != 0:
                    size = (
                        70 * (1 - int((n_colours - 1) * (dist / max_cut)) / n_colours)
                        + 15
                    )
                    lw = 1
                else:
                    size = 150
                    lw = 1.5
                ax.scatter(
                    [[scale, 0, 0]],
                    color="lightgrey",
                    marker="o",
                    zorder=0,
                    s=size,
                    lw=lw,
                    label=f"{dist:.0f} meV/atom",
                    edgecolors="grey",
                )

            ax.legend(title="Hull distance", fontsize="medium")

    # Add colourmaps for backgrounds
    if capmap:
        capacities = dict()
        from ternary.helpers import simplex_iterator

        for (i, j, k) in simplex_iterator(scale):
            capacities[(i, j, k)] = get_generic_grav_capacity(
                [float(i) / scale, float(j) / scale, float(scale - i - j) / scale],
                hull.species,
            )

        if not isinstance(capmap, str):
            capmap = "Pastel2"

        cmap_full = plt.cm.get_cmap(capmap)
        if cmap_full.N < 255:
            cmap_full = matplotlib.colors.LinearSegmentedColormap.from_list(
                capmap, cmap_full.colors
            )

        ax.heatmap(
            capacities,
            style="hexagonal",
            cbarlabel="Gravimetric capacity (mAh/g)",
            vmin=0,
            vmax=3000,
            cmap=cmap_full,
        )
    elif efmap:
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

        if not isinstance(efmap, str):
            efmap = "BuPu_r"

        ax.heatmap(
            energies,
            style="hexagonal",
            cbarlabel="Formation energy (eV/atom)",
            vmax=0,
            cmap=efmap,
        )
    elif sampmap:
        if not isinstance(sampmap, str):
            sampmap = "plasma"

        sampling = dict()
        from ternary.helpers import simplex_iterator

        eps = 1.0 / float(scale)
        for (i, j, k) in simplex_iterator(scale):
            sampling[(i, j, k)] = np.size(
                np.where(
                    (unfiltered_concs[:, 0] <= float(i) / scale + eps)
                    * (unfiltered_concs[:, 0] >= float(i) / scale - eps)
                    * (unfiltered_concs[:, 1] <= float(j) / scale + eps)
                    * (unfiltered_concs[:, 1] >= float(j) / scale - eps)
                    * (unfiltered_concs[:, 2] <= float(k) / scale + eps)
                    * (unfiltered_concs[:, 2] >= float(k) / scale - eps)
                )
            )
        ax.heatmap(
            sampling,
            style="hexagonal",
            cbarlabel="Number of structures",
            cmap=sampmap,
        )
    elif concmap:
        concs = dict()
        from ternary.helpers import simplex_iterator

        for (i, j, k) in simplex_iterator(scale):
            colour = colour_from_ternary_concentration(
                [float(i) / scale, float(j) / scale, float(k) / scale],
                hull.species,
                alpha=0.5,
            )
            concs[(i, j, k)] = colour
        ax.heatmap(concs, style="hexagonal", use_rgba=True, colorbar=False)

    # add labels
    if labels:
        label_cursor = _get_hull_labels(hull, label_cutoff=label_cutoff, num_species=3)
        label_coords = [
            [0.0, label_offset - i * label_spacing] for i in range(len(label_cursor))
        ]

        for ind, doc in enumerate(label_cursor):
            conc = np.asarray(doc["concentration"] + [1 - sum(doc["concentration"])])
            label = get_formula_from_stoich(
                doc["stoichiometry"],
                sort=False,
                tex=True,
                latex_sub_style=r"\mathregular",
                elements=hull.species,
            )
            if label_cutoff > 0 and doc["hull_distance"] > EPS:
                label += f' ({1000 * doc["hull_distance"]:.1f} meV/atom)'
            arrowprops = dict(
                arrowstyle="-|>", lw=2, alpha=0, zorder=1, shrinkA=2, shrinkB=4
            )
            if colour_points_by == "concentration":
                text_colour = colour_from_ternary_concentration(conc, hull.species)
            else:
                text_colour = "k"

            ax.annotate(
                label,
                scale * conc,
                textcoords="data",
                xytext=[scale * val for val in label_coords[ind]],
                ha="left",
                va="bottom",
                color=text_colour,
                arrowprops=arrowprops,
            )

    plt.tight_layout(w_pad=0.2)
    # important for retaining labels if exporting to PDF
    # see https://github.com/marcharper/python-ternary/issues/36
    ax.ax.set_aspect("equal", adjustable="box")
    ax._redraw_labels()  # noqa

    if hull.savefig or any(kwargs.get(ext) for ext in SAVE_EXTS):
        fname = plot_fname or "".join(hull.species) + "_hull"
        for ext in SAVE_EXTS:
            if hull.args.get(ext) or kwargs.get(ext):
                plt.savefig(
                    "{}.{}".format(fname, ext), bbox_inches="tight", transparent=True
                )
                print("Wrote {}.{}".format(fname, ext))
    elif show:
        print("Showing plot...")
        plt.show()

    return ax


def _scatter_plot_by_source(
    hull,
    ax,
    scale,
    kwargs,
    sources=None,
    source_labels=None,
    plot_hull_points=True,
    legend_kwargs=None,
):
    """Add scatter points to the hull depending on the guessed
    provenance of a structure.

    """
    from matador.utils.cursor_utils import get_guess_doc_provenance

    if sources is None:
        sources = [
            "AIRSS",
            "GA",
            "OQMD",
            "SWAPS",
            "ICSD",
            "DOI",
            "SM",
            "MP",
            "PF",
            "Other",
        ]

    if source_labels is None:
        source_labels = sources
    else:
        assert len(source_labels) == len(sources)

    if "Other" not in sources:
        sources.append("Other")
        source_labels.append("Other")

    # hack: double length of hull colours
    hull.colours.extend(hull.colours)

    colour_choices = {
        source: hull.colours[ind + 1] for ind, source in enumerate(sources)
    }
    points_by_source = {source: defaultdict(list) for source in sources}
    hull_points_by_source = {source: defaultdict(list) for source in sources}
    sources_present = set()
    for doc in hull.cursor:
        source = get_guess_doc_provenance(doc["source"])
        if source not in sources:
            # use grey for undesired sources
            source = "Other"
        if doc["hull_distance"] <= 0 + 2e-3 and not plot_hull_points:
            hull_points_by_source[source]["concs"].append(doc["concentration"])
            hull_points_by_source[source]["energies"].append(
                doc["formation_{}".format(hull.energy_key)]
            )
        else:
            points_by_source[source]["concs"].append(doc["concentration"])
            points_by_source[source]["energies"].append(
                doc["formation_{}".format(hull.energy_key)]
            )

        sources_present.add(source)

    alpha = kwargs.get("alpha")
    if alpha is None:
        alpha = 0.2

    legend_sources = {}

    for source in sources:
        if "concs" not in points_by_source[source]:
            continue

        concs = points_by_source[source]["concs"]
        energies = points_by_source[source]["energies"]
        ax.scatter(
            concs,
            energies,
            c=colour_choices[source],
            alpha=alpha,
            s=scale * 20,
            lw=0,
            zorder=100,
            rasterized=True,
        )

        legend_sources[source] = colour_choices[source]

    hull_point_options = dict(edgecolor="k", alpha=1, s=scale * 40, lw=1.5, zorder=1e5)

    if not plot_hull_points:
        for source in sources:
            if "concs" not in hull_points_by_source[source]:
                continue

            concs = hull_points_by_source[source]["concs"]
            energies = hull_points_by_source[source]["energies"]
            ax.scatter(
                concs, energies, facecolor=colour_choices[source], **hull_point_options
            )

            legend_sources[source] = colour_choices[source]

    for ind, source in enumerate(sources):
        if source in legend_sources:
            ax.scatter(
                1e10,
                1e10,
                facecolor=legend_sources[source],
                label=source_labels[ind],
                **hull_point_options,
            )

    if legend_kwargs is not None:
        legend = ax.legend(**legend_kwargs)
    else:
        legend = ax.legend(ncol=2)
    legend.set_zorder(1e20)

    return ax
