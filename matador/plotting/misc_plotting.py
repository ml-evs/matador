from collections import defaultdict

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

from matador.plotting.plotting import get_linear_cmap, plotting_function


@plotting_function
def plot_relative_diffs(
    structure_map,
    structure_comparator,
    categories,
    field,
    ax=None,
    x_field="concentration",
    category_labels=None,
    field_label=None,
    x_field_label=None,
    colour_by_index=1,
    cbar_label=None,
    diff_type="rel",
    invert=False,
    colour_by=None,
    figsize=None,
    colourbar=True,
):
    """Plot relative and absolute differences between fields for a set of structures
    computed at different levels of theory/accuracy.

    """
    if ax is None:
        _, axes = plt.subplots(
            1, len(categories) - 1, figsize=figsize, squeeze=True, sharey=True
        )
    else:
        axes = [ax]
    if category_labels is None:
        category_labels = categories
    if field_label is None:
        field_label = f"change in {field.replace('_', ' ')}"
        if diff_type == "rel":
            field_label = "% " + field_label
    if x_field_label is None:
        x_field_label = x_field.replace("_", " ")

    if len(categories) > 2:
        raise RuntimeError("Cannot plot more than two categories simultaneously.")

    try:
        len(axes)
    except Exception:
        axes = [axes]

    for ax in axes:
        ax.set_ylabel(field_label)
        ax.set_xlabel(x_field_label)
        ax.axhline(0, ls="--", c="grey", alpha=0.5)
    x = []
    for s in structure_comparator:
        x.append(structure_map[s][categories[0]][x_field])
    if colour_by is not None:
        colour_by_x = []
        for s in structure_comparator:
            colour_by_x.append(
                structure_comparator[s]
                .get(categories[colour_by_index], {})
                .get(colour_by)
            )
    else:
        colour_by_x = None
    ys = defaultdict(list)
    if diff_type in ("abs", "rel"):
        for s in structure_comparator:
            for label in categories[1:]:
                ys[label].append(
                    structure_comparator[s].get(label, {}).get(f"{diff_type}_{field}")
                )
    else:
        for s in structure_comparator:
            for label in categories[1:]:
                ys[label].append(structure_comparator[s].get(label, {}).get(field))

    for ind, (label, y) in enumerate(ys.items()):

        if diff_type == "rel":
            y = [100 * yv if yv is not None else None for yv in y]
        if invert:
            y = [-1 * yv if yv is not None else None for yv in y]

        colour_by_x = np.array(colour_by_x, dtype=np.float64)
        x = np.array(x)[np.where(~np.isnan(colour_by_x))]
        y = np.array(y)[np.where(~np.isnan(colour_by_x))]
        colour_by_x = colour_by_x[np.where(~np.isnan(colour_by_x))]

        order = np.argsort(colour_by_x)[::-1]
        y = y[order]
        x = x[order]
        colour_by_x = colour_by_x[order]

        norm = None
        ticks = None
        default_cmap = "viridis"
        extend = None
        if colour_by == "hull_distance":
            norm = matplotlib.colors.LogNorm(0.01, 1, clip=True)
            ticks = [0.0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64]
            colours = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])
            default_cmap = get_linear_cmap(colours[1:4], list_only=False)
            extend = "both"
            hull_points = np.where(colour_by_x <= 0.0)

            if x_field == "concentration":
                x = np.array(x).flatten()
                y = np.array(y).flatten()

                axes[ind].plot(
                    x[hull_points][np.argsort(x[hull_points])],
                    y[hull_points][np.argsort(x[hull_points])],
                    c="k",
                    zorder=1e10,
                    lw=2,
                    markeredgecolor="k",
                    markerfacecolor=colours[1],
                    marker="o",
                    markersize=np.sqrt(40),
                    markeredgewidth=1.5,
                    alpha=1,
                )

            else:
                axes[ind].scatter(
                    np.array(x)[hull_points],
                    np.array(y)[hull_points],
                    c=colours[1],
                    zorder=1e10,
                    lw=1.5,
                    edgecolor="k",
                    norm=norm,
                    alpha=1,
                )

        scatter = axes[ind].scatter(
            x,
            y,
            cmap=default_cmap,
            c=colour_by_x,
            zorder=-1,
            alpha=0.8,
        )
        axes[ind].set_title(f"{category_labels[0]} vs {category_labels[ind+1]}")

    if colour_by_x is not None and colourbar:
        cbar = plt.colorbar(scatter, aspect=30, pad=0.02, extend=extend, ticks=ticks)
        if ticks:
            cbar.ax.set_yticklabels(ticks)
            cbar.ax.tick_params(which="minor", length=0)
        if cbar_label is None and colour_by is not None:
            cbar_label = colour_by.replace("_", " ")
        cbar.set_label(cbar_label)

    if len(axes) == 1:
        axes = axes[0]

    return axes
