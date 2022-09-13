# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot battery-relevant curves,
such as voltages and volume expansions.

"""

from typing import List, Optional, Dict, Union

import numpy as np
from matador.utils.chem_utils import (
    get_formula_from_stoich,
    get_subscripted_formula_tex,
)
from matador.plotting.plotting import plotting_function, SAVE_EXTS
from matador.battery import VoltageProfile

EPS = 1e-12


__all__ = ["plot_voltage_curve", "plot_volume_curve"]


@plotting_function
def plot_voltage_curve(
    profiles: Union[List[VoltageProfile], VoltageProfile],
    ax=None,
    show: bool = False,
    labels: bool = False,
    savefig: Optional[str] = None,
    curve_labels: Optional[Union[str, List[str]]] = None,
    line_kwargs: Optional[Union[Dict, List[Dict]]] = None,
    shade_average_voltage: bool = False,
    colour_by_stoichiometry: bool = False,
    plot_average_voltage: bool = False,
    plot_max_capacity: bool = False,
    legend: bool = True,
    label_suffix: str = None,
    expt: Optional[str] = None,
    expt_label: Optional[str] = None,
):
    """Plot voltage curve calculated for phase diagram.

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
        shade_average_voltage (bool): whether to shade the region spanned from average
            voltage to max capacity.
        plot_average_voltage (bool): whether to draw a guideline at the average voltage.
        plot_max_capacity (bool): whether to draw a guideline at the maximum capacity.
        legend (bool): whether to draw the legend.
        label_suffix (str): A suffix to add to every label on the plot (useful when combining axes)
        colour_by_stoichiometry (bool): whether to use mixed VESTA colours for each curve
            based on the starting formula of the electrode.
        expt (str): string to a filename of a csv Q, V to add to the plot.
        expt_label (str): label for any experimental profile passed to the plot.

    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    if not isinstance(profiles, list):
        profiles = [profiles]
    if curve_labels is not None and not isinstance(curve_labels, list):
        curve_labels = [curve_labels]
    if line_kwargs is not None and not isinstance(line_kwargs, list):
        line_kwargs = len(profiles) * [line_kwargs]

    if curve_labels is not None and len(curve_labels) != len(profiles):
        raise RuntimeError(
            "Wrong number of labels passed for number of profiles: {} vs {}".format(
                len(curve_labels), len(profiles)
            )
        )

    if line_kwargs is not None and len(line_kwargs) != len(profiles):
        raise RuntimeError(
            "Wrong number of line kwargs passed for number of profiles: {} vs {}".format(
                len(line_kwargs), len(profiles)
            )
        )

    dft_label = None
    if expt is not None:
        expt_data = np.loadtxt(expt, delimiter=",")
        if expt_label:
            ax.plot(
                expt_data[:, 0], expt_data[:, 1], c="k", lw=2, ls="-", label=expt_label
            )
        else:
            ax.plot(
                expt_data[:, 0],
                expt_data[:, 1],
                c="k",
                lw=2,
                ls="-",
                label="Experiment",
            )

        if len(profiles) == 1:
            dft_label = "DFT (this work)"

    line_colours = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])

    for ind, profile in enumerate(profiles):
        if dft_label is None and curve_labels is None:
            stoich_label = get_formula_from_stoich(
                profile.starting_stoichiometry, tex=True
            )
        else:
            stoich_label = None

        _label = stoich_label if dft_label is None else dft_label
        if curve_labels is not None and len(curve_labels) > ind:
            _label = curve_labels[ind]

        if label_suffix:
            _label += f" {label_suffix}"

        _line_kwargs = {"c": line_colours[(ind % len(line_colours)) + 2]}
        if line_kwargs is not None:
            _line_kwargs.update(line_kwargs[ind])

        if colour_by_stoichiometry:
            from matador.utils.viz_utils import formula_to_colour

            colour = formula_to_colour(
                get_formula_from_stoich(profile.starting_stoichiometry)
            )
            _line_kwargs["c"] = colour
        else:
            colour = _line_kwargs["c"]

        _add_voltage_curve(
            profile.capacities, profile.voltages, ax, label=_label, **_line_kwargs
        )

        if shade_average_voltage:
            ax.fill_between(
                np.linspace(
                    0,
                    np.nanmax(profile.capacities),
                    2,
                ),
                np.linspace(0, 0, 2),
                [profile.average_voltage, profile.average_voltage],
                lw=1,
                edgecolor=colour,
                facecolor=colour + "20",
                zorder=-ind + 2,
                ls="--",
            )

        if plot_average_voltage:
            vline_kwargs = _line_kwargs
            vline_kwargs["ls"] = "--"
            vline_kwargs["alpha"] = 0.8
            vline_kwargs["zorder"] = 0
            ax.axhline(profile.average_voltage, **vline_kwargs)

        if plot_max_capacity:
            vline_kwargs = _line_kwargs
            vline_kwargs["ls"] = "-."
            vline_kwargs["alpha"] = 0.8
            vline_kwargs["zorder"] = 0
            ax.axvline(np.nanmax(profile.capacities), **vline_kwargs)

    if labels:
        if len(profiles) > 1:
            print("Only labelling first voltage profile.")
        for ind, reaction in enumerate(profiles[0].reactions[:-1]):
            _labels = []
            for phase in reaction:
                if phase[0] is None or phase[0] == 1.0:
                    _label = ""
                else:
                    _label = "{:.1f} ".format(phase[0])
                _label += "{}".format(get_subscripted_formula_tex(phase[1]))
                _labels.append(_label)
            _label = "+".join(_labels)
            _position = (
                profiles[0].capacities[ind + 1],
                profiles[0].voltages[ind + 1] + max(profiles[0].voltages) * 0.01,
            )
            ax.annotate(
                _label, xy=_position, textcoords="data", ha="center", zorder=9999
            )

    if expt or len(profiles) > 1:

        if plot_average_voltage:
            ax.axvline(-1, c="k", ls="--", alpha=0.8, label="Average voltages")
        if plot_max_capacity:
            ax.axvline(
                -1, c="k", ls="-.", alpha=0.8, label="Maximum gravimetric capacity"
            )

        if legend:
            ax.legend()

    ax.set_ylabel("Voltage (V) vs {ion}$^+/${ion}".format(ion=profile.active_ion))
    ax.set_xlabel("Gravimetric cap. (mAh/g)")

    _, end = ax.get_ylim()
    from matplotlib.ticker import MultipleLocator

    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.set_ylim(0, 1.1 * end)
    _, end = ax.get_xlim()
    ax.set_xlim(0, 1.1 * end)
    ax.grid(False)
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)

    if savefig:
        plt.savefig(savefig)
        print("Wrote {}".format(savefig))

    elif show:
        plt.show()

    return ax


def _add_voltage_curve(capacities, voltages, ax, label=None, **kwargs):
    """Add the voltage curves stored under hull['voltage_data'] to the plot.

    Parameters:
        capacities (list): list or numpy array of capacities.
        voltages (list): list or numpy array of voltages.
        ax (matplotlib.axes.Axes): an existing axis object on which to plot.

    Keyword arguments:
        **kwargs (dict): to pass to matplotlib, using abbreviated names (e.g. 'c' not 'color')
        label (str): if present, add this label to the first segment of the voltage curve
            so it can be added to the legend.
        alpha (float): transparency of line

    """
    line_kwargs = {
        "lw": 2,
        "alpha": 1,
    }
    line_kwargs.update(kwargs)

    for i in range(1, len(voltages) - 1):
        if i == 1 and label is not None:
            ax.plot(
                [capacities[i - 1], capacities[i]],
                [voltages[i], voltages[i]],
                label=label,
                **line_kwargs,
            )
        else:
            ax.plot(
                [capacities[i - 1], capacities[i]],
                [voltages[i], voltages[i]],
                **line_kwargs,
            )
        if i != len(voltages) - 2:
            ax.plot(
                [capacities[i], capacities[i]],
                [voltages[i], voltages[i + 1]],
                **line_kwargs,
            )


@plotting_function
def plot_volume_curve(
    hull,
    ax=None,
    show=True,
    legend=False,
    as_percentages=False,
    label=None,
    label_suffix=None,
    exclude_elemental=False,
    line_kwargs: Optional[Union[Dict, List[Dict]]] = None,
    colour_by_stoichiometry: bool = False,
    end_marker_only: bool = True,
    **kwargs,
):
    """Plot volume curve calculated for phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        ax (matplotlib.axes.Axes): an existing axis on which to plot.
        show (bool): whether or not to display plot in X-window.
        legend (bool): whether to add the legend.
        label (str/list): Either a list of labels or a single label for the volume curves.
        label_suffix (str): A suffix to add to every label on the plot (useful when combining axes)
        as_percentages (bool): whether to show the expansion as a percentage.
        exclude_elemental (bool): whether to include any elemental electrodes
            when considering ternary phase diagrams.
        line_kwargs (list or dict): parameters to pass to the curve plotter,
            if a list then the line kwargs will be passed to each line individually.
        colour_by_stoichiometry (bool): whether to use mixed VESTA colours for each curve
            based on the starting formula of the electrode.
        end_marker_only (bool): Only place a marker on the worst volume expansion for each curve.

    """
    import matplotlib.pyplot as plt

    if ax is None:
        if hull.savefig or any([kwargs.get(ext) for ext in SAVE_EXTS]):
            fig = plt.figure()
        else:
            fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        ax = ax

    if label and not isinstance(label, list):
        label = [label]

    num_curves = len(hull.volume_data["Q"])

    if line_kwargs is not None and not isinstance(line_kwargs, list):
        line_kwargs = num_curves * [line_kwargs]

    if as_percentages:
        volume_key = "volume_expansion_percentage"
    else:
        volume_key = "volume_ratio_with_bulk"

    line_colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for j in range(num_curves):
        if exclude_elemental and len(hull.volume_data["endstoichs"][j]) == 1:
            continue

        _line_kwargs = {"c": line_colours[(j % len(line_colours)) + 2]}
        if line_kwargs is not None:
            _line_kwargs.update(line_kwargs[j])

        if colour_by_stoichiometry:
            from matador.utils.viz_utils import formula_to_colour

            colour = formula_to_colour(
                get_formula_from_stoich(hull.volume_data["endstoichs"][j])
            )
            _line_kwargs["c"] = colour
        else:
            colour = _line_kwargs["c"]

        lw = _line_kwargs.pop("lw", 2)
        marker = _line_kwargs.pop("marker", "o")
        markeredgewidth = _line_kwargs.pop("markeredgewidth", 1.5)
        markeredgecolor = _line_kwargs.pop("markeredgecolor", "k")

        stable_hull_dist = hull.volume_data["hull_distances"][j]
        if len(stable_hull_dist) != len(hull.volume_data["Q"][j]):
            raise RuntimeError("This plot does not support --hull_cutoff.")

        if end_marker_only:
            ax.plot(
                [0, np.max(hull.volume_data["Q"][j])],
                [1, np.max(hull.volume_data[volume_key][j])],
                marker=marker,
                markeredgewidth=markeredgewidth,
                markeredgecolor=markeredgecolor,
                zorder=1000,
                lw=0,
                **_line_kwargs,
            )

        else:
            ax.plot(
                [
                    q
                    for ind, q in enumerate(hull.volume_data["Q"][j])
                    if stable_hull_dist[ind] <= EPS
                ],
                [
                    v
                    for ind, v in enumerate(hull.volume_data[volume_key][j])
                    if stable_hull_dist[ind] <= EPS
                ],
                marker=marker,
                markeredgewidth=markeredgewidth,
                markeredgecolor=markeredgecolor,
                zorder=1000,
                lw=0,
                **_line_kwargs,
            )

        if label:
            _label = label[j]
        else:
            _label = get_formula_from_stoich(
                hull.volume_data["endstoichs"][j], tex=True
            )
        if label_suffix:
            _label += f" {label_suffix}"

        ax.plot(
            [
                q
                for ind, q in enumerate(hull.volume_data["Q"][j])
                if stable_hull_dist[ind] <= EPS
            ],
            [
                v
                for ind, v in enumerate(hull.volume_data[volume_key][j])
                if stable_hull_dist[ind] <= EPS
            ],
            label=_label,
            lw=lw,
            **_line_kwargs,
        )

    ax.set_xlabel("Gravimetric capacity (mAh/g)")
    if as_percentages:
        ax.set_ylabel("Volume expansion (%)")
    else:
        if len(hull.volume_data["endstoichs"]) == 1:
            ax.set_ylabel(
                f"Volume ratio with {get_formula_from_stoich(hull.volume_data['endstoichs'][0], tex=True)} starting electrode"
            )
        else:
            ax.set_ylabel("Volume ratio with starting electrode")

    if legend or len(hull.volume_data["Q"]) > 1:
        ax.legend()

    return ax
