# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements plotting routines specifically
for the PXRD objects defined in the
matador.fingerprints.pxrd module.

"""


from matador.plotting.plotting import plotting_function
from matador.utils.cell_utils import get_space_group_label_latex
from matador.crystal import Crystal


__all__ = ["plot_pxrd"]


@plotting_function
def plot_pxrd(
    pxrds,
    two_theta_range=None,
    text_loc="right",
    rug=False,
    rug_height=0.05,
    rug_offset=0.04,
    offset=None,
    ax=None,
    labels=None,
    figsize=None,
    text_offset=0.1,
    colour_labels=True,
    filename=None,
    **kwargs,
):
    """Plot PXRD or PXRDs.

    Parameters:
        pxrds (list or matador.fingerprints.pxrd.PXRD): the PXRD
            or list of PXRDs to plot.

    Keyword arguments:
        two_theta_range (tuple): plotting limits for 2theta
        rug (bool): whether to provide a rug plot of all peaks.
        rug_height (float): size of rug ticks.
        rug_offset (float): offset of rug ticks.
        offset (float): extra space added between patterns (as fraction
            of max intensity). Default 0.1.
        labels (list of str): list of labels to plot alongside pattern.
        figsize (tuple): specify a figure size, the default
            scales with the number of PXRDs to be plotted.
        ax (matplotlib.Axis): optional axis object to plot on.
        text_offset (float): amount by which to offset the labels.
        text_loc (str): location of text, either 'left' or 'right'.
        filename (str): optional filename for saving.


    """

    import matplotlib.pyplot as plt

    if not isinstance(pxrds, list):
        pxrds = [pxrds]
    if labels is not None and not isinstance(labels, list):
        labels = [labels]

    if figsize is None and ax is None:
        _user_default_figsize = plt.rcParams.get("figure.figsize", (8, 6))
        height = len(pxrds) * max(0.5, _user_default_figsize[1] / 1.5 / len(pxrds))
        figsize = (_user_default_figsize[0], height)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if offset is None and rug:
        offset = 0.2
    elif offset is None:
        offset = 0.1

    colour_cycle = ax._get_lines.prop_cycler

    for ind, pxrd in enumerate(pxrds):
        if isinstance(pxrd, Crystal):
            pxrd = pxrd.pxrd
        elif isinstance(pxrd, dict) and "pxrd" in pxrd:
            pxrd = pxrd["pxrd"]

        c = next(colour_cycle).get("color")

        if labels:
            label = labels[ind]
        else:
            label = f"{pxrd.formula}-{get_space_group_label_latex(pxrd.spg)}"

        ax.plot(pxrd.two_thetas, (1 - offset) * pxrd.pattern + ind, c=c)

        if colour_labels:
            text_colour = c
        else:
            text_colour = None

        if text_loc == "right":
            ax.text(
                0.95,
                ind + text_offset,
                label,
                transform=ax.get_yaxis_transform(),
                horizontalalignment="right",
                color=text_colour,
            )

        if text_loc == "left":
            ax.text(
                0.05,
                ind + text_offset,
                label,
                transform=ax.get_yaxis_transform(),
                horizontalalignment="left",
                color=text_colour,
            )

        if rug:
            import numpy as np

            peaks = np.unique(pxrd.peak_positions)
            for peak in peaks:
                ax.plot(
                    [peak, peak],
                    [ind - rug_height - rug_offset, ind - rug_offset],
                    c=c,
                    alpha=0.5,
                )

    ax.set_ylim(-0.2, len(pxrds) + 0.1)
    if two_theta_range is not None:
        ax.set_xlim(*two_theta_range)
    ax.set_ylabel("Relative intensity")
    ax.set_xlabel("$2\\theta$ (degrees)")

    return ax
