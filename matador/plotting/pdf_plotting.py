# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements plotting routines specifically
for the PDF and PDFOverlap objects defined in the
matador.fingerprints.pdf module.

"""

import numpy as np

from matador.crystal import Crystal
from matador.plotting.plotting import plotting_function
from matador.utils.cell_utils import get_space_group_label_latex
from matador.utils.chem_utils import get_formula_from_stoich

__all__ = [
    "plot_pdf",
    "plot_projected_pdf",
    "plot_diff_overlap",
    "plot_projected_diff_overlap",
]


@plotting_function
def plot_pdf(
    pdfs,
    labels=None,
    r_min=None,
    r_max=None,
    offset=1.2,
    text_offset=(0.0, 0.0),
    legend=False,
    annotate=True,
    figsize=None,
    ax=None,
    projected=False,
    colour_labels=True,
    **kwargs,
):
    """Plot PDFs.

    Parameters:
        pdfs (list of matador.fingerprints.pdf.PDF or matador.crystal.Crystal or dict):
            the PDFs to plot, as a list of PDF or Crystal objects, or a matador document.

    Keyword arguments:
        labels (list of str): labels to add to the PDF plot.
        offset (float): amount by which to separate the PDFs in the plot. A value of 1
            will separate by the maximum intensity across the PDFs. Default is 1.5.
        text_offset (tuple of float): two float values to move annotations around relative
            to the base of their corresponding PDF, in units of (Angstrom, max_gr).
        r_max (float): the radius to plot out to. Default is the minmax(radius across
            all PDFs).
        annotate (bool): whether or not to apply the PDF labels as an annotation.
        legend (bool): whether or not to apply the PDF labels as a legend.
        figsize (tuple of float): matplotlib figure size. Default scales with number of PDFs.
        ax (matplotlib.Axis): optional axis object to plot on.
        projected (list(str)): if provided or True, will plot the PDFs projected onto the given keys.
        colour_labels (bool): whether to colour the labels based on the PDF's colour.

    Returns:
        matplotlib.pyplot.Axes: axis object which can be modified further.

    """

    import matplotlib.pyplot as plt
    from matador.utils.viz_utils import get_element_colours

    if not isinstance(pdfs, list):
        pdfs = [pdfs]
    if labels is not None and not isinstance(labels, list):
        labels = [labels]

    if figsize is None and ax is None:
        _user_default_figsize = plt.rcParams.get("figure.figsize", (8, 6))
        height = len(pdfs) * max(0.5, _user_default_figsize[1] / 1.5 / len(pdfs))
        figsize = (_user_default_figsize[0], height)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if labels is not None and len(labels) != len(pdfs):
        raise RuntimeError("Wrong number of labels {} for PDFs.".format(labels))

    if isinstance(pdfs[0], Crystal):
        gr_max = max(np.max(pdf.pdf.gr) for pdf in pdfs)
        _r_max = min(np.max(pdf.pdf.r_space) for pdf in pdfs)
    elif isinstance(pdfs[0], dict):
        gr_max = max(np.max(pdf["pdf"].gr) for pdf in pdfs)
        _r_max = min(np.max(pdf["pdf"].r_space) for pdf in pdfs)
    else:
        gr_max = max(np.max(pdf.gr) for pdf in pdfs)
        _r_max = min(np.max(pdf.r_space) for pdf in pdfs)
    abs_offset = offset * gr_max

    if r_max is None:
        r_max = _r_max
    if r_min is None:
        r_min = 0.0

    ax.set_ylabel("Pair distribution function, $g(r)$")
    ax.get_yaxis().set_ticks([])
    ax.set_xlim(r_min, r_max + 0.5)

    projected_keys = set()
    labelled_keys = set()

    for ind, pdf in enumerate(pdfs):

        if isinstance(pdf, Crystal):
            pdf = pdf.pdf
        elif isinstance(pdf, dict) and "pdf" in pdf:
            pdf = pdf["pdf"]

        if labels:
            label = labels[ind]
        else:
            label = f"{get_formula_from_stoich(pdf.stoichiometry, tex=True)}-{get_space_group_label_latex(pdf.spg)}"

        if projected:
            if isinstance(projected, bool):
                keys = [key for key in pdf.elem_gr]
                [projected_keys.add(key) for key in keys]
            else:
                for key in projected:
                    if key in pdf.elem_gr:
                        projected_keys.add(key)
                    elif ([key[1], key[0]]) in pdf.elem_gr:
                        projected_keys.add((key[1], key[0]))

            for key in projected_keys:
                if key in pdf.elem_gr:
                    _label = None
                    if key not in labelled_keys:
                        if len(key) == 2:
                            _label = f"{key[0]}-{key[1]}"
                        else:
                            _label = f"{key[0]}-{key[0]}"
                        labelled_keys.add(key)
                    if len(key) == 2:
                        color = (
                            np.array(get_element_colours()[key[0]])
                            + np.array(get_element_colours()[key[1]])
                        ) / 2
                    else:
                        color = get_element_colours()[key[0]]
                    ax.plot(
                        pdf.r_space,
                        pdf.elem_gr[key] + abs_offset * ind,
                        label=_label,
                        c=color,
                    )

        else:
            color = next(ax._get_lines.prop_cycler)["color"]
            ax.plot(pdf.r_space, pdf.gr + abs_offset * ind, label=label, color=color)
        if text_offset is not None:
            text_x = text_offset[0] + r_min
        if text_offset is not None:
            text_y = abs_offset * ind + text_offset[1] * gr_max
        if label is not None and annotate:
            text_color = None
            if colour_labels:
                text_color = color
            ax.text(text_x, text_y, label, color=text_color)

    ax.set_ylim(-gr_max * 0.2, offset * gr_max * len(pdfs))

    ax.set_xlabel("$r$ ($\\AA$)")

    if legend or projected:
        legend = ax.legend(c="upper right")

    return ax


@plotting_function
def plot_projected_pdf(*args, **kwargs):
    """DEPRECATED"""
    return plot_pdf(*args, **kwargs)


@plotting_function
def plot_diff_overlap(pdf_overlap):
    """Simple plot for comparing two PDFs.

    Parameters:
        pdf_overlap (matador.fingerprints.pdf.PDFOverlap): the
        overlap object to plot.

    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import numpy as np

    plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    gs.update(hspace=0)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)

    ax2.set_xlabel("$r$ (\\AA)")
    ax1.set_ylabel("$g(r)$")
    ax2.set_ylabel("$g_a(r) - g_b(r)$")
    ax2.axhline(0, ls="--", c="k", lw=0.5)
    ax1.set_xlim(0, np.max(pdf_overlap.fine_space))

    ax1.plot(
        pdf_overlap.fine_space, pdf_overlap.fine_gr_a, label=pdf_overlap.pdf_a.label
    )
    ax1.plot(
        pdf_overlap.fine_space, pdf_overlap.fine_gr_b, label=pdf_overlap.pdf_b.label
    )

    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.set_ylim(-0.5 * ax1.get_ylim()[1], 0.5 * ax1.get_ylim()[1])

    ax1.legend(loc=0)
    ax2.plot(pdf_overlap.fine_space, pdf_overlap.overlap_fn, ls="-")
    ax2.set_ylim(ax1.get_ylim()[1], ax1.get_ylim()[1])


@plotting_function
def plot_projected_diff_overlap(pdf_overlap):
    """Simple plot for comparing two PDFs.

    Parameters:
        pdf_overlap (matador.fingerprints.pdf.PDFOverlap): the
            overlap object to plot.

    """
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.gridspec as gridspec

    plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    gs.update(hspace=0)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.set_xlabel("$r$ (\\AA)")
    ax1.set_ylabel("$g(r)$")
    ax2.set_ylabel("$g_a(r) - g_b(r)$")
    ax2.axhline(0, ls="--", c="k", lw=0.5)
    ax1.set_xlim(0, np.max(pdf_overlap.fine_space))
    for _, key in enumerate(pdf_overlap.fine_elem_gr_a):
        ax1.plot(
            pdf_overlap.fine_space,
            pdf_overlap.fine_elem_gr_a[key],
            label="-".join(key) + " {}".format(pdf_overlap.pdf_a.label),
        )
        ax1.plot(
            pdf_overlap.fine_space,
            pdf_overlap.fine_elem_gr_b[key],
            label="-".join(key) + " {}".format(pdf_overlap.pdf_b.label),
            ls="--",
        )
        ax2.plot(
            pdf_overlap.fine_space,
            pdf_overlap.fine_elem_gr_a[key] - pdf_overlap.fine_elem_gr_b[key],
            label="-".join(key) + " diff",
        )
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.set_ylim(ax1.get_ylim()[1], ax1.get_ylim()[1])
    ax1.legend(loc=0)
    ax2.legend(loc=2)
