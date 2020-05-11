# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements plotting routines specifically
for the PDF and PDFOverlap objects defined in the
matador.fingerprints.pdf module.

"""

import numpy as np

from matador.fingerprints.pdf import PDF
from matador.crystal import Crystal
from matador.plotting.plotting import plotting_function
from matador.utils.cell_utils import get_space_group_label_latex
from matador.utils.chem_utils import get_formula_from_stoich

__all__ = ['plot_pdf', 'plot_projected_pdf', 'plot_diff_overlap', 'plot_projected_diff_overlap']


@plotting_function
def plot_pdf(pdfs,
             labels=None, r_min=None, r_max=None,
             offset=1.2, text_offset=(0.0, 0.0),
             legend=False, annotate=True, figsize=None,
             filename=None,
             **kwargs):
    """ Plot PDFs.

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

    Returns:
        matplotlib.pyplot.Axes: axis object which can be modified further.

    """

    import matplotlib.pyplot as plt

    if not isinstance(pdfs, list):
        pdfs = [pdfs]
    if labels is not None and not isinstance(labels, list):
        labels = [labels]

    if figsize is None:
        _user_default_figsize = plt.rcParams.get('figure.figsize', (8, 6))
        height = len(pdfs) * max(0.5, _user_default_figsize[1] / 1.5 / len(pdfs))
        figsize = (_user_default_figsize[0], height)

    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(111)

    if labels is not None and len(labels) != len(pdfs):
        raise RuntimeError("Wrong number of labels {} for PDFs.".format(labels))

    if isinstance(pdfs[0], Crystal):
        gr_max = max(np.max(pdf.pdf.gr) for pdf in pdfs)
        _r_max = min(np.max(pdf.pdf.r_space) for pdf in pdfs)
    elif isinstance(pdfs[0], dict):
        gr_max = max(np.max(pdf['pdf'].gr) for pdf in pdfs)
        _r_max = min(np.max(pdf['pdf'].r_space) for pdf in pdfs)
    else:
        gr_max = max(np.max(pdf) for pdf in pdfs)
        _r_max = min(np.max(pdf.r_space) for pdf in pdfs)
    abs_offset = offset * gr_max

    if r_max is None:
        r_max = _r_max
    if r_min is None:
        r_min = 0.0

    ax1.set_ylabel('Pair distribution function, $g(r)$')
    ax1.get_yaxis().set_ticks([])
    ax1.set_xlim(r_min, r_max+0.5)

    for ind, pdf in enumerate(pdfs):

        if isinstance(pdf, Crystal):
            pdf = pdf.pdf
        elif isinstance(pdf, dict) and 'pdf' in pdf:
            pdf = pdf['pdf']

        if labels:
            label = labels[ind]
        else:
            label = get_space_group_label_latex(pdf.spg) + '-' + get_formula_from_stoich(pdf.stoichiometry, tex=True)

        ax1.plot(pdf.r_space, pdf.gr + abs_offset * ind, label=label)
        if text_offset is not None:
            text_x = text_offset[0]
        if text_offset is not None:
            text_y = abs_offset*ind + text_offset[1]*gr_max
        if label is not None and annotate:
            ax1.text(text_x, text_y, label)

    ax1.set_ylim(-gr_max * 0.2, offset * gr_max * len(pdfs))

    ax1.set_xlabel('$r$ ($\\AA$)')

    if legend:
        legend = ax1.legend()

    if any([kwargs.get('pdf'), kwargs.get('svg'), kwargs.get('png')]):
        bbox_extra_artists = None
        if filename is None:
            filename = '-'.join([get_formula_from_stoich(pdf.stoichiometry) for pdf in pdfs]) + '_pdf'

        if kwargs.get('pdf'):
            plt.savefig('{}.pdf'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
        if kwargs.get('svg'):
            plt.savefig('{}.svg'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
        if kwargs.get('png'):
            plt.savefig('{}.png'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)

    return ax1


@plotting_function
def plot_projected_pdf(pdf, keys=None, other_pdfs=None, vlines=None):
    """ Plot projected PDFs.

    Parameters:
        pdf (matador.fingerprints.pdf.PDF): the main PDF to plot.

    Keyword arguments:
        keys (list): plot only a subset of projections, e.g. [('K', )].
        other_pdfs (list of PDF): other PDFs to plot.
        vlines (list of float): plot vertical lines at these points.

    """
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    if keys is None:
        keys = [key for key in pdf.elem_gr]
    for key in keys:
        if key not in pdf.elem_gr:
            key = (key[1], key[0])
        ax1.plot(pdf.r_space, pdf.elem_gr[key], label='-'.join(key) + ' {}'.format(pdf.label))
    if other_pdfs is not None:
        if isinstance(other_pdfs, PDF):
            other_pdfs = [other_pdfs]
        for _pdf in other_pdfs:
            if isinstance(_pdf, PDF):
                for key in keys:
                    ax1.plot(_pdf.r_space, _pdf.elem_gr[key], ls='--',
                             label='-'.join(key) + ' {}'.format(_pdf.label))
            elif isinstance(pdf, tuple):
                ax1.plot(_pdf[0], _pdf[1], alpha=1, ls='--')
            else:
                raise RuntimeError

    if vlines is not None:
        for line in vlines:
            ax1.axvline(line, ls='--', alpha=0.8, c='grey')
    ax1.legend(loc=1)
    ax1.set_ylabel('$g(r)$')
    ax1.set_xlabel('$r$ (Angstrom)')


@plotting_function
def plot_diff_overlap(pdf_overlap):
    """ Simple plot for comparing two PDFs.

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

    ax2.set_xlabel('$r$ (\\AA)')
    ax1.set_ylabel('$g(r)$')
    ax2.set_ylabel('$g_a(r) - g_b(r)$')
    ax2.axhline(0, ls='--', c='k', lw=0.5)
    ax1.set_xlim(0, np.max(pdf_overlap.fine_space))

    ax1.plot(pdf_overlap.fine_space, pdf_overlap.fine_gr_a, label=pdf_overlap.pdf_a.label)
    ax1.plot(pdf_overlap.fine_space, pdf_overlap.fine_gr_b, label=pdf_overlap.pdf_b.label)

    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.set_ylim(-0.5 * ax1.get_ylim()[1], 0.5 * ax1.get_ylim()[1])

    ax1.legend(loc=0)
    ax2.plot(pdf_overlap.fine_space, pdf_overlap.overlap_fn, ls='-')
    ax2.set_ylim(ax1.get_ylim()[1], ax1.get_ylim()[1])


@plotting_function
def plot_projected_diff_overlap(pdf_overlap):
    """ Simple plot for comparing two PDFs.

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
    ax2.set_xlabel('$r$ (\\AA)')
    ax1.set_ylabel('$g(r)$')
    ax2.set_ylabel('$g_a(r) - g_b(r)$')
    ax2.axhline(0, ls='--', c='k', lw=0.5)
    ax1.set_xlim(0, np.max(pdf_overlap.fine_space))
    for _, key in enumerate(pdf_overlap.fine_elem_gr_a):
        ax1.plot(pdf_overlap.fine_space, pdf_overlap.fine_elem_gr_a[key],
                 label='-'.join(key) + ' {}'.format(pdf_overlap.pdf_a.label))
        ax1.plot(pdf_overlap.fine_space, pdf_overlap.fine_elem_gr_b[key],
                 label='-'.join(key) + ' {}'.format(pdf_overlap.pdf_b.label),
                 ls='--')
        ax2.plot(pdf_overlap.fine_space, pdf_overlap.fine_elem_gr_a[key] - pdf_overlap.fine_elem_gr_b[key],
                 label='-'.join(key) + ' diff')
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.set_ylim(ax1.get_ylim()[1], ax1.get_ylim()[1])
    ax1.legend(loc=0)
    ax2.legend(loc=2)
