# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements plotting routines specifically
for the PDF and PDFOverlap objects defined in the
matador.fingerprints.pdf module.

"""


from matador.fingerprints.pdf import PDF
from matador.plotting.plotting import plotting_function


@plotting_function
def plot_pdf(pdf, other_pdfs=None):
    """ Plot PDFs.

    Parameters:
        pdf (matador.fingerprints.pdf.PDF): the main PDF to plot.

    Keyword arguments:
        other_pdfs (list of PDF): other PDFs to add to the plot.

    """
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot(111)
    ax1.plot(pdf.r_space, pdf.gr, label=pdf.label)
    ax1.set_ylabel('Pair distribution function, $g(r)$')
    ax1.set_xlim(0, pdf.rmax)
    if other_pdfs is not None:
        if not isinstance(other_pdfs, list):
            other_pdfs = [other_pdfs]
        for _pdf in other_pdfs:
            # check if we have a Crystal object without importing
            try:
                _pdf = _pdf.pdf
            except AttributeError:
                pass
            # check if we have a normal matador doc
            if isinstance(_pdf, dict) and 'pdf' in _pdf:
                _pdf = _pdf['pdf']
            if isinstance(_pdf, PDF):
                ax1.plot(_pdf.r_space, _pdf.gr, label=_pdf.label, ls='--', alpha=1)
            elif isinstance(_pdf, tuple):
                ax1.plot(_pdf[0], _pdf[1], alpha=1, ls='--')
            else:
                raise RuntimeError('Wrong PDF format specified, please either pass a PDF object or (r, g(r)) tuple.')
    ax1.set_xlabel('$r$ (Angstrom)')
    ax1.legend()


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
    fig = plt.figure(figsize=(10, 6))
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

    plt.figure(figsize=(8, 6))
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

    plt.figure(figsize=(8, 6))
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
