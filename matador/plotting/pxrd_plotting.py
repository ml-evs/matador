# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements plotting routines specifically
for the PXRD objects defined in the
matador.fingerprints.pxrd module.

"""


from matador.fingerprints.pxrd import PXRD
from matador.plotting.plotting import plotting_function
from matador.utils.cell_utils import get_space_group_label_latex
from matador.crystal import Crystal


__all__ = ['plot_pxrd']


@plotting_function
def plot_pxrd(pxrds, two_theta_range=(10, 70), labels=None, figsize=None, text_offset=0.1, **kwargs):
    """ Plot PXRD or PXRDs.

    Parameters:
        pxrds (list or matador.fingerprints.pxrd.PXRD): the PXRD
            or list of PXRDs to plot.

    Keyword arguments:
        labels (list of str): list of labels to plot alongside pattern.
        two_theta_range (tuple): plotting limits for 2theta
        figsize (tuple): specify a figure size, the default
            scales with the number of PXRDs to be plotted.

    """

    import matplotlib.pyplot as plt

    if not isinstance(pxrds, list):
        pxrds = [pxrds]
    if labels is not None and not isinstance(labels, list):
        labels = [labels]

    if figsize is None:
        height = len(pxrds) * max(0.5, 5/len(pxrds))
        figsize = (8, height)

    fig, ax = plt.subplots(figsize=figsize)

    for ind, pxrd in enumerate(pxrds):
        if isinstance(pxrd, Crystal):
            pxrd = pxrd.pxrd
        elif isinstance(pxrd, dict) and 'pxrd' in pxrd:
            pxrd = pxrd['pxrd']

        if labels:
            label = labels[ind]
        else:
            label = get_space_group_label_latex(pxrd.spg) + '-' + pxrd.formula

        ax.plot(pxrd.two_thetas, 0.9*pxrd.spectrum + ind)

        ax.text(0.95, ind+text_offset, label,
                transform=ax.get_yaxis_transform(),
                horizontalalignment='right')

    ax.set_yticks([])
    ax.set_ylim(0 - len(pxrds)*0.01, len(pxrds))
    ax.set_xlim(*two_theta_range)
    ax.set_ylabel('Relative intensity')
    ax.set_xlabel('$2\\theta$ (degrees)')

    if any([kwargs.get('pdf'), kwargs.get('svg'), kwargs.get('png')]):
        bbox_extra_artists = None
        filename = 'pxrd_plot'
        if kwargs.get('pdf'):
            plt.savefig('{}.pdf'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
        if kwargs.get('svg'):
            plt.savefig('{}.svg'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
        if kwargs.get('png'):
            plt.savefig('{}.png'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
