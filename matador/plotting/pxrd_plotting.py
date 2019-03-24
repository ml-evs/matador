# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements plotting routines specifically
for the PXRD objects defined in the
matador.similarity.pxrd module.

"""


from matador.similarity.pxrd import PXRD
from matador.plotting.plotting import plotting_function


@plotting_function
def plot_pxrd(pxrds):
    """ Plot PXRD or PXRDs.

    Parameters:
        pxrds (list or matador.similarity.pxrd.PXRD): the PXRD
            or list of PXRDs to plot.

    """
    if isinstance(pxrds, PXRD):
        pxrds = [pxrds]
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, 6))
    for ind, pxrd in enumerate(pxrds):
        ax = fig.add_subplot(111)
        ax.plot(pxrd.two_thetas, pxrd.spectrum + ind, label='{} ({})'.format(pxrd.formula, pxrd.spg))
    ax.set_ylabel('Relative intensity')
    ax.set_xlabel('$2\\theta$')
    ax.legend()
