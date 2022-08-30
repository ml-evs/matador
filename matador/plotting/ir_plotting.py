""" This submodule implements functions useful for plotting
the results of infrared and Raman spectroscopy calculations.

"""

import numpy as np
from matador.scrapers import phonon2dict
from matador.plotting.plotting import plotting_function
from matador.utils.chem_utils import INVERSE_CM_TO_EV


@plotting_function
def plot_ir_spectrum(seed, bin_width=1.0, ax=None, show=True, **kwargs):
    """This function plots the IR/Raman spectrum found in the given .phonon file.

    Parameters:
        seed (str): the filename to scrape and plot.

    Keyword arguments:
        bin_width (float): the bin width for the IR plot.
        ax (matplotlib.axes.Axes): an existing axis on which to plot.
        show (bool): whether or not to display the plot in an X window

    Raises:
        RuntimeError: if unable to scrape IR data.

    Returns:
        matplotlib.axes.Axes: matplotlib axis with plotted data.

    """

    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax_ir = plt.subplots(figsize=kwargs.get("figsize"))
    else:
        ax_ir = ax

    ir_data, s = phonon2dict(seed)
    if not s or "infrared_intensity" not in ir_data:
        raise RuntimeError("Error scraping file: no IR intensities. {}".format(ir_data))

    wavenumbers = ir_data["eigenvalues_q"] / INVERSE_CM_TO_EV

    plotting_raman = "raman_intensity" in ir_data

    max_wavenumber = np.max(wavenumbers)
    bins = np.linspace(
        0, 1.2 * max_wavenumber, num=int((1 / bin_width) * 1.2 * max_wavenumber)
    )

    ir_spectrum, bin_edges = np.histogram(
        wavenumbers, bins=bins, weights=ir_data["infrared_intensity"], density=True
    )
    # normalize so the highest peak is at 1
    ir_spectrum /= np.max(ir_spectrum)
    if plotting_raman:
        raman_spectrum, _ = np.histogram(
            wavenumbers, bins=bins, weights=ir_data["raman_intensity"], density=True
        )
        raman_spectrum /= np.max(raman_spectrum)

    # get bin centres
    bins = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    ax_ir.plot(bins, ir_spectrum, color="#EE3425")
    ax_ir.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax_ir.set_ylabel("Relative IR intensity", color="#EE3425")
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    ir_max = np.max(ir_spectrum)
    ax_ir.set_ylim(2.1 * ir_max, -ir_max * 0.05)
    ax_ir.set_yticks(np.linspace(0, 1, 6))

    if plotting_raman:
        ax_raman = ax_ir.twinx()
        raman_max = np.max(raman_spectrum)
        ax_raman.set_yticks(np.linspace(0, 1, 6))
        ax_raman.set_ylim(-raman_max * 0.05, 2.1 * raman_max)
        ax_raman.plot(bins, raman_spectrum, color="#236DE8")
        ax_raman.set_ylabel("Relative Raman activity", color="#236DE8")

    plt.title(seed)

    if any([kwargs.get("pdf"), kwargs.get("svg"), kwargs.get("png")]):
        filename = seed.split("/")[-1].replace(".phonon", "") + "_ir"
        if kwargs.get("pdf"):
            plt.savefig(
                "{}.pdf".format(filename), bbox_inches="tight", transparent=True
            )
        if kwargs.get("svg"):
            plt.savefig(
                "{}.svg".format(filename), bbox_inches="tight", transparent=True
            )
        if kwargs.get("png"):
            plt.savefig(
                "{}.png".format(filename), bbox_inches="tight", transparent=True
            )

    elif show:
        plt.show()

    return ax
