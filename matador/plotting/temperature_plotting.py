# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot temperature-dependent
quantities for particular structures, e.g. free energies and heat
capacities.

"""

from matador.plotting.plotting import plotting_function

__all__ = ["plot_free_energy"]


@plotting_function
def plot_free_energy(model, temperatures=None, ax=None, line_kwargs=None, **kwargs):
    """Plot G(T) on the array of given temperatures. Default T is [0, 800].

    Keyword arguments:
        temperatures (list/np.ndarray): list or array of temperatures to plot.
            If the array/list has length 2, use these as the start and endpoints
            with 21 plotting points.
        ax (matplotlib.pyplot.Axis): axis object to plot onto.

    """

    import numpy as np
    import matplotlib.pyplot as plt

    legend = isinstance(line_kwargs, dict) and "label" in line_kwargs and ax is not None

    if line_kwargs is None:
        line_kwargs = {}

    num_points = 21
    start = 0
    end = 800
    if temperatures is not None and len(temperatures) == 2:
        start = temperatures[0]
        end = temperatures[1]

    if temperatures is None or len(temperatures) == 2:
        _temperatures = np.linspace(start, end, num_points)
    elif temperatures is not None:
        _temperatures = temperatures

    if ax is None:
        fig, ax = plt.subplots()

    t, e = model.vibrational_free_energy(temperatures=_temperatures)
    ax.plot(t, e, **line_kwargs)
    ax.set_ylabel("Vibrational free energy (eV/atom)")
    ax.set_xlabel("Temperature (K)")
    if legend:
        ax.legend()

    return ax
