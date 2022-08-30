# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot chemical shifts
and chemical shieldings from NMR calculations.

"""

from typing import List, Optional, Dict, Union, Tuple

import numpy as np
from matador.crystal import Crystal
from matador.plotting.plotting import plotting_function
from matador.fingerprints import Fingerprint

__all__ = ["plot_magres"]


@plotting_function
def plot_magres(
    magres: Union[List[Crystal], Crystal],
    species: str,
    magres_key: str = "chemical_shielding_iso",
    xlabel: str = None,
    broadening_width: float = 1,
    text_offset: float = 0.2,
    ax=None,
    figsize: Tuple[float, float] = None,
    show: bool = False,
    savefig: Optional[str] = None,
    signal_labels: Optional[Union[str, List[str]]] = None,
    signal_limits: Tuple[float, float] = None,
    line_kwargs: Optional[Union[Dict, List[Dict]]] = None,
    invert_xaxis: bool = True,
    species_label: Optional[str] = None,
    label_fontsize: int = None,
):
    """Plot magnetic resonance data for a set of crystals.

    Parameters:
        magres (Union[Crystal, List[Crystal]]): list of :class:`Crystal` containing
            magres data.
        species (str | Tuple[str, int]): the species to plot the shifts of, either as
            a simple element symbol, or an element and isotope tuple used for labelling,
            e.g., "P" or ("P", 31).

    Keyword arguments:
        ax (matplotlib.axes.Axes): an existing axis on which to plot.
        magres_key (str): the data key for which the magres site data is stored under.
        show (bool): whether to show plot in an X window.
        figsize (Tuple[float]): overrides the default size for the matplotlib figure.
        broadening_width (float): the Lorentzian width to apply to the shifts.
        xlabel (str): a custom label for the x-axis.
        savefig (str): filename to use to save the plot.
        signal_labels (list): optional list of labels for the curves in
            the magres list.
        signal_limits (Tuple[float]): values at which to clip the magres signals. Defaults
            to the maximum and minimum shifts across all passed structures.
        line_kwargs (list or dict): parameters to pass to the curve plotter,
            if a list then the line kwargs will be passed to each line individually.
        label_fontsize (int): font size override for the labels.

    """
    import matplotlib.pyplot as plt

    if not isinstance(magres, list):
        magres = [magres]
    if signal_labels is not None and not isinstance(signal_labels, list):
        signal_labels = [signal_labels]
    if line_kwargs is not None and not isinstance(line_kwargs, list):
        line_kwargs = len(magres) * [line_kwargs]

    if figsize is None:
        _user_default_figsize = plt.rcParams.get("figure.figsize", (8, 6))
        height = len(magres) * max(0.5, _user_default_figsize[1] / 1.5 / len(magres))
        figsize = (_user_default_figsize[0], height)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

    if species is None:
        raise RuntimeError("You must provide the species of interest plotting.")
    if isinstance(species, tuple):
        species_label = f"$^{{ {species[1]} }}${species[0]}"
        species = species[0]
    else:
        species_label = species

    if signal_labels is not None and len(signal_labels) != len(magres):
        raise RuntimeError(
            "Wrong number of labels passed for number of magres: {} vs {}".format(
                len(signal_labels), len(magres)
            )
        )

    _magres = []

    if signal_limits is not None:
        min_shielding, max_shielding = sorted(signal_limits)
        if signal_limits[0] > signal_limits[1]:
            ax.invert_xaxis()
    else:
        min_shielding, max_shielding = (1e20, -1e20)

    for ind, doc in enumerate(magres):

        if isinstance(doc, dict):
            _doc = Crystal(doc)
        else:
            _doc = doc

        _magres.append(_doc)

        relevant_sites = [atom for atom in _doc if atom.species == species]
        if relevant_sites:
            shielding = [atom[magres_key] for atom in relevant_sites]
            if signal_limits is None:
                min_shielding = min(np.min(shielding), min_shielding)
                max_shielding = max(np.max(shielding), max_shielding)

    if min_shielding > 1e19 and max_shielding < -1e19:
        raise RuntimeError(
            f"No sites of {species} found in any of the passed crystals."
        )

    _buffer = 0.2 * np.abs(min_shielding - max_shielding)
    s_space = np.linspace(min_shielding - _buffer, max_shielding + _buffer, num=1000)

    _padded_colours = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    _padded_colours = (1 + (len(magres) // len(_padded_colours))) * _padded_colours

    if line_kwargs is not None and len(line_kwargs) != len(magres):
        raise RuntimeError(
            "Wrong number of line kwargs passed for number of magres: {} vs {}".format(
                len(line_kwargs), len(magres)
            )
        )

    for ind, doc in enumerate(_magres):
        if signal_labels is None:
            stoich_label = doc.formula_tex
        else:
            stoich_label = None

        _label = stoich_label
        if signal_labels is not None and len(signal_labels) > ind:
            _label = signal_labels[ind]

        _line_kwargs = {"c": _padded_colours[ind]}
        if line_kwargs is not None:
            _line_kwargs.update(line_kwargs[ind])

        relevant_sites = [site for site in doc if site.species == species]
        if not relevant_sites:
            print(
                f"No sites of {species} found in {doc.root_source}, signal will be empty."
            )
            signal = np.zeros_like(s_space)

        else:
            shifts = [site[magres_key] for site in relevant_sites]

            hist, bins = np.histogram(shifts, bins=s_space)

            if broadening_width > 0:
                signal = Fingerprint._broadening_unrolled(
                    hist, s_space, broadening_width, broadening_type="lorentzian"
                )
            else:
                signal = np.array(hist, dtype=np.float64)
                bin_centres = s_space[:-1] + (s_space[1] - s_space[0]) / 2
                s_space = bin_centres

            if np.max(signal) > 1e-10:
                signal /= np.max(signal)
            else:
                signal *= 0

        ax.plot(s_space, signal + (ind * 1.1), **_line_kwargs)

        if _label is not None:
            ax.text(
                0.95,
                (ind * 1.1) + text_offset,
                _label,
                transform=ax.get_yaxis_transform(),
                horizontalalignment="right",
                fontsize=label_fontsize,
            )

    if xlabel is None:
        unit = set(doc.get("magres_units", {}).get("ms", "ppm") for doc in magres)
        if len(unit) > 1:
            raise RuntimeError(
                f"Multiple incompatible units found for chemical shift: {unit}"
            )
        unit = list(unit)[0]
        if magres_key == "chemical_shielding_iso":
            xlabel = f"{species_label}: Isotropic chemical shielding $\\sigma_\\mathrm{{iso}}$ ({unit})"
        elif magres_key == "chemical_shift_iso":
            xlabel = f"{species_label}: Isotropic chemical shift $\\delta_\\mathrm{{iso}}$ ({unit})"
        elif magres_key == "chemical_shift_aniso":
            xlabel = f"{species_label}: Anisotropic chemical shift $\\Delta $ ({unit})"
        elif magres_key == "chemical_shift_asymmetry":
            xlabel = f"{species_label}: Chemial shift asymmetry, $\\eta$"

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Intensity (arb. units)")

    if len(magres) > 1:
        ax.set_yticks([])
    else:
        ax.set_yticks(np.linspace(0, 1, 5, endpoint=True))

    ax.set_ylim(-0.1, 1.1 * len(magres))
    if invert_xaxis:
        ax.invert_xaxis()

    if savefig:
        plt.savefig(savefig)
        print("Wrote {}".format(savefig))

    elif show:
        plt.show()

    return ax
