# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot densities of states and
bandstructures for electronic and vibrational calculations.

"""


import os
import copy
import numpy as np

import matplotlib.pyplot as plt

from matador.utils.viz_utils import get_element_colours
from matador.plotting.plotting import plotting_function
from matador.scrapers import optados2dict, phonon2dict, bands2dict, phonon_dos2dict
from matador.scrapers import cell2dict, res2dict
from matador.orm.spectral import (
    ElectronicDispersion,
    ElectronicDOS,
    VibrationalDispersion,
    VibrationalDOS,
    Dispersion,
    DensityOfStates,
)
from matador.utils.chem_utils import INVERSE_CM_TO_EV

__all__ = ["plot_spectral", "dos_plot", "dispersion_plot"]

PROJECTOR_MIN = 1e-5


@plotting_function
def plot_spectral(seeds, **options):
    """This function wraps all of the spectral plotting capability of matador.
    When provided with a seed, or seeds, several files will be checked:

        - <seed>.bands: CASTEP bandstructure expected, not DOS,
        - <seed>.adaptive.dat: OptaDOS total DOS,
        - <seed>.pdos.dat: OptaDOS pDOS,
        - <seed>.pdis.dat: OptaDOS projected dispersion,

    or, if the "phonons" flag is passed, or if a <seed>.phonon file is detected,

        - <seed>.phonon: CASTEP phonon dispersion curves expected,
        - <seed>.phonon_dos: CASTEP phonon DOS.

    This function will then automatically choose what to plot, favouring a bandstructure
    with "the-most-projected" DOS it can find.

    Parameters:
        seeds (list): list of filenames of bands/phonon files

    Keyword Arguments:
        plot_bandstructure (bool): whether to plot bandstructure, if available
        plot_dos (bool): whether to plot density of states, if available
        plot_pdos (bool): whether or not to plot projected DOS, if available
        plot_pdis (bool): whether to plot projected dispersion, if available
        dos (str): separate seed name for pDOS/DOS data
        phonons (bool): whether to plot phonon or electronic data
        labels (list): list of strings for legend labels for multiple bandstructures
        gap (bool): whether to draw on the band gap
        external_efermi (float or list): replace scraped Fermi energy with this value (eV) (can be
            specified per spin channel).
        highlight_bands (list): list of integer indices, colour the bands
            with these indices in red
        band_colour (str): if passed "occ", bands will be coloured using
            cmap depending on whether they lie above or below the Fermi
            level. Otherwise, override all colour options with
            matplotlib-interpretable colour (e.g. hexcode or html colour
            name) to use for all bands (DEFAULT: 'occ').
        band_alpha (float): transparency of plotted bands.
        filename (str): filename for figure saving.
        cmap (str): matplotlib colourmap name to use for the bands
        cmap_limits (tuple): fraction of cmap to use (DEFAULT: (0.2, 0.8)).
        unstacked_pdos (bool): whether to plot projected DOS as stack or overlapping.
        spin_only (str): either 'up' or 'down' to only plot one spin channel.
        preserve_kspace_distance (bool): whether to preserve distances in reciprocal space when
            linearising the kpoint path. If False, bandstructures of different lattice parameters
            with the same Bravais lattice can be more easily compared. If True, bandstructures may
            appear rarefied or compressed in particular regions.
        pdis_interpolation_factor (float): multiple by which to interpolate pDIS bands
        pdis_point_scale (float): size of points in pDIS (DEFAULT: 50).
        colours (list of str): list of matplotlib colours to override default colour cycle for projectors and otherwise.
        spin_up_colour (str): matplotlib colour to override the spin-up DOS colour.
        spin_down_colour (str): matplotlib colour to override the spin-up DOS colour.
        projectors_to_plot (str): comma-separted list of projectors to
            plot in the PDIS or PDOS, provided as element:orbital, e.g.
            "K:s,K:p,P" will plot s and p orbitals for K, and all orbitals for P.
        band_reorder (bool): try to reorder bands based on local gradients (DEFAULT: True for phonons, otherwise False).
        title (str): optional plot title
        pdos_hide_sum (bool): whether or not to plot the total DOS on a PDOS plot; this is to hide
            regions where the PDOS is negative (leading to total DOS lower than stacked PDOS) (DEFAULT: False).

    """
    from cycler import cycler

    # set defaults and update class with desired values
    prop_defaults = {
        "plot_bandstructure": True,
        "plot_dos": True,
        "plot_pdos": True,
        "plot_pdis": True,
        "phonons": False,
        "gap": False,
        "external_efermi": None,
        "labels": None,
        "cmap": None,
        "cmap_limits": (0.2, 0.8),
        "band_colour": None,
        "spin_only": None,
        "figsize": None,
        "filename": None,
        "pdis_interpolation_factor": 2,
        "pdis_point_scale": 25,
        "projectors_to_plot": None,
        "projector_colours": None,
        "colours": None,
        "unstacked_pdos": False,
        "preserve_kspace_distance": False,
        "band_reorder": False,
        "title": None,
        "show": True,
        "verbosity": 0,
        "highlight_bands": None,
        "pdos_hide_sum": True,
        "spin_up_colour": "firebrick",
        "spin_down_colour": "dodgerblue",
    }

    for key in options:
        if options[key] is not None:
            prop_defaults[key] = options[key]
    options = prop_defaults

    options["projectors_to_plot"] = _parse_projectors_list(
        options["projectors_to_plot"]
    )
    options["colour_cycle"] = tuple(plt.rcParams["axes.prop_cycle"].by_key()["color"])

    if options["projector_colours"] is not None:
        options["colours"] = options["projector_colours"]

    if options["colours"] is not None:
        options["colour_cycle"] = options["colours"]

    if options.get("cmap") is None:
        plt.rcParams["axes.prop_cycle"] = cycler("color", options["colour_cycle"])
    else:
        print("Adjusting colour palette... to {}".format(options.get("cmap")))
        try:
            options["colours"] = plt.cm.get_cmap(options.get("cmap")).colors
            plt.rcParams["axes.prop_cycle"] = cycler("color", options["colours"])
        except AttributeError:
            options["colours"] = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        options["_mpl_cmap"] = plt.get_cmap(options.get("cmap"))

    if (
        options.get("phonons")
        and options.get("cmap") is None
        and options.get("colours") is None
    ):
        options["band_colour"] = options.get("band_colour", "grey")
        options["band_alpha"] = options.get("band_alpha", 0.8)

    if not isinstance(seeds, list):
        seeds = [seeds]

    if len(seeds) > 1:
        if options["plot_pdis"] or options["plot_dos"]:
            options["plot_pdos"] = False
            options["plot_pdis"] = False
            print("Disabling projections as mutiple seeds requested.")

    if options.get("plot_window") is not None:
        if isinstance(options.get("plot_window"), list):
            if len(options.get("plot_window")) == 1:
                options["plot_window"] = (
                    -options["plot_window"][0],
                    options["plot_window"][0],
                )
            elif len(options.get("plot_window")) != 2:
                raise RuntimeError(
                    f"`plot_window` must have length 2 or be a single number, not {options['plot_window']}"
                )
        else:
            options["plot_window"] = (-options["plot_window"], options["plot_window"])

        options["plot_window"] = sorted(options.get("plot_window"))

    else:
        options["plot_window"] = None

    if all(isinstance(seed, str) for seed in seeds):

        if options["plot_dos"]:
            # check an optados file exists
            exts = [
                "pdos.dat",
                "adaptive.dat",
                "fixed.dat",
                "linear.dat",
                "jdos.dat",
                "phonon_dos",
                "bands_dos",
            ]
            options["plot_dos"] = any(
                [
                    any([os.path.isfile("{}.{}".format(seed, ext)) for ext in exts])
                    for seed in seeds
                ]
            )

        if options["plot_pdos"]:
            exts = ["pdos.dat", "phonon_dos"]
            options["plot_pdos"] = any(
                [
                    any([os.path.isfile("{}.{}".format(seed, ext)) for ext in exts])
                    for seed in seeds
                ]
            )

    figsize = options["figsize"]

    if options["plot_bandstructure"] and not options["plot_dos"]:
        if figsize is None:
            figsize = (7, 6)
        fig, ax_dispersion = plt.subplots(figsize=figsize)
    elif options["plot_bandstructure"] and options["plot_dos"]:
        if figsize is None:
            figsize = (10, 6)
        fig, ax_grid = plt.subplots(
            1,
            3,
            figsize=figsize,
            sharey=True,
            gridspec_kw={"width_ratios": [4, 2, 1], "wspace": 0.1, "left": 0.15},
        )
        ax_dispersion = ax_grid[0]
        ax_dos = ax_grid[1]
        ax_grid[2].axis("off")
    elif not options["plot_bandstructure"] and options["plot_dos"]:
        if figsize is None:
            figsize = (9, 4)
        fig, ax_dos = plt.subplots(1, figsize=figsize)
    else:
        raise RuntimeError(
            "No plots requested, please set either plot_dos or plot_bandstructure to True!"
        )

    options["valence"] = options["colour_cycle"][0]
    options["conduction"] = options["colour_cycle"][-1]
    options["crossing"] = options["colour_cycle"][int(len(options["colour_cycle"]) / 2)]

    if len(seeds) > 1:
        options["ls"] = ["-"] * len(seeds)
        if options.get("labels") is None:
            try:
                options["labels"] = [
                    seed.split("/")[-1].split(".")[0] for seed in seeds
                ]
            except AttributeError:
                options["labels"] = [seed.root_source for seed in seeds]

        if len(options.get("labels", [])) != len(seeds):
            raise RuntimeError(
                f"Invalid number of labels provided for {len(seeds)} seeds: {options.get('labels')}. "
                "Multiple labels should be comma separated."
            )

        options["labels"] = [label.strip() for label in options["labels"]]

    options["ls"] = []
    for i in range(len(seeds)):
        if i % 3 == 0:
            options["ls"].append("-")
        elif i % 3 == 1:
            options["ls"].append("--")
        elif i % 3 == 2:
            options["ls"].append("-.")

    bbox_extra_artists = []
    if options["plot_bandstructure"]:
        ax_dispersion = dispersion_plot(
            seeds, ax_dispersion, options, bbox_extra_artists
        )

    if options["plot_dos"]:
        ax_dos = dos_plot(seeds, ax_dos, options, bbox_extra_artists)

    if options.get("title") is not None:
        fig.suptitle(options.get("title"))

    if any([options.get("pdf"), options.get("svg"), options.get("png")]):
        if not bbox_extra_artists:
            bbox_extra_artists = None
        filename = options.get("filename")
        if filename is None:
            filename = (
                seeds[0].split("/")[-1].replace(".bands", "").replace(".phonon", "")
                + "_spectral"
            )
        if options.get("pdf"):
            plt.savefig(
                "{}.pdf".format(filename),
                bbox_inches="tight",
                transparent=True,
                bbox_extra_artists=bbox_extra_artists,
            )
        if options.get("svg"):
            plt.savefig(
                "{}.svg".format(filename),
                bbox_inches="tight",
                transparent=True,
                bbox_extra_artists=bbox_extra_artists,
            )
        if options.get("png"):
            plt.savefig(
                "{}.png".format(filename),
                bbox_inches="tight",
                transparent=True,
                bbox_extra_artists=bbox_extra_artists,
            )

    else:
        plt.tight_layout()
        if options["show"]:
            print("Displaying plot...")
            plt.show()


@plotting_function
def dispersion_plot(seeds, ax_dispersion, options, bbox_extra_artists=None):
    """Plot a dispersion/bandstructure on the given axis. Will detect
    and projected dispersion data automatically.

    Parameters:
        seeds (str or list): the seedname(s) of the data to plot.
        ax_dispersion (matplotlib.Axes): the axis to plot on.
        options (dict): any plotting keywords (from e.g. dispersion script).
        bbox_extra_artists (list): a list to which to append legends.

    Returns:
        matplotlib.Axes: the axis that was plotted on.

    """
    plotted_pdis = False

    if not isinstance(seeds, list):
        seeds = [seeds]

    if bbox_extra_artists is None:
        bbox_extra_artists = []

    for seed_ind, seed in enumerate(seeds):

        if isinstance(seed, dict):
            if options.get("phonons"):
                dispersion = VibrationalDispersion(seed)
            else:
                dispersion = ElectronicDispersion(seed)

        elif isinstance(seed, Dispersion):
            dispersion = seed

        else:
            seed = seed.replace(".bands", "").replace(".phonon", "")
            if os.path.isfile("{}.phonon".format(seed)):
                dispersion, s = phonon2dict(
                    seed + ".phonon", verbosity=options.get("verbosity")
                )
                if not s:
                    raise RuntimeError(dispersion)

                dispersion = VibrationalDispersion(dispersion)

            elif os.path.isfile("{}.bands".format(seed)):
                dispersion, s = bands2dict(
                    seed + ".bands", verbosity=options.get("verbosity")
                )
                if not s:
                    raise RuntimeError(dispersion)

                if os.path.isfile("{}.pdis.dat".format(seed)) and options.get(
                    "plot_pdis"
                ):
                    pdis_data, s = optados2dict("{}.pdis.dat".format(seed))
                    if not s:
                        raise RuntimeError(pdis_data)
                else:
                    pdis_data = None

                dispersion = ElectronicDispersion(dispersion, projection_data=pdis_data)

            else:
                raise RuntimeError("{}.bands/.phonon not found.".format(seed))

        eigs = np.array(dispersion.eigs, copy=True)
        if options.get("phonons"):
            # convert from internal eV frequencies to cm^-1
            eigs /= INVERSE_CM_TO_EV

        if options.get("plot_window") is None:
            if options.get("phonons"):
                options["plot_window"] = [min(-10, np.min(eigs) - 10), np.max(eigs)]
            else:
                options["plot_window"] = [-10, 10]

        # try to match bands if requested
        if options.get("band_reorder"):
            if options.get("plot_pdis"):
                raise RuntimeError(
                    "PDIS not currently compatible with post-hoc band reordering."
                )
            print("Reordering bands based on local gradients...")
            eigs = Dispersion.get_band_reordering(eigs, dispersion.kpoint_branches)

        path = dispersion.linearise_path(
            preserve_kspace_distance=options.get("preserve_kspace_distance")
        )

        if (
            dispersion.projectors
            and len(seeds) == 1
            and options.get("plot_pdis")
            and not options.get("phonons")
        ):
            ax_dispersion = projected_bandstructure_plot(
                dispersion,
                ax_dispersion,
                path,
                bbox_extra_artists,
                eigs=eigs,
                **options,
            )
            options["band_colour"] = "grey"
            plotted_pdis = True

        # loop over branches and plot
        if not plotted_pdis:
            if options.get("external_efermi") is None:
                spin_fermi_energy = dispersion.spin_fermi_energy
            else:
                spin_fermi_energy = options.get("external_efermi")
            if len(spin_fermi_energy) == 1 and dispersion.num_spins != 1:
                spin_fermi_energy = [spin_fermi_energy] * dispersion.num_spins

            if options.get("cmap") is not None:
                cmap_limits = options.get("cmap_limits", (0.2, 0.8))
                options["_mpl_cmap"] = plt.cm.get_cmap(options.get("cmap"))(
                    np.linspace(*cmap_limits, num=dispersion.num_bands)
                )

            # loop over branches and plot
            for branch_ind, branch in enumerate(dispersion.kpoint_branches):
                for ns in range(dispersion.num_spins):
                    if ns == 1 and options.get("spin_only") == "up":
                        continue
                    elif ns == 0 and options.get("spin_only") == "down":
                        continue

                    for nb in range(dispersion.num_bands):
                        colour, alpha, label = _get_lineprops(
                            dispersion,
                            spin_fermi_energy,
                            nb,
                            ns,
                            branch,
                            branch_ind,
                            seed_ind,
                            options,
                            eigs=eigs,
                        )

                        ax_dispersion.plot(
                            path[(np.asarray(branch) - branch_ind).tolist()],
                            eigs[ns][nb][branch] - spin_fermi_energy[ns],
                            c=colour,
                            ls=options.get("ls", len(seeds) * ["-"])[seed_ind],
                            alpha=alpha,
                            label=label,
                        )

    if len(seeds) > 1 and options.get("legend", True):
        disp_legend = ax_dispersion.legend(loc="upper center")
        bbox_extra_artists.append(disp_legend)

    ax_dispersion.axhline(0, ls="--", lw=1, c="grey")
    ax_dispersion.set_ylim(options.get("plot_window"))
    if options.get("phonons"):
        ylabel = "Wavenumber (cm$^{-1}$)"
    else:
        ylabel = "Energy (eV)"
    ax_dispersion.set_ylabel(ylabel)
    ax_dispersion.set_xlim(-0.05, 1.05)
    _add_path_labels(seeds[-1], dispersion, ax_dispersion, path, 0, options)

    return ax_dispersion


@plotting_function
def dos_plot(seeds, ax_dos, options, bbox_extra_artists=None):
    """Plot a density of states on the given axis. Will detect
    pDOS and spin-dependent DOS data automatically.

    Parameters:
        seeds (list): the seednames of the data to plot.
        ax_dos (matplotlib.Axes): the axis to plot on.
        options (dict): any plotting keywords (from e.g. dispersion script).
        bbox_extra_artists (list): a list to which to append legends.

    Returns:
        matplotlib.Axes: the axis that was plotted on.

    """

    if bbox_extra_artists is None:
        bbox_extra_artists = []

    if not isinstance(seeds, list):
        seeds = [seeds]

    for seed_ind, seed in enumerate(seeds):
        if not options.get("phonons"):
            dos_data = _load_electronic_dos(seed, options)

            if options.get("plot_window") is None:
                options["plot_window"] = [-10, 10]
        else:
            dos_data = _load_phonon_dos(seed, options)
            max_density = np.max(dos_data["dos"])

        if options.get("plot_pdos") and "pdos" in dos_data:
            pdos_data = dos_data["pdos"]

        energies = np.copy(dos_data["energies"])
        # change unit of phonon energies and set plot window
        if options.get("phonons"):
            energies /= INVERSE_CM_TO_EV
            if "pdos" in dos_data:
                dos_data["pdos"]["energies"] /= INVERSE_CM_TO_EV

            if options.get("plot_window") is None:
                options["plot_window"] = [
                    np.min(energies[np.where(dos_data["dos"] > 1e-3)]) - 10,
                    np.max(energies[np.where(dos_data["dos"] > 1e-3)]),
                ]

        dos = dos_data["dos"]

        if "spin_dos" in dos_data:
            max_density = max(
                np.max(
                    np.abs(
                        dos_data["spin_dos"]["down"][
                            np.where(energies > options["plot_window"][0])
                        ]
                    )
                ),
                np.max(
                    np.abs(
                        dos_data["spin_dos"]["up"][
                            np.where(energies > options["plot_window"][0])
                        ]
                    )
                ),
            )
        else:
            max_density = np.max(
                dos_data["dos"][
                    np.where(
                        np.logical_and(
                            energies < options["plot_window"][1],
                            energies > options["plot_window"][0],
                        )
                    )
                ]
            )

        # plotting pdos depends on these other factors too
        plotting_pdos = (
            options.get("plot_pdos")
            and len(seeds) == 1
            and not (options.get("phonons") and len(dos_data.get("pdos", [])) <= 1)
        )

        if options.get("phonons"):
            ylabel = "Phonon DOS"
            xlabel = "Wavenumber (cm$^{{-1}}$)"
        else:
            if "dos_unit_label" in dos_data:
                ylabel = dos_data["dos_unit_label"].replace("A^3", "Å$^{3}$")
            else:
                if options.get("plot_bandstructure"):
                    ylabel = "DOS"
                else:
                    ylabel = "DOS (eV$^{{-1}}$Å$^{{-3}}$)"
            xlabel = "Energy (eV)"

        if options.get("plot_bandstructure"):
            ax_dos.set_xlabel(ylabel)
            ax_dos.axhline(0, c="grey", ls="--", lw=1)
            if "spin_dos" in dos_data:
                ax_dos.set_xlim(-max_density * 1.2, max_density * 1.2)
            else:
                ax_dos.set_xlim(0, max_density * 1.2)
            ax_dos.set_ylim(options.get("plot_window"))
            ax_dos.axvline(0, c="grey", lw=1)
            ax_dos.xaxis.set_ticks_position("none")

            if "spin_dos" not in dos_data:
                ax_dos.plot(
                    dos,
                    energies,
                    ls=options.get("ls", len(seeds) * ["-"])[seed_ind],
                    color="grey",
                    zorder=1e10,
                    label="Total DOS",
                )
                if not options.get("plot_pdos"):
                    ax_dos.fill_betweenx(
                        energies[np.where(energies > 0)],
                        0,
                        dos[np.where(energies > 0)],
                        alpha=options.get("fill_alpha", 0.2),
                        color=options.get("conduction"),
                    )
                    ax_dos.fill_betweenx(
                        energies[np.where(energies <= 0)],
                        0,
                        dos[np.where(energies <= 0)],
                        alpha=options.get("fill_alpha", 0.2),
                        color=options.get("valence"),
                    )
        else:
            ax_dos.set_xlabel(xlabel)
            ax_dos.set_ylabel(ylabel)
            ax_dos.axvline(0, c="grey", lw=1, ls="--")
            if "spin_dos" in dos_data:
                ax_dos.set_ylim(-max_density * 1.2, max_density * 1.2)
            else:
                ax_dos.set_ylim(0, max_density * 1.2)
            ax_dos.set_xlim(options.get("plot_window"))
            ax_dos.axhline(0, c="grey", lw=1)

            if "spin_dos" not in dos_data:

                dos_colour = options["colour_cycle"][seed_ind]
                if len(seeds) > 1:
                    label = options.get("labels")[seed_ind]
                else:
                    label = "Total DOS"

                ax_dos.plot(
                    energies,
                    dos,
                    ls=options.get("ls", len(seeds) * ["-"])[seed_ind],
                    alpha=1,
                    c=dos_colour,
                    lw=options.get("lw", 1),
                    zorder=1e10,
                    label=label,
                )

                if not plotting_pdos:
                    ax_dos.fill_between(
                        energies,
                        0,
                        dos,
                        alpha=options.get("fill_alpha", 0.2),
                        color=dos_colour,
                    )

        if "spin_dos" in dos_data and not options.get("pdos_hide_sum"):
            if options.get("plot_bandstructure"):
                if options.get("spin_only") in [None, "up"]:
                    if not plotting_pdos:
                        ax_dos.fill_betweenx(
                            energies,
                            0,
                            dos_data["spin_dos"]["up"],
                            alpha=options.get("fill_alpha", 0.2),
                            color=options["spin_up_colour"],
                        )
                    ax_dos.plot(
                        dos_data["spin_dos"]["up"],
                        energies,
                        ls=options.get("ls", len(seeds) * ["-"])[seed_ind],
                        color=options["spin_up_colour"],
                        zorder=1e10,
                        label="$\\uparrow$",
                    )
                if options.get("spin_only") in [None, "down"]:
                    if not plotting_pdos:
                        ax_dos.fill_betweenx(
                            energies,
                            0,
                            -dos_data["spin_dos"]["down"],
                            alpha=0.2,
                            color=options["spin_down_colour"],
                        )
                    ax_dos.plot(
                        -dos_data["spin_dos"]["down"],
                        energies,
                        ls=options.get("ls", len(seeds) * ["-"])[seed_ind],
                        color=options["spin_down_colour"],
                        zorder=1e10,
                        label="$\\downarrow$",
                    )
            else:
                if options.get("spin_only") in [None, "up"]:
                    ax_dos.plot(
                        energies,
                        dos_data["spin_dos"]["up"],
                        ls=options.get("ls", len(seeds) * ["-"])[seed_ind],
                        color=options["spin_up_colour"],
                        zorder=1e10,
                        label="$\\uparrow$",
                    )
                    if not plotting_pdos:
                        ax_dos.fill_between(
                            energies,
                            0,
                            dos_data["spin_dos"]["up"],
                            alpha=0.2,
                            color=options["spin_up_colour"],
                        )
                if options.get("spin_only") in [None, "down"]:
                    ax_dos.plot(
                        energies,
                        -dos_data["spin_dos"]["down"],
                        ls=options.get("ls", len(seeds) * ["-"])[seed_ind],
                        color=options["spin_down_colour"],
                        zorder=1e10,
                        label="$\\downarrow$",
                    )
                    if not plotting_pdos:
                        ax_dos.fill_between(
                            energies,
                            0,
                            -dos_data["spin_dos"]["down"],
                            alpha=0.2,
                            color=options["spin_down_colour"],
                        )

        if plotting_pdos:

            if options.get("projectors_to_plot") is not None:
                pdos = dict()
                for projector in pdos_data["pdos"]["projectors"]:
                    if projector in options.get("projectors_to_plot"):
                        pdos[projector] = pdos_data["pdos"][projector]
            else:
                pdos = pdos_data["pdos"]

            pdos_energies = pdos_data["energies"]

            stacks = dict()
            projector_labels, dos_colours = _get_projector_info(
                list(pdos.keys()),
                colours_override=options["colour_cycle"]
                if options.get("colours")
                else None,
            )
            unique_labels = set()
            for ind, projector in enumerate(pdos):

                # don't break PDOS label down by spin
                if projector_labels[ind] in unique_labels:
                    projector_labels[ind] = ""
                else:
                    unique_labels.add(projector_labels[ind])

                # split stacked pdos by spin channel
                stack_key = None
                if len(projector) > 2:
                    stack_key = projector[2]

                if stack_key not in stacks:
                    stacks[stack_key] = np.zeros_like(pdos[projector])

                stack = stacks[stack_key]
                if options.get("unstacked_pdos"):
                    stack = 0
                else:
                    stack = stacks[stack_key]

                if not options.get("unstacked_pdos"):
                    alpha = 0.8
                else:
                    alpha = 0.7

                # flip sign of down spin energies for spin polarised plot
                if "down" in projector:
                    pdos[projector] *= -1

                if not np.max(np.abs(pdos[projector])) < 1e-8:
                    if options.get("plot_bandstructure"):
                        label = None
                        if not options.get("unstacked_pdos"):
                            ax_dos.fill_betweenx(
                                pdos_energies,
                                stack,
                                stack + pdos[projector],
                                alpha=alpha,
                                lw=0,
                                label=projector_labels[ind],
                                color=dos_colours[ind],
                            )
                            projector_outline_alpha = 0
                        else:
                            label = projector_labels[ind]
                            projector_outline_alpha = 1

                        ax_dos.plot(
                            stack + pdos[projector],
                            pdos_energies,
                            alpha=projector_outline_alpha,
                            color=dos_colours[ind],
                            label=label,
                        )
                    else:
                        label = None
                        if not options.get("unstacked_pdos"):
                            ax_dos.fill_between(
                                pdos_energies,
                                stack,
                                stack + pdos[projector],
                                lw=0,
                                alpha=alpha,
                                label=projector_labels[ind],
                                color=dos_colours[ind],
                            )
                            projector_outline_alpha = 0

                        ax_dos.plot(
                            pdos_energies,
                            stack + pdos[projector],
                            alpha=projector_outline_alpha,
                            color=dos_colours[ind],
                            label=label,
                        )

                    stacks[stack_key] += pdos[projector]

            if not options.get("pdos_hide_sum") and options.get("unstacked_pdos"):
                for stack_key in stacks:
                    if stack_key is None:
                        label = "Sum pDOS"
                    else:
                        label = "Sum pDOS: spin-{}".format(stack_key)
                    if options.get("plot_bandstructure"):
                        ax_dos.plot(
                            stacks[stack_key],
                            pdos_energies,
                            ls="--",
                            alpha=1,
                            color="black",
                            zorder=1e9,
                            label=label,
                        )
                    else:
                        ax_dos.plot(
                            pdos_energies,
                            stacks[stack_key],
                            ls="--",
                            alpha=1,
                            color="black",
                            zorder=1e9,
                            label=label,
                        )

    if len(seeds) == 1 and (plotting_pdos or "spin_dos" in dos_data):
        if options.get("plot_bandstructure"):
            dos_legend = ax_dos.legend(bbox_to_anchor=(1, 1))

        else:
            dos_legend = ax_dos.legend(bbox_to_anchor=(1, 0.5), loc="center left")

        bbox_extra_artists.append(dos_legend)

    elif len(seeds) > 1 and not plotting_pdos:
        if options.get("plot_bandstructure"):
            dos_legend = ax_dos.legend(bbox_to_anchor=(1, 1))
        else:
            dos_legend = ax_dos.legend(loc="upper right")

    return ax_dos


def projected_bandstructure_plot(
    dispersion,
    ax,
    path,
    bbox_extra_artists,
    eigs=None,
    pdis_interpolation_factor=2,
    pdis_point_scale=25,
    projectors_to_plot=None,
    **options,
):
    """Plot projected bandstructure with weightings from OptaDOS pdis.dat file.

    Parameters:
        dispersion (matador.orm.spectral.ElectronicDispersion): scraped
            data for bandstructure and pdis.
        seed (str): seed name for files to scrape.
        ax (matplotlib.pyplot.Axes): axis to plot on.
        bbox_extra_artists (list): list to append any legends too.

    Keyword arguments:
        eigs (np.ndarray): eigenvalues for the associated Dispesion object,
            passed separately to allow for reordering.
        interpolation_factor (float): amount by which to interpolate bands.
        point_scale (float): rescale points by this amount
        projectors_to_plot (list(tuple)): list of projectors to plot.

    Returns:
        matplotlib.pyplot.Axes: the axis that was plotted on.

    """

    if eigs is None:
        eigs = dispersion.eigs_s_k

    if projectors_to_plot is not None:
        if not any(
            projector in dispersion.projectors for projector in projectors_to_plot
        ):
            raise RuntimeError(
                "None of the desired projectors {} could be found in {}".format(
                    projectors_to_plot, dispersion.projectors
                )
            )

        _projectors_to_plot = []
        _projector_inds = []
        for ind, projector in enumerate(dispersion.projectors):
            if projector in projectors_to_plot:
                _projectors_to_plot.append(projector)
                _projector_inds.append(ind)

        pdis = np.zeros(
            (dispersion.num_kpoints, dispersion.num_bands, len(_projectors_to_plot))
        )
        for jnd, ind in enumerate(_projector_inds):
            pdis[:, :, jnd] = dispersion.projector_weights[:, :, ind]

        projectors = _projectors_to_plot

    else:
        pdis = np.array(dispersion.projector_weights, copy=True)
        projectors = copy.deepcopy(dispersion.projectors)

    pdis[pdis < 0] = 0
    pdis[pdis > 1] = 1

    keep_inds = []
    for ind, _ in enumerate(projectors):
        if np.max(pdis[:, :, ind]) > 1e-8:
            keep_inds.append(ind)

    projector_labels, dos_colours = _get_projector_info(
        projectors,
        colours_override=options["colour_cycle"] if options.get("colours") else None,
    )

    fermi_energy = options.get("external_efermi") or dispersion.fermi_energy

    _ordered_scatter(
        path,
        eigs[0].T - fermi_energy,
        pdis,
        dispersion.kpoint_branches,
        interpolation_factor=pdis_interpolation_factor,
        point_scale=pdis_point_scale,
        ax=ax,
        colours=dos_colours,
    )

    for ind, _ in enumerate(projectors):
        if ind in keep_inds:
            ax.scatter(
                1e20, 0, facecolor=dos_colours[ind], label=projector_labels[ind], lw=0
            )

    if not options.get("no_legend", False):
        legend = ax.legend(loc=1)
        legend.set_zorder(1e20)
        bbox_extra_artists.append(legend)

    return ax


def _ordered_scatter(
    path,
    eigs,
    pdis,
    branches,
    ax=None,
    colours=None,
    interpolation_factor=2,
    point_scale=25,
):
    """Plots an ordered scatter plot of a projected bandstructure.

    Parameters:
        path (np.ndarray): linearised [0, 1] kpoint path array.
        eigs (np.ndarray): (num_kpoints x num_bands) array containing eigenvalues
        pdis (np.ndarray): (num_kpoints x num_bands x num_projectors) array containing
            projector weights.
        branches (list): list of branch indices, e.g. for two branches [[0,1,2], [3, 4]].

    Keyword arguments:
        ax (matplotlib.Axes): axis to plot on
        colours (list): colours assigned for each projector.
        interpolation_factor (float): multiplier for fineness of band interpolation.

    """
    from scipy.interpolate import interp1d

    flat_pts_k = []
    flat_pts_e = []
    flat_sizes = []
    flat_colours = []
    flat_zorders = []

    for nb in range(len(eigs[0])):
        for branch_ind, branch in enumerate(branches):
            k = path[(np.asarray(branch) - branch_ind).tolist()]
            projections = pdis[branch, nb]
            ek_fn = interp1d(k, eigs[branch, nb])
            k_interp = np.linspace(
                np.min(k), np.max(k), num=int(interpolation_factor * len(k))
            )
            ek_interp = ek_fn(k_interp)
            projections = projections.T
            interp_projections = []
            for i, _ in enumerate(projections):
                interp_projections.append(interp1d(k, projections[i])(k_interp))
            projections = np.asarray(interp_projections).T
            pts = np.array([k_interp, ek_interp]).T.reshape(-1, 1, 2)

            if colours is not None:
                plot_colours = [colours[i] for i in range(len(projections[0]))]
            else:
                plot_colours = [None for i in range(len(projections[0]))]
            for i, _ in enumerate(projections):
                # use masked arrays to exclude the small projections
                sizes = np.ma.masked_where(
                    projections[i] <= PROJECTOR_MIN, np.cumsum(projections[i])
                )
                # zorders should be large and negative in order to pass rasterization condition on axis
                zorders = 1000 * (-100 * nb - sizes) - 1e7

                # this loop is slow, but will still be orders of magnitude faster than the matplotlib rendering
                for j in range(len(projections[i])):
                    flat_pts_k.append(pts[i, 0, 0])
                    flat_pts_e.append(pts[i, 0, 1])
                    size = sizes[j]
                    flat_sizes.append(point_scale * (size) ** 2)
                    flat_colours.append(plot_colours[j])
                    flat_zorders.append(zorders[j])

            # plot all bands in light grey as a skeleton
            ax.plot(pts[:, 0, 0], pts[:, 0, 1], lw=1, alpha=0.5, c="grey", zorder=0)

    flat_zorders = np.asarray(flat_zorders)
    flat_pts_k = np.asarray(flat_pts_k)[np.argsort(flat_zorders)]
    flat_pts_e = np.asarray(flat_pts_e)[np.argsort(flat_zorders)]
    flat_sizes = np.asarray(flat_sizes)[np.argsort(flat_zorders)]
    flat_colours = np.asarray(flat_colours)[np.argsort(flat_zorders)]

    ax.scatter(
        flat_pts_k,
        flat_pts_e,
        s=flat_sizes,
        c=flat_colours,
        marker="o",
        rasterized=True,
    )


def _get_lineprops(
    dispersion,
    spin_fermi_energy,
    nb,
    ns,
    branch,
    branch_ind,
    seed_ind,
    options,
    eigs=None,
):
    """Get the properties of the line to plot."""
    if seed_ind is None:
        seed_ind = 0
    colour = options.get("colour_cycle")[seed_ind]
    alpha = 1
    label = None

    if isinstance(dispersion, ElectronicDispersion):
        if eigs is None:
            eigs = dispersion.eigs
        if dispersion.num_spins == 2:
            if ns == 0:
                colour = "red"
                alpha = 0.8
            else:
                colour = "blue"
                alpha = 0.8

    if options.get("band_colour") is not None:
        colour = options.get("band_colour")

    if options.get("_mpl_cmap") is not None:
        colour = options["_mpl_cmap"][nb]

    if options.get("band_alpha") is not None:
        alpha = options["band_alpha"]

    if options.get("highlight_bands") is not None:
        if nb in options.get("highlight_bands"):
            colour = "red"
        else:
            alpha = 0.5

    if branch_ind == 0 and ns == 0 and nb == 0 and options.get("labels") is not None:
        label = options.get("labels")[seed_ind]

    return colour, alpha, label


def _add_path_labels(seed, dispersion, ax_dispersion, path, seed_ind, options):
    """Scrape k-point path labels from cell file and seekpath, then add them to the plot."""
    from matador.utils.cell_utils import doc2spg, get_seekpath_kpoint_path

    xticks = []
    xticklabels = []
    shear_planes = []
    labelled = []
    path_labels = dict()

    # first, try to grab them from the cell file
    if isinstance(seed, str) and os.path.isfile(seed + ".cell"):
        doc, success = cell2dict(
            seed + ".cell",
            db=False,
            verbosity=options.get("verbosity", 0),
            lattice=True,
            positions=True,
        )
        if not success:
            print(f"Unable to scrape {seed}.cell: {doc}")
            doc = {}
    else:
        doc = seed

    if options.get("phonons"):
        key = "phonon_fine_kpoint_path"
    else:
        key = "spectral_kpoints_path"
    if key in doc and key + "_labels" in doc:
        for label, point in zip(doc.get(key + "_labels", []), doc.get(key, None)):
            path_labels[label] = point
        print("Detected path labels from cell file")

    if not path_labels:
        # try to get dispersion path labels from spglib/seekpath
        spg_structure = None
        if isinstance(dispersion, Dispersion):
            try:
                spg_structure = doc2spg(dispersion)
            except (KeyError, RuntimeError) as exc:
                print(
                    f"Unable to create spglib structure from input data: skipping path labels: {exc}."
                )

        if not spg_structure:
            res = False
            cell = False
            if isinstance(seed, str):
                if os.path.isfile(seed + ".res"):
                    res = True
                elif os.path.isfile(seed + ".cell"):
                    cell = True
                else:
                    print(
                        "Failed to find {}.cell or {}.res, will not be able to generate labels.".format(
                            seed, seed
                        )
                    )

            success = False
            if cell:
                doc, success = cell2dict(
                    seed + ".cell",
                    db=False,
                    verbosity=options.get("verbosity", 0),
                    lattice=True,
                    positions=True,
                )
            if res and not success:
                doc, success = res2dict(
                    seed + ".res", db=False, verbosity=options.get("verbosity", 0)
                )

            if cell or res:
                if success:
                    spg_structure = doc2spg(doc)
                else:
                    print(
                        "Failed to scrape {}.cell/.res, will not be able to generate labels.".format(
                            seed
                        )
                    )

        if spg_structure:
            _, _, seekpath_results = get_seekpath_kpoint_path(
                spg_structure, standardize=False, explicit=False
            )
            path_labels = seekpath_results["point_coords"]

    for branch_ind, branch in enumerate(dispersion.kpoint_branches):
        for sub_ind, ind in enumerate(branch):
            kpt = dispersion.kpoint_path[ind]
            for label, point in path_labels.items():
                if np.allclose(point, kpt):
                    if ind - branch_ind not in labelled:
                        label = label.replace("GAMMA", r"\Gamma")
                        label = label.replace("SIGMA", r"\Sigma")
                        label = label.replace("DELTA", r"\Delta")
                        label = label.replace("LAMBDA", r"\Lambda")
                        if sub_ind == len(branch) - 1:
                            if branch_ind < len(dispersion.kpoint_branches) - 1:
                                _tmp = dispersion.kpoint_path
                                next_point = _tmp[
                                    dispersion.kpoint_branches[branch_ind + 1][0]
                                ]
                                for new_label, new_point in path_labels.items():
                                    new_label = new_label.replace("GAMMA", r"\Gamma")
                                    new_label = new_label.replace("SIGMA", r"\Sigma")
                                    new_label = new_label.replace("DELTA", r"\Delta")
                                    new_label = new_label.replace("LAMBDA", r"\Lambda")
                                    # import matplotlib
                                    if np.allclose(new_point, next_point):
                                        label = "\\dfrac{{{}}}{{{}}}".format(
                                            label, new_label
                                        )
                                        ax_dispersion.axvline(
                                            path[ind - branch_ind],
                                            ls="-",
                                            c="grey",
                                            zorder=1,
                                            lw=0.5,
                                        )
                                        labelled.append(ind - branch_ind)
                                        shear_planes.append(ind)
                        label = "${}$".format(label.replace("$", ""))
                        ax_dispersion.axvline(
                            path[ind - branch_ind], ls="--", c="grey", zorder=0, lw=0.5
                        )
                        xticklabels.append(label)
                        xticks.append(path[ind - branch_ind])
                        break

    if isinstance(dispersion, ElectronicDispersion) and options.get("gap"):
        if dispersion.num_spins != 1:
            raise NotImplementedError(
                "Band gap summary not implemented for multiple spin channels."
            )
        if dispersion.band_gap > 0:
            vbm_pos = dispersion["band_gap_path_inds"][1]
            vbm = dispersion["valence_band_min"] - dispersion.fermi_energy
            cbm_pos = dispersion["band_gap_path_inds"][0]
            cbm = dispersion["conduction_band_max"] - dispersion.fermi_energy
            if vbm_pos != cbm_pos:
                vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                ax_dispersion.plot(
                    [path[vbm_pos - vbm_offset], path[cbm_pos - cbm_offset]],
                    [vbm, cbm],
                    c="blue",
                    lw=5,
                    alpha=0.3,
                    zorder=0,
                    label="indirect gap {:3.3f} eV".format(cbm - vbm),
                )

            vbm_pos = dispersion["direct_gap_path_inds"][1]
            vbm = dispersion["direct_valence_band_min"] - dispersion.fermi_energy
            cbm_pos = dispersion["direct_gap_path_inds"][0]
            cbm = dispersion["direct_conduction_band_max"] - dispersion.fermi_energy
            vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
            cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
            ax_dispersion.plot(
                [path[vbm_pos - vbm_offset], path[cbm_pos - cbm_offset]],
                [vbm, cbm],
                c="red",
                lw=5,
                alpha=0.3,
                zorder=0,
                label="direct gap {:3.3f} eV".format(cbm - vbm),
            )
            if not options.get("no_legend", False):
                ax_dispersion.legend(
                    loc="upper center",
                    bbox_to_anchor=(0.5, 1.1),
                    fancybox=True,
                    shadow=True,
                    ncol=2,
                    handlelength=1,
                )

    if seed_ind == 0:
        ax_dispersion.set_xticks(xticks)
        ax_dispersion.set_xticklabels(xticklabels)
        ax_dispersion.grid(False)


def _get_projector_info(projectors, colours_override=None):
    """Grab appropriate colours and labels from a list of projectors.

    Parameters:
        projectors (list): list containing (element_str, l_channel) tuples.

    Returns:
        list: list of projector labels, e.g. {element_str}-${l_channel}$.
        list: list of colours for density of states, derived from vesta colours.

    """

    element_colours = get_element_colours()
    projector_labels = []
    dos_colours = {}

    all_species = set(proj[0] for proj in projectors)

    for ind, projector in enumerate(projectors):

        # pad out projectors for e.g. phonon case
        species = projector[0]
        if len(projector) > 1:
            ang_mom = projector[1]
        else:
            ang_mom = None
        if len(projector) > 2:
            spin = projector[2]
        else:
            spin = None

        # (species, None, None)
        if species is not None and ang_mom is None and spin is None:
            projector_label = species
        # (None, ang mom, None)
        if species is None and ang_mom is not None and spin is None:
            projector_label = "${}$".format(ang_mom)
        # (None, None, spin)
        elif species is None and ang_mom is None and spin is not None:
            projector_label = ""
        # (species, ang_mom, None/spin)
        elif species is not None and ang_mom is not None:
            projector_label = "{} (${}$)".format(species, ang_mom)
        # (species, None, None/spin)
        elif species is not None and ang_mom is None:
            projector_label = "{}".format(species)
        # (None, ang_mom, None/spin)
        elif species is None and ang_mom is not None:
            projector_label = "${}$".format(ang_mom)
        # (species, ang_mom, None/spin)
        else:
            projector_label = "{} (${}$)".format(species, ang_mom)

        projector_labels.append(projector_label)

        dos_colour = None

        # if species-projected only, then use VESTA colours
        if species is not None and ang_mom is None:
            dos_colour = element_colours.get(projector[0])
        # if species_ang-projected, then use VESTA colours but lightened
        elif len(all_species) > 1 and species is not None and ang_mom is not None:
            dos_colour = copy.deepcopy(element_colours.get(projector[0]))
            multi = ["s", "p", "d", "f"].index(projector[1]) - 1
            for jind, _ in enumerate(dos_colour):
                dos_colour[jind] = max(min(dos_colour[jind] + multi * 0.2, 1), 0)
        # otherwise if just ang-projected, use colour_cycle
        if dos_colour is None:
            dos_colour = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])[ind]

        if colours_override:
            try:
                dos_colour = colours_override[ind]
            except IndexError:
                pass

        # always want to use the same colour for spin projectors,
        # so collect them in a spin-agnostic way and then unfold into
        # list before returning
        spinless_projector = (species, ang_mom)
        if spinless_projector not in dos_colours:
            dos_colours[spinless_projector] = dos_colour

    dos_colours = [
        dos_colours[(species, ang_mom)] for (species, ang_mom, _) in projectors
    ]

    return projector_labels, dos_colours


def _load_electronic_dos(seed, options):
    """Try to obtain electronic DOS data, either from files, or as
    a dictionary.

    Parameters:
        seed (str/dict): either a filename or dictionary containing dos
            data.
        options (dict): plotting options.

    Returns:
        ElectronicDOS object containing scraped data.

    """
    if isinstance(seed, dict):
        return ElectronicDOS(seed)

    if isinstance(seed, ElectronicDOS):
        return seed

    seed = seed.replace(".bands", "")
    if options.get("dos") is None:
        # look for dat files, and just use the first
        # bands_dos should only be used as a last resort
        exts = ["adaptive.dat", "fixed.dat", "linear.dat", "pdos.dat", "bands_dos"]
        for ext in exts:
            if os.path.isfile("{}.{}".format(seed, ext)):
                dos_seed = "{}.{}".format(seed, ext)
                break
        else:
            raise RuntimeError("No total DOS files found.")
    else:
        dos_seed = options.get("dos")

    # If bands_dos exists, do some manual broadening:
    # .bands_dos is a file written by run3 when doing a
    # full spectral calculation, it is simply the .bands
    # file output from a DOS calculation
    if dos_seed.endswith(".bands_dos"):
        dos_data, s = bands2dict(dos_seed)
        gaussian_width = options.get("gaussian_width", 0.1)
        dos_data["dos"], dos_data["energies"] = DensityOfStates.bands_as_dos(
            dos_data, gaussian_width=gaussian_width
        )

        if isinstance(dos_data["dos"], dict):
            dos_data["spin_dos"] = dos_data["dos"]
            del dos_data["dos"]

    else:
        dos_data, s = optados2dict(dos_seed, verbosity=0)
        if dos_seed.endswith("pdos.dat"):
            _dos_data = {}
            _dos_data["pdos"] = dos_data
            dos_data = _dos_data

    # if a pdos.dat file is found, add it to the dos_data under the pdos key
    if not dos_seed.endswith("pdos.dat") and os.path.isfile(f"{seed}.pdos.dat"):
        pdos_data, s = optados2dict(f"{seed}.pdos.dat", verbosity=0)
        if not s:
            raise RuntimeError(pdos_data)
        dos_data["pdos"] = pdos_data

    if not s:
        raise RuntimeError(dos_data)

    return ElectronicDOS(dos_data)


def _load_phonon_dos(seed, options):
    """Try to obtain phonon DOS data, either from files, or as
    a dictionary.

    Parameters:
        seed (str/dict): either a filename or dictionary containing dos
            data.
        options (dict): plotting options.

    Returns:
        VibrationalDOS object containing scraped data.

    """

    if isinstance(seed, dict):
        return VibrationalDOS(seed)
    if isinstance(seed, DensityOfStates):
        return seed

    # otherwise, just read the phonon_dos file
    dos_data, s = phonon_dos2dict(seed + ".phonon_dos")
    if not s:
        raise RuntimeError(dos_data)

    return VibrationalDOS(dos_data)


def _parse_projectors_list(projectors):
    """Convert CLI args into the appropriate projector, ignoring
    spin channels.

    Parameters:
        projectors (str): a string of comma-separated element:orbital
            pairs. If the colon is omitted, all orbitals will be used.

    Returns:
        list(tuple): list of projectors in format [(element, orbital, spin)].

    """
    if projectors is None:
        return None
    _projectors = []
    orbitals = ["s", "p", "d", "f"]

    for projector in projectors.split(","):
        if ":" not in projector:
            element = projector
            for orbital in orbitals:
                _projectors.append((element, orbital, None))
            _projectors.append((element, None, None))
        else:
            element = projector.split(":")[0]
            orbital = projector.split(":")[1]
            _projectors.append((element, orbital, None))

    return _projectors
