# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot convergence test data
for a series of calculations.

"""

import glob
from decimal import Decimal, ROUND_UP
from collections import defaultdict
from traceback import print_exc
import numpy as np
from matador.plotting.plotting import plotting_function, SAVE_EXTS
from matador.scrapers.castep_scrapers import castep2dict
from matador.utils.chem_utils import get_formula_from_stoich

__all__ = ["plot_cutoff_kpt_grid"]


@plotting_function
def plot_cutoff_kpt_grid(
    data,
    figsize=None,
    forces=False,
    max_energy=None,
    max_force=None,
    legend=True,
    colour_by="formula",
    log=False,
    **kwargs,
):
    """Create a composite plot of cutoff energy and kpoint convergence.

    Parameters:
        data (dict): dictionary of convergence data.

    Keyword arguments:
        forces (bool): whether or not to plot forces.
        legend (bool): whether or not to show a legend.
        max_energy (float): max ylimit for convergence data in meV.

    """
    import matplotlib.pyplot as plt

    if figsize is None:
        if forces:
            figsize = (8, 7)
        else:
            figsize = (8, 4)

    if forces:
        _, axes = plt.subplots(2, 2, figsize=figsize, sharex="col", sharey="row")
    else:
        _, axes = plt.subplots(1, 2, figsize=figsize, sharey="row")

    if forces:
        ax = axes[0][0]
    else:
        ax = axes[0]

    k_label = "Reduced $k$-point spacing (1/Å)"
    e_label = "Planewave cutoff energy (eV)"

    _, lines, labels = plot_field(
        data,
        field="formation_energy_per_atom",
        parameter="cut_off_energy",
        ax=ax,
        colour_by=colour_by,
        max_y=max_energy,
        log=log,
        label_x=e_label if not forces else False,
        label_y=True,
        reference="last",
    )

    if len(labels) < 30 and legend:
        plt.figlegend(lines, labels, loc="upper center", ncol=min(len(lines) // 2, 8))

    if forces:
        ax = axes[0][1]
    else:
        ax = axes[1]

    plot_field(
        data,
        field="formation_energy_per_atom",
        parameter="kpoints_mp_spacing",
        ax=ax,
        colour_by=colour_by,
        max_y=max_energy,
        log=log,
        label_x=k_label if not forces else False,
        label_y=False,
        reference="last",
    )

    if forces:
        plot_field(
            data,
            field="forces",
            parameter="cut_off_energy",
            ax=axes[1][0],
            max_y=max_force,
            label_x=e_label,
            log=log,
            label_y=True,
            reference="last",
        )
        plot_field(
            data,
            field="forces",
            parameter="kpoints_mp_spacing",
            max_y=max_force,
            log=log,
            ax=axes[1][1],
            label_x=k_label,
            label_y=False,
            reference="last",
        )

    plt.tight_layout()

    if any(kwargs.get(ext) for ext in SAVE_EXTS):
        fname = "conv"
        for ext in SAVE_EXTS:
            if kwargs.get(ext):
                plt.savefig(
                    "{}.{}".format(fname, ext), bbox_inches="tight", transparent=False
                )
                print("Wrote {}.{}".format(fname, ext))

    if kwargs.get("plot_fname") or any([kwargs.get(ext) for ext in SAVE_EXTS]):
        import os

        fname = kwargs.get("plot_fname") or "conv"
        for ext in SAVE_EXTS:
            if kwargs.get(ext):
                fname_tmp = fname
                ind = 0
                while os.path.isfile("{}.{}".format(fname_tmp, ext)):
                    ind += 1
                    fname_tmp = fname + str(ind)

                fname = fname_tmp
                plt.savefig(
                    "{}.{}".format(fname, ext), bbox_inches="tight", transparent=True
                )
                print("Wrote {}.{}".format(fname, ext))

    if kwargs.get("show"):
        plt.show()


def plot_field(
    data,
    field="formation_energy_per_atom",
    parameter="cut_off_energy",
    ax=None,
    reference="last",
    max_y=None,
    log=False,
    colour_by="formula",
    label_x=True,
    label_y=True,
    plot_points=True,
):
    """Plot the convergence fields for each structure in `data` at
    each value of `parameter`.

    Parameters:
        data (dict): nested dictionary with keys for each structure
            that stores the convergence data to be plotted. e.g.
            {'structure_A':
                {'formation_energy_per_atom': [1, 2, 3],
                 'cutoff': [500, 400, 300]
                }
            }.

    Keyword arguments:
        field (str): the key under which convergence data is stored. This
            string is used for the y-axis label.
        parameter (str): the key under which the convergence parameter
            values are stored. This string is used for the x-axis label.
        ax (matplotlib.Axis): optional axis object to use for plot.
        reference (str): if 'last' then use the last value of the
            convergence data as the zero.
        plot_points (bool): whether or not to show markers

    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()

    lines = []
    labels = []
    labelled = []

    colourmap = {}
    if colour_by == "formula":
        formulae = sorted(
            list(
                {
                    get_formula_from_stoich(data[key]["stoichiometry"], tex=True)
                    for key in data
                }
            )
        )
        for formula in formulae:
            colourmap[formula] = next(ax._get_lines.prop_cycler)["color"]

    for ind, key in enumerate(data):

        formula = get_formula_from_stoich(data[key]["stoichiometry"], tex=True)
        try:
            values, parameters = get_convergence_values(
                data[key], parameter, field, reference=reference, log=log
            )
        except Exception:
            print("Missing data for {}->{}->{}".format(key, parameter, field))
            continue

        label = None
        if colour_by == "formula":
            c = colourmap[formula]
            if formula not in labelled:
                label = formula
                labelled.append(formula)
        else:
            c = next(ax._get_lines.prop_cycler)["color"]
            label = key

        # only plot points if there's less than 25 lines
        if plot_points:
            ax.plot(
                parameters,
                values,
                "o",
                markersize=5,
                alpha=1,
                label=label,
                lw=0,
                zorder=1000,
                c=c,
            )
            ax.plot(parameters, values, "-", alpha=0.2, lw=1, zorder=1000, c=c)
        else:
            ax.plot(
                parameters, values, "-", alpha=0.5, lw=1, zorder=1000, label=label, c=c
            )

        if label is not None:
            lines, labels = ax.get_legend_handles_labels()

    if label_x:
        if isinstance(label_x, bool):
            label_x = parameter.replace("_", " ")
        ax.set_xlabel(label_x)

    if label_y:
        if "force" in field.lower():
            ylabel = "$\\Delta F$ (eV/Å)"
        elif "formation" in field.lower():
            ylabel = "$\\Delta E$ (meV/atom)"
        else:
            ylabel = field.replace("_", " ")

        ax.set_ylabel(ylabel)

    if max_y is not None:
        ax.set_ylim(None, max_y)

    return ax, lines, labels


def round(n, prec):
    """Replace default (bankers) rounding with "normal" rounding."""
    if prec is None:
        return n
    else:
        return float(Decimal(str(n)).quantize(Decimal("0.05"), rounding=ROUND_UP))


def get_convergence_files(path, only=None):
    """Find all CASTEP files in the directory."""
    structure_files = defaultdict(list)
    files = glob.glob(path + "/*.castep")
    for file in files:
        if only is None or only in file:
            castep_dict, success = castep2dict(file, db=False)
            if not success:
                print("Failure to read castep file {}".format(file))
            else:
                source = castep_dict["source"][0].split("/")[-1]
                source = source.replace(".castep", "")
                source = "".join(source.split("_")[:-1])
                structure_files[source].append(castep_dict)

    return structure_files


def get_convergence_data(
    structure_files, conv_parameter="cut_off_energy", species=None
):
    """Parse cutoff energy/kpt spacing convergence calculations from list of files.

    Parameters:
        structure_files (list): list of filenames.

    Keyword arguments:
        conv_parameter (str): field for convergence parameter.
        species (list): only include structures with these species included.

    """
    scraped_from_filename = None

    form_set = {
        get_formula_from_stoich(
            structure_files[structure][0]["stoichiometry"], tex=True
        )
        for structure in structure_files
    }
    if len(form_set) == 1:
        chempot_mode = False
        single = True

    else:
        chempot_mode = True
        single = False
        chempots_dict = defaultdict(dict)

    if conv_parameter == "kpoints_mp_spacing":
        rounding = None
    else:
        rounding = None

    data = {}

    if chempot_mode:
        for key in structure_files:
            skip = False
            if species is not None:
                for elem in structure_files[key][0]["stoichiometry"]:
                    if elem[0] not in species:
                        skip = True
            if skip:
                continue
            for doc in structure_files[key]:
                scraped_from_filename = None
                if not single and len(doc["stoichiometry"]) == 1:
                    if conv_parameter == "kpoints_mp_spacing":
                        try:
                            scraped_from_filename = float(
                                doc["source"][0]
                                .split("/")[-1]
                                .split("_")[-1]
                                .split("A")[0]
                            )
                        except ValueError:
                            print(
                                f"Unable to determine kpoints label from {doc['source'][0]}, skipping..."
                            )
                            continue

                    if scraped_from_filename is not None:
                        rounded_field = round(scraped_from_filename, rounding)
                    else:
                        rounded_field = round(doc[conv_parameter], rounding)
                    energy_key = "total_energy_per_atom"
                    chempots_dict[str(rounded_field)][doc["atom_types"][0]] = doc[
                        energy_key
                    ]

    elems = set()

    if chempot_mode:
        for value in chempots_dict:
            for elem in chempots_dict[value]:
                elems.add(elem)

        for value in chempots_dict:
            for elem in elems:
                if elem not in chempots_dict[value]:
                    print(
                        "WARNING: {} chemical potential missing at {} = {} eV, skipping this value.".format(
                            elem, conv_parameter, value
                        )
                    )

    for key in structure_files:
        skip = False
        if species is not None:
            for elem in structure_files[key][0]["stoichiometry"]:
                if elem[0] not in species:
                    skip = True
        if skip:
            continue

        data[key] = {}
        data[key][conv_parameter] = defaultdict(list)

        for doc in structure_files[key]:
            if conv_parameter == "kpoints_mp_spacing":
                try:
                    scraped_from_filename = float(
                        doc["source"][0].split("/")[-1].split("_")[-1].split("A")[0]
                    )
                except ValueError:
                    print(
                        f"Unable to determine kpoints label from {doc['source'][0]}, skipping..."
                    )
                    continue
            try:
                doc["formation_energy_per_atom"] = doc["total_energy_per_atom"]
                if scraped_from_filename is not None:
                    rounded_field = round(scraped_from_filename, rounding)
                else:
                    rounded_field = round(doc[conv_parameter], rounding)
                if chempot_mode and len(doc["stoichiometry"]) > 1:
                    for atom in doc["atom_types"]:
                        doc["formation_energy_per_atom"] -= chempots_dict[
                            str(rounded_field)
                        ][atom] / len(doc["atom_types"])

                data[key][conv_parameter][conv_parameter].append(rounded_field)
                data[key]["stoichiometry"] = doc["stoichiometry"]
                data[key][conv_parameter]["formation_energy_per_atom"].append(
                    doc["formation_energy_per_atom"]
                )

                if "forces" in doc:
                    data[key][conv_parameter]["forces"].append(
                        np.sqrt(np.sum(np.asarray(doc["forces"]) ** 2, axis=-1))
                    )

            except Exception:
                print_exc()
                print("Error with {}".format(key))

    if conv_parameter == "kpoints_mp_spacing":
        reverse = True
    else:
        reverse = False

    for key in data:
        inds = [
            ind
            for ind, _ in sorted(
                enumerate(data[key][conv_parameter][conv_parameter]),
                key=lambda x: x[1],
                reverse=reverse,
            )
        ]
        for field in data[key][conv_parameter]:
            if field != "stoichiometry":
                data[key][conv_parameter][field] = [
                    data[key][conv_parameter][field][ind] for ind in inds
                ]

    return data


def get_convergence_values(data, parameter, field, reference="last", log=False):
    """Extract the data to plot for the given dictionary."""

    values = data[parameter][field]
    parameters = data[parameter][parameter]
    if len(values) != len(parameters):
        raise RuntimeError(
            "Mismatched convergence data for {}->{}: {}".format(field, parameter, data)
        )

    values = np.asarray(values)
    parameters = np.asarray(parameters)

    if reference == "last":
        values -= values[-1]
        if "force" not in field:
            values *= 1000
        values = np.abs(values)

    if log:
        with np.errstate(divide="ignore"):
            values = np.log10(values)

    return values, parameters


def combine_convergence_data(data_A, data_B):
    """Combine dictionaries with potentially overlapping keys."""
    data = {}
    for key in data_A:
        data[key] = data_A[key]
        if key in data_B:
            data[key].update(data_B[key])
    for key in data_B:
        if key not in data:
            data[key] = data_B[key]

    return data
