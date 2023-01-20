# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the scraper functions for CASTEP-related
inputs and outputs.

"""


from time import strptime
from collections import defaultdict
import os
import glob
import warnings
from pathlib import Path

import numpy as np
from matador.utils.cell_utils import (
    abc2cart,
    calc_mp_spacing,
    cart2volume,
    wrap_frac_coords,
    cart2abc,
    cart2frac,
)
from matador.utils.chem_utils import get_stoich, INVERSE_CM_TO_EV
from matador.scrapers.utils import (
    DFTError,
    ComputationError,
    scraper_function,
    f90_float_parse,
    get_flines_extension_agnostic,
)


@scraper_function
def res2dict(fname, db=True, **kwargs):
    """Extract available information from .res file; preferably
    used in conjunction with cell or param file.

    Parameters:
        fname (str or list): filename or list of filenames of res file(s)
            (with or without extension).

    Keyword arguments:
        db (bool): whether to fail if unable to scrape energies.

    Returns:
        (tuple): containing either dict/str containing data or error, and a bool stating
            if the scrape was successful.

    """

    flines, fname = get_flines_extension_agnostic(fname, "res")
    res = dict()

    # add .res to source
    res["source"] = [fname]
    # grab file owner username
    try:
        import pwd

        res["user"] = pwd.getpwuid(os.stat(fname).st_uid).pw_name
    except Exception:
        res["user"] = "xxx"

    try:
        get_seed_metadata(res, fname)
    except Exception as exc:
        warnings.warn(
            f"There was an error scraping provenance from filename {fname}: {exc}"
        )

    # alias special lines in res file
    titl = None
    cell = None
    remark = None
    for line in flines:
        if line.startswith("TITL") and db:
            # if not db, then don't read title
            titl = line.split()
            if len(titl) != 12:
                raise RuntimeError(f"Missing some TITL info in {fname}")
        elif line.startswith("CELL"):
            cell = line.split()
        elif line.startswith("REM"):
            remark = line.split()
        elif line.startswith("SFAC"):
            break

    if not cell:
        raise RuntimeError(f"Missing CELL line from {fname}")
    if not titl and db:
        raise RuntimeError(f"Missing TITL line in {fname}")
    if db:
        res["pressure"] = f90_float_parse(titl[2])
        res["cell_volume"] = f90_float_parse(titl[3])
        res["enthalpy"] = f90_float_parse(titl[4])
        res["num_atoms"] = int(titl[7])
        res["space_group"] = titl[8].strip("()")
        res["enthalpy_per_atom"] = res["enthalpy"] / res["num_atoms"]

    res["lattice_abc"] = [
        list(map(f90_float_parse, cell[2:5])),
        list(map(f90_float_parse, cell[5:8])),
    ]
    # calculate lattice_cart from abc
    res["lattice_cart"] = abc2cart(res["lattice_abc"])
    if "cell_volume" not in res:
        res["cell_volume"] = cart2volume(res["lattice_cart"])
    res["atom_types"] = []
    res["positions_frac"] = []
    res["site_occupancy"] = []
    for line_no, line in enumerate(flines):
        if "SFAC" in line:
            i = 1
            while "END" not in flines[line_no + i] and line_no + i < len(flines):
                # check if we don't have some other SHELX keyword in the way, e.g. "UNIT"
                if len(flines[line_no + i]) >= 4 and all(
                    [char.isupper() for char in flines[line_no + i][0:4]]
                ):
                    i += 1
                    continue
                cursor = flines[line_no + i].split()
                res["atom_types"].append(cursor[0])
                res["positions_frac"].append(list(map(f90_float_parse, cursor[2:5])))
                try:
                    res["site_occupancy"].append(f90_float_parse(cursor[5]))
                except IndexError:
                    res["site_occupancy"].append(1.0)
                i += 1

            break

    res["positions_frac"] = wrap_frac_coords(res["positions_frac"])
    res["num_atoms"] = len(res["atom_types"])
    # Parse any remark regarding implicit nanotube encapsulation
    if remark:
        if "NTPROPS" in remark:
            res["cnt_chiral"] = [0, 0]
            res["encapsulated"] = True
            for ind, entry in enumerate(remark):
                if "chiralN" in entry:
                    res["cnt_chiral"][0] = int(remark[ind + 1].replace(",", ""))
                if "chiralM" in entry:
                    res["cnt_chiral"][1] = int(remark[ind + 1].replace(",", ""))
                if entry == "'r':":
                    res["cnt_radius"] = f90_float_parse(
                        remark[ind + 1].replace(",", "")
                    )
                if entry == "'z':":
                    temp_length = remark[ind + 1].replace(",", "")
                    temp_length = temp_length.replace("\n", "")
                    temp_length = temp_length.replace("}", "")
                    res["cnt_length"] = f90_float_parse(temp_length)

    res["stoichiometry"] = get_stoich(res["atom_types"])
    res["num_fu"] = len(res["atom_types"]) / sum(
        [elem[1] for elem in res["stoichiometry"]]
    )

    return res, True


@scraper_function
def cell2dict(fname, db=False, lattice=True, positions=True, **kwargs):
    """Extract available information from .cell file; probably
    to be merged with another dict from a .param or .res file.

    Parameters:
        fname (str/list): filename or list of filenames of cell file(s)
            to scrape, with or without extension.

    Keyword arguments:
        db (bool): scrape database quality file
        lattice (bool): scrape lattice vectors
        positions (bool): scrape positions

    Returns:
        (tuple): containing either dict/str containing data or error, and a bool stating
            if the scrape was successful.

    """
    cell = dict()
    flines, fname = get_flines_extension_agnostic(fname, "cell")

    # add cell file to source
    cell["source"] = [fname]

    for line_no, line in enumerate(flines):
        if line.startswith(("#", "!")):
            continue
        if "#" or "!" in line:
            line = line.split("#")[0].split("!")[0]
        if "%block lattice_cart" in line.lower() and lattice:
            cell["lattice_cart"] = []
            i = 1
            while "endblock" not in flines[line_no + i].lower():
                if not flines[line_no + i].strip()[0].isalpha():
                    cell["lattice_cart"].append(
                        list(map(f90_float_parse, flines[line_no + i].split()))
                    )
                    if not len(cell["lattice_cart"][-1]) == 3:
                        raise RuntimeError(
                            "Lattice vector does not have enough elements!"
                        )
                i += 1
            if not len(cell["lattice_cart"]) == 3:
                raise RuntimeError("Wrong number of lattice vectors!")
        elif "%block lattice_abc" in line.lower() and lattice:
            cell["lattice_abc"] = []
            i = 1
            while "endblock" not in flines[line_no + i].lower():
                if not flines[line_no + i].strip()[0].isalpha():
                    cell["lattice_abc"].append(
                        list(map(f90_float_parse, flines[line_no + i].split()))
                    )
                    if not len(cell["lattice_abc"][-1]) == 3:
                        raise RuntimeError(
                            "Lattice vector does not have enough elements!"
                        )
                i += 1
            if not len(cell["lattice_abc"]) == 2:
                raise RuntimeError("Wrong specification of lattice_abc")

        elif "%block species_pot" in line.lower():
            cell["species_pot"] = dict()
            i = 1
            while "endblock" not in flines[line_no + i].lower():
                # handle blank lines in species pot
                split_line = flines[line_no + i].split()
                if not split_line:
                    i += 1
                    continue
                if len(split_line) == 2:
                    species = split_line[0]
                    pspot_string = split_line[1].split("/")[-1]
                    cell["species_pot"][species] = pspot_string.replace(
                        "()", ""
                    ).replace("[]", "")
                elif db:
                    raise RuntimeError(
                        f"Cannot parse `species_pot` block line {flines[line_no+i]} with `db=True`, "
                        "expected a (species, pspot) pair. "
                        "Try using `db=False` if specifying a pspot library for all species."
                    )
                elif len(split_line) == 1:
                    cell["species_pot"]["library"] = split_line[0].upper()

                i += 1
            if not cell["species_pot"]:
                cell.pop("species_pot")
        elif "%block cell_constraints" in line.lower():
            cell["cell_constraints"] = []
            for j in range(2):
                cell["cell_constraints"].append(
                    list(map(int, flines[line_no + j + 1].split()))
                )
            if (
                any(len(cell["cell_constraints"][i]) != 3 for i in range(2))
                or len(cell["cell_constraints"]) != 2
            ):
                raise RuntimeError("Invalid cell constraints block.")
        elif "%block hubbard_u" in line.lower():
            cell["hubbard_u"] = defaultdict(list)
            i = 0
            while "endblock" not in flines[line_no + i].lower():
                line = flines[line_no + i]
                if line == "eV" or len(line.split()) < 3:
                    i += 1
                    continue
                else:
                    atom = line.split()[0]
                    orbital = line.split()[1].replace(":", "")
                    shift = f90_float_parse(line.split()[-1])
                    atom = line.split()[0]
                    cell["hubbard_u"][atom] = dict()
                    cell["hubbard_u"][atom][orbital] = shift
                    i += 1
        elif "%block external_pressure" in line.lower():
            cell["external_pressure"] = np.zeros((3, 3))
            i = 1
            j = 0
            while "endblock" not in flines[line_no + i].lower():
                if not flines[line_no + i].strip()[0].isalpha():
                    flines[line_no + i] = flines[line_no + i].replace(",", "")
                    vals = list(map(f90_float_parse, flines[line_no + i].split()))
                    if len(vals) != (3 - j):
                        raise RuntimeError(
                            "External pressure should be specified as upper triangular matrix."
                        )
                    cell["external_pressure"][j] = np.asarray(j * [0.0] + vals).reshape(
                        3
                    )
                    j += 1
                i += 1
            cell["external_pressure"] = cell["external_pressure"].tolist()

        elif "%block external_efield" in line.lower():
            cell["external_efield"] = [
                f90_float_parse(e) for e in flines[line_no + 1].split()
            ]
            if len(cell["external_efield"]) != 3:
                raise RuntimeError(
                    f"EXTERNAL_EFIELD block has wrong shape, should be 3-D not: {cell['external_efield']}"
                )

        elif "%block ionic_constraints" in line.lower():
            cell["ionic_constraints"] = []
            i = 1
            while "endblock" not in flines[line_no + i].lower():
                cell["ionic_constraints"].append(flines[line_no + i].strip())
                i += 1
        # parse kpoints
        elif (
            "kpoints_mp_spacing" in line.lower() or "kpoint_mp_spacing" in line.lower()
        ):
            if (
                "spectral_kpoints_mp_spacing" in line.lower()
                or "spectral_kpoint_mp_spacing" in line.lower()
            ):
                cell["spectral_kpoints_mp_spacing"] = f90_float_parse(line.split()[-1])
            elif (
                "bs_kpoints_mp_spacing" in line.lower()
                or "bs_kpoint_mp_spacing" in line.lower()
            ):
                cell["spectral_kpoints_mp_spacing"] = f90_float_parse(line.split()[-1])
            elif (
                "supercell_kpoints_mp_spacing" in line.lower()
                or "supercell_kpoint_mp_spacing" in line.lower()
            ):
                cell["supercell_kpoints_mp_spacing"] = f90_float_parse(line.split()[-1])
            elif (
                "phonon_kpoints_mp_spacing" in line.lower()
                or "phonon_kpoint_mp_spacing" in line.lower()
            ):
                cell["phonon_kpoint_mp_spacing"] = f90_float_parse(line.split()[-1])
            elif (
                "phonon_fine_kpoints_mp_spacing" in line.lower()
                or "phonon_fine_kpoint_mp_spacing" in line.lower()
            ):
                cell["phonon_fine_kpoint_mp_spacing"] = f90_float_parse(
                    line.split()[-1]
                )
            else:
                cell["kpoints_mp_spacing"] = f90_float_parse(line.split()[-1])
        elif "kpoints_mp_grid" in line.lower() or "kpoint_mp_grid" in line.lower():
            if (
                "spectral_kpoints_mp_grid" in line.lower()
                or "spectral_kpoint_mp_grid" in line.lower()
            ):
                cell["spectral_kpoints_mp_grid"] = list(map(int, line.split()[-3:]))
            if (
                "bs_kpoints_mp_grid" in line.lower()
                or "bs_kpoint_mp_grid" in line.lower()
            ):
                cell["spectral_kpoints_mp_grid"] = list(map(int, line.split()[-3:]))
            # these two keywords do not have a corresponding "kpoints" alias (instead of kpoint)
            # so must remain unpluralised (see below for other phonon keyword exceptions)
            elif (
                "phonon_kpoints_mp_grid" in line.lower()
                or "phonon_kpoint_mp_grid" in line.lower()
            ):
                cell["phonon_kpoint_mp_grid"] = list(map(int, line.split()[-3:]))
            elif (
                "phonon_fine_kpoints_mp_grid" in line.lower()
                or "phonon_fine_kpoint_mp_grid" in line.lower()
            ):
                cell["phonon_fine_kpoint_mp_grid"] = list(map(int, line.split()[-3:]))
            else:
                cell["kpoints_mp_grid"] = list(map(int, line.split()[-3:]))
        elif "kpoints_mp_offset" in line.lower() or "kpoint_mp_offset" in line.lower():
            if (
                "spectral_kpoints_mp_offset" in line.lower()
                or "spectral_kpoint_mp_offset" in line.lower()
            ):
                cell["spectral_kpoints_mp_offset"] = list(
                    map(f90_float_parse, line.split()[-3:])
                )
            if (
                "bs_kpoints_mp_offset" in line.lower()
                or "bs_kpoint_mp_offset" in line.lower()
            ):
                cell["spectral_kpoints_mp_offset"] = list(
                    map(f90_float_parse, line.split()[-3:])
                )
            elif (
                "phonon_kpoints_mp_offset" in line.lower()
                or "phonon_kpoint_mp_offset" in line.lower()
            ):
                # this is a special case where phonon_kpointS_mp_offset doesn't exist
                cell["phonon_kpoint_mp_offset"] = list(
                    map(f90_float_parse, line.split()[-3:])
                )
            elif (
                "phonon_fine_kpoints_mp_offset" in line.lower()
                or "phonon_fine_kpoint_mp_offset" in line.lower()
            ):
                cell["phonon_fine_kpoint_mp_offset"] = list(
                    map(f90_float_parse, line.split()[-3:])
                )
            else:
                cell["kpoints_mp_offset"] = list(
                    map(f90_float_parse, line.split()[-3:])
                )
        elif "%block" in line.lower() and (
            "_kpoints_path" in line.lower() or "_kpoint_path" in line.lower()
        ):
            if "spectral_" in line.lower() or "bs_" in line.lower():
                key = "spectral_kpoints_path"
            elif "phonon_fine_" in line.lower():
                key = "phonon_fine_kpoint_path"
            else:
                raise RuntimeError(f"Found unknown kpoint path key in line: {line}.")
            labels = []
            found_labels = False
            comment_delims = ["!", "%", "#"]
            cell[key] = []
            for line in flines[line_no + 1 :]:
                if "%endblock" in line.lower():
                    break

                # ignore CASTEP BREAK keyword and let dispersion script figure out discontinuities
                if "BREAK" in line.strip().upper():
                    continue

                cell[key].append(list(map(f90_float_parse, line.split()[:3])))

                for delim in comment_delims:
                    if delim in line:
                        labels.append(line.split(delim)[-1].strip())
                        found_labels = True
                        break
                else:
                    labels.append(None)

            else:
                raise RuntimeError(
                    f"Unable to find closing of block {key}, missing %endblock."
                )

            if found_labels:
                cell[key + "_labels"] = labels

        elif any(
            block in line.lower()
            for block in [
                "%block spectral_kpoints_list",
                "%block spectral_kpoint_list",
                "%block bs_kpoint_list",
                "%block bs_kpoints_list",
            ]
        ):
            i = 1
            cell["spectral_kpoints_list"] = []
            while "%endblock" not in flines[line_no + i].lower():
                cell["spectral_kpoints_list"].append(
                    list(map(f90_float_parse, flines[line_no + i].split()[:4]))
                )
                i += 1
        elif (
            "%block phonon_fine_kpoints_list" in line.lower()
            or "%block phonon_fine_kpoint_list" in line.lower()
        ):
            i = 1
            # this is a special case where phonon_fine_kpointS_list doesn't exist
            cell["phonon_fine_kpoint_list"] = []
            while "%endblock" not in flines[line_no + i].lower():
                cell["phonon_fine_kpoint_list"].append(
                    list(map(f90_float_parse, flines[line_no + i].split()[:4]))
                )
                i += 1
        elif "%block phonon_supercell_matrix" in line.lower():
            cell["phonon_supercell_matrix"] = []
            i = 1
            while "endblock" not in flines[line_no + i].lower():
                cell["phonon_supercell_matrix"].append(
                    list(map(int, flines[line_no + i].split()))
                )
                if not len(cell["phonon_supercell_matrix"][-1]) == 3:
                    raise RuntimeError(
                        "Supercell matrix row does not have enough elements!"
                    )
                i += 1
            if not len(cell["phonon_supercell_matrix"]) == 3:
                raise RuntimeError("Wrong supercell matrix shape!")

        elif not db:
            if "%block positions_frac" in line.lower():
                atomic_init_spins = []
                i = 1
                if positions:
                    cell["atom_types"] = []
                    cell["positions_frac"] = []
                while "%endblock positions_frac" not in flines[line_no + i].lower():
                    line = flines[line_no + i].split()
                    if positions:
                        cell["atom_types"].append(line[0])
                        cell["positions_frac"].append(
                            list(map(f90_float_parse, line[1:4]))
                        )
                    if "spin=" in flines[line_no + i].lower():
                        split_line = flines[line_no + i].split()
                        atomic_init_spins.append(
                            float(split_line[-1].lower().replace("spin=", ""))
                        )
                    else:
                        atomic_init_spins.append(None)
                    i += 1
                if any(atomic_init_spins):
                    cell["atomic_init_spins"] = atomic_init_spins
                    if len(cell["atomic_init_spins"]) != len(cell["positions_frac"]):
                        raise RuntimeError("Atomic init spins do not match positions")
                if positions:
                    cell["positions_frac"] = wrap_frac_coords(cell["positions_frac"])
                    cell["num_atoms"] = len(cell["atom_types"])
            elif "%block positions_abs" in line.lower():
                atomic_init_spins = []
                i = 1
                # avoid units
                if len(flines[line_no + i].split()) < 3:
                    i += 1
                if positions:
                    cell["atom_types"] = []
                    cell["positions_abs"] = []
                while "%endblock positions_abs" not in flines[line_no + i].lower():
                    line = flines[line_no + i].split()
                    if positions:
                        cell["atom_types"].append(line[0])
                        cell["positions_abs"].append(
                            list(map(f90_float_parse, line[1:4]))
                        )
                    if "spin=" in flines[line_no + i].lower():
                        split_line = flines[line_no + i].split()
                        atomic_init_spins.append(
                            float(split_line[-1].lower().replace("spin=", ""))
                        )
                    else:
                        atomic_init_spins.append(None)
                    i += 1
                if any(atomic_init_spins):
                    cell["atomic_init_spins"] = atomic_init_spins
                    if len(cell["atomic_init_spins"]) != len(cell["positions_frac"]):
                        raise RuntimeError("Atomic init spins do not match positions")

            elif "fix_com" in line.lower():
                cell["fix_com"] = bool(line.split()[-1])
            elif "fix_all_ions" in line.lower():
                cell["fix_all_ions"] = bool(line.split()[-1])
            elif "fix_all_cell" in line.lower():
                cell["fix_all_cell"] = bool(line.split()[-1])
            elif "fix_vol" in line.lower():
                cell["fix_vol"] = bool(line.split()[-1])
            elif "symmetry_generate" in line.lower():
                cell["symmetry_generate"] = True
            elif "symmetry_tol" in line.lower():
                cell["symmetry_tol"] = f90_float_parse(line.split()[-1])
            elif "snap_to_symmetry" in line.lower():
                cell["snap_to_symmetry"] = True
            elif "quantisation_axis" in line.lower():
                cell["quantisation_axis"] = list(map(int, line.split()[1:]))
            elif "positions_noise" in line.lower():
                cell["positions_noise"] = f90_float_parse(line.split()[-1])
            elif "cell_noise" in line.lower():
                cell["cell_noise"] = f90_float_parse(line.split()[-1])
            elif (
                "kpoints_path" in line.lower()
                or "kpoint_path" in line.lower()
                and "%block" not in line.lower()
            ):
                if (
                    "spectral_kpoints_path_spacing" in line.lower()
                    or "spectral_kpoint_path_spacing" in line.lower()
                ):
                    cell["spectral_kpoints_path_spacing"] = f90_float_parse(
                        line.split()[-1]
                    )
                if (
                    "bs_kpoints_path_spacing" in line.lower()
                    or "bs_kpoint_path_spacing" in line.lower()
                ):
                    cell["spectral_kpoints_path_spacing"] = f90_float_parse(
                        line.split()[-1]
                    )
                elif (
                    "phonon_fine_kpoints_path_spacing" in line.lower()
                    or "phonon_fine_kpoint_path_spacing" in line.lower()
                ):
                    cell["phonon_fine_kpoint_path_spacing"] = f90_float_parse(
                        line.split()[-1]
                    )
                elif (
                    "kpoints_path_spacing" in line.lower()
                    or "kpoint_path_spacing" in line.lower()
                ):
                    cell["kpoints_path_spacing"] = f90_float_parse(line.split()[-1])

    if "external_pressure" not in cell or not cell["external_pressure"]:
        cell["external_pressure"] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    if lattice:
        if "lattice_cart" not in cell and "lattice_abc" in cell:
            cell["lattice_cart"] = abc2cart(cell["lattice_abc"])
        elif "lattice_cart" in cell and "lattice_abc" not in cell:
            cell["lattice_abc"] = cart2abc(cell["lattice_cart"])
        if "lattice_cart" in cell:
            cell["cell_volume"] = cart2volume(cell["lattice_cart"])

    if positions:
        if "positions_frac" not in cell and "positions_abs" in cell:
            cell["positions_frac"] = cart2frac(
                cell["lattice_cart"], cell["positions_abs"]
            )
            cell["positions_frac"] = wrap_frac_coords(cell["positions_frac"])

    if db:
        for species in cell["species_pot"]:
            if "OTF" in cell["species_pot"][species].upper():
                pspot_seed = (
                    "/".join(fname.split("/")[:-1]) + "/" + cell["species_pot"][species]
                )
                if os.path.isfile(pspot_seed):
                    cell["species_pot"].update(usp2dict(pspot_seed))

    return cell, True


@scraper_function
def param2dict(fname, db=True, **kwargs):
    """Extract available information from .param file; probably
    to be merged with other dicts from other files.

    Parameters:
        fname (str/list): param filename or list of filenames with or
            without file extension

    Keyword arguments:
        db (bool): if True, only scrape relevant info, otherwise scrape all

    Returns:
        (tuple): containing either dict/str containing data or error, and a bool stating
            if the scrape was successful.

    """
    from matador.utils.castep_params import CASTEP_PARAMS

    flines, fname = get_flines_extension_agnostic(fname, "param")
    param = {}
    param["source"] = [fname]
    # exclude some useless info if importing to db
    scrub_list = [
        "checkpoint",
        "write_bib",
        "mix_history_length",
        "fix_occupancy",
        "page_wvfns",
        "num_dump_cycles",
        "backup_interval",
        "geom_max_iter",
        "fixed_npw",
        "write_cell_structure",
        "bs_write_eigenvalues",
        "calculate_stress",
        "opt_strategy",
        "max_scf_cycles",
    ]
    splitters = [":", "=", "\t", " "]
    unrecognised = []
    devel_lines = []
    for line_no, line in enumerate(flines):
        if "#" or "!" in line:
            line = line.split("#")[0].split("!")[0]
        line = line.lower()
        # skip blank lines and comments
        if line.startswith(("#", "!")) or not line.strip() or line_no in devel_lines:
            continue
        else:
            # if scraping to db, ignore "rubbish"
            if db:
                if [rubbish for rubbish in scrub_list if rubbish in line]:
                    continue
            # read all other parameters in
            for splitter in splitters:
                if splitter in line:
                    if not line.startswith("%"):
                        keyword = line.split(splitter)[0].strip()
                        value = line.split(splitter)[-1].strip()
                        if keyword.lower() not in CASTEP_PARAMS:
                            unrecognised.append(keyword.lower())
                        param[keyword] = value
                    # deal with edge cases
                    if "%block devel_code" in line:
                        i = 1
                        while "%endblock devel_code" not in flines[line_no + i].lower():
                            if i + line_no >= len(flines):
                                raise RuntimeError("Found unclosed %block devel_code.")
                            if "magres" in flines[line_no + i].lower():
                                line = flines[line_no + i].upper()
                            else:
                                line = flines[line_no + i].lower()
                            devel_lines.append(line_no + i)
                            if "devel_code" not in param:
                                param["devel_code"] = ""
                            param["devel_code"] += line
                            i += 1
                        break

                    if "true" in value.lower():
                        param[keyword] = True
                    elif "false" in value.lower():
                        param[keyword] = False

                    if "spin_polarised" in line:
                        param["spin_polarized"] = param["spin_polarised"]
                        if "spin_polarised" in param:
                            del param["spin_polarised"]

                    if keyword == "geom_max_iter":
                        param["geom_max_iter"] = int(param["geom_max_iter"])

                    if "cut_off_energy" in line and "mix_cut_off_energy" not in line:
                        temp_cut_off = param["cut_off_energy"].split()
                        if len(temp_cut_off) > 1:
                            if temp_cut_off[1] == "ev":
                                param["cut_off_energy"] = f90_float_parse(
                                    temp_cut_off[0]
                                )
                            elif db:
                                raise RuntimeError(
                                    "cut_off_energy units must be eV or blank in db mode, not {}".format(
                                        temp_cut_off[1]
                                    )
                                )
                            else:
                                param["cut_off_energy"] = "{} {}".format(
                                    temp_cut_off[0], temp_cut_off[1]
                                )
                        else:
                            param["cut_off_energy"] = f90_float_parse(temp_cut_off[0])
                    # standardize tasks that have alternative spellings
                    elif (
                        "task" in line
                        and "spectral_task" not in line
                        and "magres_task" not in line
                    ):
                        if "singlepoint" in param["task"]:
                            param["task"] = "singlepoint"
                        elif "geometry" in param["task"]:
                            param["task"] = "geometryoptimization"
                    elif "xc_functional" in line:
                        param["xc_functional"] = param["xc_functional"].upper()
                    elif "spin" in line and "polar" not in line and "spin" in param:
                        param["spin"] = int(param["spin"])
                    elif "perc_extra_bands" in line:
                        param["perc_extra_bands"] = f90_float_parse(
                            param["perc_extra_bands"]
                        )
                    elif "geom_force_tol" in line:
                        param["geom_force_tol"] = f90_float_parse(
                            param["geom_force_tol"]
                        )
                    elif "elec_energy_tol" in line:
                        temp = param["elec_energy_tol"].lower().replace("ev", "")
                        temp = temp.strip()
                        param["elec_energy_tol"] = f90_float_parse(temp)
                    break

    if unrecognised:
        raise RuntimeError(
            "Found several unrecognised parameters: {}".format(unrecognised)
        )

    return param, True


@scraper_function
def castep2dict(fname, db=True, intermediates=False, timings=False, **kwargs):
    """From seed filename, create dict of the most relevant
    information about a calculation.

    Parameters:
        fname (str/list): filename or list of filenames of castep file(s)

    Keyword arguments:
        db (bool): whether to error on missing relaxation info
        intermediates (bool): instead of a single dict containing the relaxed structure
            return a list of snapshots found in .castep file
        timings (bool): Run through the CASTEP file one extra time to calculate total time
            taken.

    Returns:
        (tuple): containing either dict/str containing data or error, and a bool stating
            if the scrape was successful.

    """
    flines, fname = get_flines_extension_agnostic(
        fname, ["castep", "history", "history.gz"]
    )

    castep = dict()
    # set source tag to castep file
    castep["source"] = [fname]

    # grab file owner
    try:
        import pwd

        castep["user"] = pwd.getpwuid(os.stat(fname).st_uid).pw_name
    except Exception:
        castep["user"] = "xxx"

    try:
        get_seed_metadata(castep, fname)
    except Exception as exc:
        warnings.warn(
            f"There was an error scraping provenance from filename {fname}: {exc}"
        )

    # wrangle castep file for parameters in 3 passes:
    # once forwards to get number and types of atoms
    _castep_scrape_atoms(flines, castep)
    # once backwards to get the final parameter set for the calculation
    _castep_scrape_final_parameters(flines, castep)

    # task specific options
    if db and "geometry" not in castep["task"]:
        raise RuntimeError("CASTEP file does not contain GO calculation")

    if not db and "thermo" in castep["task"].lower():
        _castep_scrape_thermo_data(flines, castep)

    if (
        not db
        and "thermo" in castep["task"].lower()
        or "phonon" in castep["task"].lower()
    ):
        _castep_scrape_phonon_frequencies(flines, castep)

    # only scrape snapshots/number of intermediates if requested
    try:
        snapshots, castep["geom_iter"] = _castep_scrape_all_snapshots(
            flines, intermediates=intermediates
        )
        if intermediates:
            castep["intermediates"] = snapshots
    except RuntimeError as exc:
        if intermediates:
            raise RuntimeError("Failed to scrape intermediates: {}".format(exc))

    if timings:
        _castep_scrape_metadata(flines, castep)

    # once more forwards, from the final step, to get the final structure
    _castep_scrape_final_structure(flines, castep, db=db)

    # scrape any BEEF post-processing
    _castep_scrape_beef(flines, castep)

    # scrape any AJM group-specific devel codes
    _castep_scrape_devel_code(flines, castep)

    if "positions_frac" not in castep or not castep["positions_frac"]:
        raise ComputationError("Could not find positions")

    # unfortunately CASTEP does not write forces when there is only one atom
    if (
        "forces" not in castep
        and castep["num_atoms"] == 1
        and "geometry" in castep["task"]
    ):
        castep["forces"] = [[0, 0, 0]]

    # finally check for pseudopotential files if OTF is present in species_pot
    if db:
        for species in castep["species_pot"]:
            if "OTF" in castep["species_pot"][species].upper():
                pspot_seed = str(
                    Path(fname).parent.joinpath(castep["species_pot"][species])
                )
                for globbed in glob.glob(pspot_seed):
                    if os.path.isfile(globbed):
                        castep["species_pot"].update(usp2dict(globbed))

    # check that any optimized results were saved and raise errors if not
    if not castep.get("optimised"):
        castep["optimised"] = False
        if db:
            # if importing to db, skip unconverged structure
            raise DFTError("CASTEP GO failed to converge.")

    return castep, True


@scraper_function
def bands2dict(fname, **kwargs):
    """Parse a CASTEP bands file into a dictionary, which can be used as input to
    an :obj:`matador.orm.spectral.ElectronicDispersion` object.

    Parameters:
        fname (str/list): filename of list of filenames to be scraped.

    Returns:
        (tuple): containing either dict/str containing data or error, and a bool stating
            if the scrape was successful.

    """
    from matador.utils.chem_utils import HARTREE_TO_EV, BOHR_TO_ANGSTROM

    bs = dict()
    flines, fname = get_flines_extension_agnostic(fname, "bands")

    header = flines[:9]
    data = flines[9:]

    bs["source"] = [fname]

    bs["num_kpoints"] = int(header[0].split()[-1])
    bs["num_spins"] = int(header[1].split()[-1])
    bs["num_electrons"] = f90_float_parse(header[2].split()[-1])
    bs["num_bands"] = int(header[3].split()[-1])

    bs["spin_fermi_energy_Ha"] = [
        f90_float_parse(val) for val in header[4].split()[-bs["num_spins"] :]
    ]
    bs["spin_fermi_energy"] = [
        HARTREE_TO_EV * val for val in bs["spin_fermi_energy_Ha"]
    ]
    bs["fermi_energy"] = np.mean(bs["spin_fermi_energy"])

    bs["lattice_cart"] = []
    for i in range(3):
        bs["lattice_cart"].append(
            [BOHR_TO_ANGSTROM * f90_float_parse(elem) for elem in header[6 + i].split()]
        )
    bs["kpoint_path"] = np.zeros((bs["num_kpoints"], 3))
    bs["kpoint_weights"] = np.zeros((bs["num_kpoints"]))
    bs["eigenvalues_k_s"] = np.empty(
        (bs["num_spins"], bs["num_bands"], bs["num_kpoints"])
    )

    for nk in range(bs["num_kpoints"]):
        kpt_ind = nk * (bs["num_spins"] * bs["num_bands"] + bs["num_spins"] + 1)
        bs["kpoint_path"][int(data[kpt_ind].split()[1]) - 1] = np.asarray(
            [f90_float_parse(elem) for elem in data[kpt_ind].split()[-4:-1]]
        )
        bs["kpoint_weights"][int(data[kpt_ind].split()[1]) - 1] = f90_float_parse(
            data[kpt_ind].split()[-1]
        )
        for nb in range(bs["num_bands"]):
            for ns in range(bs["num_spins"]):
                line_number = kpt_ind + 2 + ns + (ns * bs["num_bands"]) + nb
                bs["eigenvalues_k_s"][ns][nb][
                    int(data[kpt_ind].split()[1]) - 1
                ] = f90_float_parse(data[line_number].strip())

    bs["eigenvalues_k_s"] *= HARTREE_TO_EV
    bs["eigs_s_k"] = bs["eigenvalues_k_s"]

    return bs, True


@scraper_function
def arbitrary2dict(fname, **kwargs):
    """Read arbitrary CASTEP-style input files into
    a dictionary.

    Parameters:
        fname (str/list): filename or list of filenames.

    Returns:
        (tuple): containing either dict/str containing data or error, and a bool stating
            if the scrape was successful.

    """

    flines, fname = get_flines_extension_agnostic(fname, None)

    splitters = [":", "=", "\t", " "]
    comment_delims = ["!", "%", "#", "/*", "*/"]

    result = {}
    result["source"] = [fname]
    for line in flines:
        comment = False
        for delim in comment_delims:
            if line.strip().startswith(delim):
                comment = True
                break
        if not comment:
            for splitter in splitters:
                if splitter in line:
                    keyword = line.split(splitter)[0].strip()
                    value = line.split(splitter)[-1].strip()
                    result[keyword.lower()] = value
                    break

    return result, True


@scraper_function
def optados2dict(fname, **kwargs):
    """Scrape optados output file (*.*.dat) or (*.pdos.*.dat)
    for DOS, projectors and projected DOS/dispersion.

    Parameters:
        fname (str/list): optados filename or list of filenames.

    Returns:
        (tuple): containing either dict/str containing data or error, and a bool stating
            if the scrape was successful.

    """
    optados = dict()
    is_pdos = False
    is_pdis = False
    is_spin_dos = False
    dos_unit_label = None

    flines, fname = get_flines_extension_agnostic(fname, None)

    header = []
    header_ind = 0
    for ind, line in enumerate(flines):
        header_ind = ind
        if not line.strip().startswith("#") or "K-point" in line:
            break
        if "Partial" in line:
            is_pdos = True
        elif "spin" in line:
            is_spin_dos = True
        elif "Projected Dispersion" in line:
            is_pdis = True
        elif "DOS (" in line and "Integrated" not in line:
            dos_unit_label = " ".join(line.split()[2:])
        else:
            header.append(line)

    flines = flines[header_ind:]

    if not is_pdis:
        data = np.loadtxt(fname, comments="#")
        optados["energies"] = data[:, 0]

    if is_pdos or is_pdis:
        optados["projectors"] = _optados_get_projector_labels(header)
        optados["num_projectors"] = len(optados["projectors"])

    if dos_unit_label is not None:
        optados["dos_unit_label"] = dos_unit_label

    if is_pdos:
        # get pdos values
        optados["pdos"] = dict()
        optados["sum_pdos"] = np.zeros_like(data[:, 0])
        for i, projector in enumerate(optados["projectors"]):
            # optados spin-down projectors are negative, unfortunately
            if "down" in projector:
                optados["pdos"][projector] = -data[:, i + 1]
            else:
                optados["pdos"][projector] = data[:, i + 1]
            optados["sum_pdos"] += data[:, i + 1]

    elif is_spin_dos:
        optados["spin_dos"] = dict()
        optados["spin_dos"]["up"] = data[:, 1]
        optados["spin_dos"]["down"] = -data[:, 2]
        optados["dos"] = np.abs(optados["spin_dos"]["up"]) + np.abs(
            optados["spin_dos"]["down"]
        )

    elif is_pdis:
        optados["projector_weights"] = []
        optados["kpoints"] = []
        optados["eigenvalues"] = []
        # get kpoints and count number of bands
        kpt_ind = -1
        for i, line in enumerate(flines):
            if "K-point" in line:
                optados["kpoints"].append(
                    [f90_float_parse(val) for val in line.split()[-3:]]
                )
                if kpt_ind == -1:
                    kpt_ind = i
                else:
                    if not optados.get("num_bands"):
                        optados["num_bands"] = i - kpt_ind - 1

        optados["num_kpoints"] = len(optados["kpoints"])

        for nk in range(optados["num_kpoints"]):
            eigs = []
            pdis = []
            for nb in range(0, optados["num_bands"]):
                eigs.append(
                    f90_float_parse(
                        flines[nk * (optados["num_bands"] + 1) + nb + 1].split()[0]
                    )
                )
                pdis.append(
                    [
                        f90_float_parse(val)
                        for val in flines[
                            nk * (optados["num_bands"] + 1) + 1 + nb
                        ].split()[1:]
                    ]
                )
            optados["eigenvalues"].append(eigs)
            optados["projector_weights"].append(pdis)

        optados["projector_weights"] = np.asarray(optados["projector_weights"])

    else:
        optados["dos"] = data[:, 1]

    return optados, True


@scraper_function
def phonon2dict(fname, **kwargs):
    """Parse a CASTEP phonon file into a dictionary.

    Parameters:
        fname (str/list): phonon filename or list of filenames.

    Returns:
        (tuple): containing either dict/str containing data or error, and a bool stating
            if the scrape was successful.

    """

    flines, fname = get_flines_extension_agnostic(fname, ["phonon", "phonon_dos"])

    ph = dict()
    ph["source"] = [fname]

    dos_present = False
    data_start = 0

    for line_no, line in enumerate(flines):
        line = line.lower()
        if "number of ions" in line:
            ph["num_atoms"] = int(line.split()[-1])
        elif "number of branches" in line:
            ph["num_modes"] = int(line.split()[-1])
            ph["num_branches"] = ph["num_modes"]
        elif "number of wavevectors" in line:
            ph["num_kpoints"] = int(line.split()[-1])
        elif "frequencies in" in line:
            ph["freq_unit"] = line.split()[-1]
        elif "unit cell vectors" in line:
            ph["lattice_cart"] = []
            for i in range(3):
                ph["lattice_cart"].append(
                    [f90_float_parse(elem) for elem in flines[line_no + i + 1].split()]
                )
            ph["lattice_abc"] = cart2abc(ph["lattice_cart"])
        elif "fractional co-ordinates" in line:
            ph["positions_frac"] = []
            ph["atom_types"] = []
            ph["atom_masses"] = []
            i = 1
            while "END header" not in flines[line_no + i]:
                ph["positions_frac"].append(
                    [f90_float_parse(elem) for elem in flines[line_no + i].split()[1:4]]
                )
                ph["atom_types"].append(flines[line_no + i].split()[-2])
                ph["atom_masses"].append(
                    f90_float_parse(flines[line_no + i].split()[-1])
                )
                i += 1
        elif "end header" in line:
            data_start = line_no + 1
        elif "begin dos" in line:
            dos_present = True
            projector_labels = flines[line_no].split()[5:]
            projector_labels = [(label, None, None) for label in projector_labels]
            begin_dos = line_no + 1

        elif "q-pt" in line:
            last_qpt_ind = int(line.split()[1])

    if dos_present:
        # extra header line with GRADIENTS written when dos is present
        data_start += 1
        # no eigenvectors written when dos is present
        line_offset = ph["num_modes"] + 1
    else:
        line_offset = ph["num_modes"] * (ph["num_atoms"] + 1) + 3

    data = flines[data_start:]

    ph["phonon_kpoint_list"] = []
    if "num_kpoints" not in ph:
        ph["num_kpoints"] = last_qpt_ind
    ph["eigenvalues_q"] = np.zeros((1, ph["num_modes"], ph["num_kpoints"]))
    raman_intensity = np.zeros_like(ph["eigenvalues_q"])
    infrared_intensity = np.zeros_like(ph["eigenvalues_q"])
    raman = False
    ir = False
    for qind in range(ph["num_kpoints"]):
        ph["phonon_kpoint_list"].append(
            [f90_float_parse(elem) for elem in data[qind * line_offset].split()[2:]]
        )
        for i in range(1, ph["num_modes"] + 1):
            line_split = data[qind * line_offset + i].split()
            ph["eigenvalues_q"][0][i - 1][qind] = f90_float_parse(line_split[1])
            if len(line_split) > 2:
                infrared_intensity[0][i - 1][qind] = f90_float_parse(line_split[2])
                ir = True
            if len(line_split) > 3:
                raman_intensity[0][i - 1][qind] = f90_float_parse(line_split[3])
                raman = True

    if ir:
        ph["infrared_intensity"] = infrared_intensity
    if raman:
        ph["raman_intensity"] = raman_intensity

    if dos_present:
        # remove header and "END"
        flines = flines[begin_dos:-1]
        raw_data = np.genfromtxt(flines)
        ph["energies"] = raw_data[:, 0] * INVERSE_CM_TO_EV
        ph["dos"] = raw_data[:, 1]
        ph["pdos"] = dict()
        ph["pdos"]["pdos"] = dict()
        ph["pdos"]["projectors"] = []
        ph["pdos"]["energies"] = ph["energies"]

        for i, label in enumerate(projector_labels):
            ph["pdos"]["pdos"][label] = raw_data[:, i + 2]
            ph["pdos"]["projectors"].append(label)

    ph["kpoint_path"] = np.asarray([qpt[0:3] for qpt in ph["phonon_kpoint_list"]])
    ph["kpoint_weights"] = [qpt[3] for qpt in ph["phonon_kpoint_list"]]
    ph["eigenvalues_q"] *= INVERSE_CM_TO_EV
    ph["softest_mode_freq"] = np.min(ph["eigenvalues_q"])

    if ph["softest_mode_freq"] < -1:
        import warnings

        warnings.warn(
            f"File {fname} has negative eigenvalue {ph['softest_mode_freq']}."
        )
    ph["eigs_q"] = ph["eigenvalues_q"]

    return ph, True


@scraper_function
def phonon_dos2dict(*args, **kwargs):
    """Wrapper for old phonon DOS scraper, which has since been merged
    with `phonon2dict`. Note that this function still has a different
    effect to `phonon2dict` when `as_model` is used as the results will
    be cast into a :class:`VibrationalDOS` object.

    """
    return phonon2dict(*args, no_wrap=True, **kwargs)


def usp2dict(fname, **kwargs):
    """Extract pseudopotential string from a CASTEP
    OTF .USP file.

    Parameters:
        fname (str/list): filename of usp file, or list of filenames.

    Returns:
        dict: partial species_pot dict from usp file.

    """
    species_pot = dict()
    with open(fname, "r") as f:
        flines = f.readlines()
        for line_no, line in enumerate(flines):
            if "Pseudopotential Report" in line:
                i = 0
                while i + line_no < len(flines) - 3:
                    if "Pseudopotential Report" in flines[line_no + i]:
                        i += 2
                        elem = flines[line_no + i].split(":")[1].split()[0]
                    elif "core correction" in flines[line_no + i]:
                        i += 2
                        species_pot[elem] = flines[line_no + i].strip().split()[1]
                        # check next line for wrapped definition
                        if flines[line_no + i + 1].strip().startswith("--------"):
                            break
                        else:
                            species_pot[elem] += (
                                flines[line_no + i + 1].strip().split()[1]
                            )
                            break
                    i += 1
    species_pot[elem] = species_pot[elem].replace('"', "")
    species_pot[elem] = species_pot[elem].replace("[]", "")
    return species_pot


def _castep_scrape_thermo_data(flines, castep):
    """Scrape the data from a CASTEP Thermodynamics claculation.

    Note:
        This only scrapes from Thermodynamics section currently,
        NOT Atomic Displacement Parameters

    Parameters:
        flines (list): list of lines contained in file
        castep (dict): dictionary to update with data

    Returns:
        dict: dictionary updated with scraped thermodynamics data
              in format doc['energy_type'] = {300:-0.01,200:-0.05}

    """
    for line_no, line in enumerate(flines):
        if "Number of temperature values" in line:
            castep["thermo_num_temp_vals"] = int(line.split(":")[-1].strip())
        elif "Initial temperature" in line:
            castep["thermo_temp_init"] = f90_float_parse(
                line.split(":")[1].strip().split(" ")[0]
            )
        elif "Final temperature" in line:
            castep["thermo_temp_final"] = f90_float_parse(
                line.split(":")[1].strip().split(" ")[0]
            )
        elif "Spacing between temperature values" in line:
            castep["thermo_temp_spacing"] = f90_float_parse(
                line.split(":")[1].strip().split(" ")[0]
            )
        elif "Zero-point energy" in line:
            castep["thermo_zero_point_energy"] = f90_float_parse(
                line.split("=")[1].strip().split(" ")[0]
            )
        elif "T(K)" and "E(eV)" in line:
            castep["thermo_temps"] = []  # temperatures calculation was done at
            castep["thermo_enthalpy"] = {}  # enthalpy E(eV)
            castep["thermo_free_energy"] = {}  # free energy F(eV)
            castep["thermo_entropy"] = {}  # entropy S(J/mol/K)
            castep["thermo_heat_cap"] = {}  # heat capacity Cv(J/mol/K)
            i = 2
            while True:
                if not flines[line_no + i + 1].strip():
                    break
                else:
                    temp_line = flines[line_no + i].split()
                    castep["thermo_temps"].append(f90_float_parse(temp_line[0]))
                    castep["thermo_enthalpy"][
                        f90_float_parse(temp_line[0])
                    ] = f90_float_parse(temp_line[1])
                    castep["thermo_free_energy"][
                        f90_float_parse(temp_line[0])
                    ] = f90_float_parse(temp_line[2])
                    castep["thermo_entropy"][
                        f90_float_parse(temp_line[0])
                    ] = f90_float_parse(temp_line[3])
                    castep["thermo_heat_cap"][
                        f90_float_parse(temp_line[0])
                    ] = f90_float_parse(temp_line[4])
                i += 1


def _castep_scrape_phonon_frequencies(flines, castep):
    """Iterate through flines to scrape the phonon frequencies
    from this CASTEP calculation. Will only scrape the *final* set
    of frequencies in the CASTEP file, ignoring any others.

    """
    phonons = {}

    phonons["phonon_fine_kpoint_list"] = []
    phonons["phonon_fine_kpoint_weights"] = []
    phonons["eigs_q"] = []

    # first find the last block of frequencies
    for line_no, line in enumerate(flines):
        if "Performing frequency calculation at " in line:
            start_line_no = line_no + 2

    q_pt_ind = 0
    for line_no, line in enumerate(flines[start_line_no:]):
        if "============" in line:
            break
        if "q-pt=" in line:
            q_pt = [
                f90_float_parse(val)
                for val in line.split("(")[-1].split(")")[0].split()
            ]
            phonons["phonon_fine_kpoint_list"].append(q_pt)
            phonons["phonon_fine_kpoint_weights"].append(
                f90_float_parse(line.split()[-2])
            )
            phonons["eigs_q"].append([])

            for _, freq_line in enumerate(flines[start_line_no:][line_no + 6 :]):
                if "........................." in freq_line:
                    break
                phonons["eigs_q"][q_pt_ind].append(
                    f90_float_parse(freq_line.split()[2])
                )
            q_pt_ind += 1

    phonons["num_modes"] = len(phonons["eigs_q"][0])
    phonons["eigs_q"] = np.asarray(phonons["eigs_q"]).T
    phonons["eigs_q"] = INVERSE_CM_TO_EV * phonons["eigs_q"].reshape(
        1, *np.shape(phonons["eigs_q"])
    )
    phonons["kpoint_path"] = phonons["phonon_fine_kpoint_list"]
    phonons["kpoint_weights"] = phonons["phonon_fine_kpoint_weights"]
    phonons["num_kpoints"] = len(phonons["phonon_fine_kpoint_list"])
    phonons["num_qpoints"] = len(phonons["phonon_fine_kpoint_list"])
    phonons["softest_mode_freq"] = np.min(phonons["eigs_q"])
    if phonons["softest_mode_freq"] < -1:
        import warnings

        warnings.warn(
            f"File {castep['source'][0]} has negative eigenvalue {phonons['softest_mode_freq']}."
        )

    for key in phonons:
        castep[key] = phonons[key]


def _castep_scrape_atoms(flines, castep):
    """Iterate forwards through flines to scrape atomic types and
    initial positions, and get a preliminary lattice.

    Parameters:
        flines (list): list of lines in file
        castep (dict): dictionary to update with scraped data

    Returns:
        dict: dictionary updated with scraped data

    """
    for line_no, line in enumerate(flines):
        if "Real Lattice" in line:
            castep["lattice_cart"] = []
            i = 1
            while True:
                if not flines[line_no + i].strip():
                    break
                else:
                    temp_line = flines[line_no + i].split()[0:3]
                    castep["lattice_cart"].append(list(map(f90_float_parse, temp_line)))
                i += 1
        if "Lattice parameters" in line:
            castep["lattice_abc"] = []
            i = 1
            castep["lattice_abc"].append(
                list(
                    map(
                        f90_float_parse,
                        [
                            flines[line_no + i].split("=")[1].strip().split(" ")[0],
                            flines[line_no + i + 1].split("=")[1].strip().split(" ")[0],
                            flines[line_no + i + 2].split("=")[1].strip().split(" ")[0],
                        ],
                    )
                )
            )
            castep["lattice_abc"].append(
                list(
                    map(
                        f90_float_parse,
                        [
                            flines[line_no + i].split("=")[-1].strip(),
                            flines[line_no + i + 1].split("=")[-1].strip(),
                            flines[line_no + i + 2].split("=")[-1].strip(),
                        ],
                    )
                )
            )
        if "Current cell volume" in line:
            castep["cell_volume"] = f90_float_parse(
                line.split("=")[1].split()[0].strip()
            )
        if "atom types" not in castep and "Cell Contents" in line:
            castep["atom_types"] = []
            castep["positions_frac"] = []
            i = 1
            atoms = False
            while True:
                if atoms:
                    if "xxxxxxxxx" in flines[line_no + i]:
                        atoms = False
                        break
                    else:
                        castep["atom_types"].append(flines[line_no + i].split()[1])
                        castep["positions_frac"].append(
                            list(
                                map(f90_float_parse, (flines[line_no + i].split()[3:6]))
                            )
                        )
                if "x------" in flines[line_no + i]:
                    atoms = True
                i += 1
            for ind, pos in enumerate(castep["positions_frac"]):
                for k in range(3):
                    if pos[k] > 1 or pos[k] < 0:
                        castep["positions_frac"][ind][k] %= 1
            castep["num_atoms"] = len(castep["atom_types"])
            castep["stoichiometry"] = get_stoich(castep["atom_types"])
            castep["num_fu"] = castep["num_atoms"] / sum(
                [elem[1] for elem in castep["stoichiometry"]]
            )
            break
    else:
        raise ComputationError("Unable to find atoms in CASTEP file.")


def _castep_scrape_final_parameters(flines, castep):
    """Scrape the DFT parameters and metadata from a CASTEP file,
    using those listed last in the file
    (i.e. those used to make the final structure).

    Parameters:
        flines (list): list of lines in the file
        castep (dict): dictionary in which to put scraped data

    Returns:
        dict: dictionary updated with scraped data

    """

    pspot_report_dict = dict()
    for line_no, line in enumerate(reversed(flines)):
        line_no = len(flines) - 1 - line_no
        if "castep_version" not in castep and "Release CASTEP version" in line:
            # version is first thing printed in file, so break once this is read
            castep["castep_version"] = line.replace("|", "").split()[-1]
            break
        elif "date" not in castep and "Run started:" in line:
            try:
                year = line.split()[5]
                month = str(strptime(line.split()[4], "%b").tm_mon)
                day = line.split()[3]
                castep["date"] = day + "-" + month + "-" + year
            except (IndexError, ValueError):
                castep["date"] = "unknown"
        elif "_compiler_architecture" not in castep and "compiled for" in line.lower():
            try:
                castep["_compiler_architecture"] = line.split()[2]
            except IndexError:
                castep["_compiler_architecture"] = "unknown"
        elif "_castep_commit" not in castep and " from code version" in line.lower():
            try:
                castep["_castep_commit"] = " ".join(
                    flines[line_no : line_no + 2]
                ).split()[3]
            except (IndexError, ValueError):
                castep["_castep_commit"] = "unknown"
        elif "task" not in castep and "type of calculation" in line:
            castep["task"] = line.split(":")[-1].strip().replace(" ", "")
        elif "xc_functional" not in castep and "functional" in line:
            # convert from .castep file xc_functional to param style
            xc_string = line.split(":")[-1].strip()
            if "Local Density Approximation" in xc_string:
                castep["xc_functional"] = "LDA"
            elif "Perdew Burke Ernzerhof" in xc_string:
                castep["xc_functional"] = "PBE"
            elif "PBE for solids" in xc_string:
                castep["xc_functional"] = "PBESol"
            elif "hybrid B3LYP" in xc_string:
                castep["xc_functional"] = "B3LYP"
            elif "hybrid HSE03" in xc_string:
                castep["xc_functional"] = "HSE03"
            elif "hybrid HSE06" in xc_string:
                castep["xc_functional"] = "HSE06"
            elif "RSCAN" in xc_string:
                castep["xc_functional"] = "RSCAN"
            elif "Screened Hartree-Fock" in xc_string:
                castep["xc_functional"] = "SHF-LDA"
            else:
                castep["xc_functional"] = xc_string.split(":")[-1]
                warnings.warn(
                    "Unrecognised functional {xc}: scraping as {xc}."
                    "This may lead to incompatible param files.".format(
                        xc=castep["xc_functional"]
                    )
                )
        elif "cut_off_energy" not in castep and "plane wave basis set" in line:
            castep["cut_off_energy"] = f90_float_parse(line.split(":")[-1].split()[0])
        elif (
            "finite_basis_corr" not in castep
            and "finite basis set correction  " in line
        ):
            castep["finite_basis_corr"] = line.split(":")[-1].strip()
        elif "kpoints_mp_grid" not in castep and "MP grid size for SCF" in line:
            castep["kpoints_mp_grid"] = list(
                map(int, list(line.split("is")[-1].split()))
            )
        elif "num_kpoints" not in castep and "Number of kpoints used" in line:
            castep["num_kpoints"] = int(line.split()[-1])
        elif "geom_force_tol" not in castep and "max ionic |force| tolerance" in line:
            castep["geom_force_tol"] = f90_float_parse(line.split()[-2])
        elif (
            "elec_energy_tol" not in castep
            and "total energy / atom convergence tol" in line
        ):
            castep["elec_energy_tol"] = f90_float_parse(line.split()[-2])
        elif (
            "sedc_apply" not in castep
            and "DFT+D: Semi-empirical dispersion correction    : on" in line
        ):
            castep["sedc_apply"] = True
            castep["sedc_scheme"] = flines[line_no + 1].split(":")[1].split()[0]
        elif "space_group" not in castep and "Space group of crystal" in line:
            space_group = line.split(":")[-1].split(",")[0].strip().replace(" ", "")
            if space_group:
                castep["space_group"] = space_group
        elif "nelectrons" not in castep and "number of  electrons" in line:
            castep["nelectrons"] = f90_float_parse(line.split(":")[-1])
        elif "nelectrons" not in castep and "number of bands" in line:
            castep["nbands"] = int(line.split(":")[-1])
        elif "Cell constraints are" in line and "cell_constraints" not in castep:
            castep["cell_constraints"] = [
                int(val) for val in line.split(":")[-1].split()
            ]
            if castep["cell_constraints"] == [1, 2, 3, 4, 5, 6]:
                del castep["cell_constraints"]
            elif all(val == 0 for val in castep["cell_constraints"]):
                castep["fix_all_cell"] = True
                del castep["cell_constraints"]
            if "cell_constraints" in castep:
                if len(castep["cell_constraints"]) != 6:
                    raise RuntimeError("Unable to read cell constraints block.")
                castep["cell_constraints"] = [
                    castep["cell_constraints"][0:3],
                    castep["cell_constraints"][3:],
                ]
        elif "external_pressure" not in castep and "External pressure/stress" in line:
            try:
                castep["external_pressure"] = []
                for i in range(3):
                    parsed_line = list(
                        map(f90_float_parse, flines[line_no + i + 1].split())
                    )
                    if len(parsed_line) != 3:
                        parsed_line = i * [0.0] + parsed_line
                    castep["external_pressure"].append(parsed_line)
            except ValueError:
                castep["external_pressure"] = [
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                ]
        elif (
            "spin_polarized" not in castep
            and "treating system as spin-polarized" in line
        ):
            castep["spin_polarized"] = True
        elif "hubbard_u" not in castep and "Hubbard U values are eV" in line:
            castep["hubbard_u"] = defaultdict(list)
            n_lines_header = 5
            i = 0
            while i < castep["num_atoms"]:
                line = flines[line_no + i + n_lines_header].strip()
                atom = line.split()[0].replace("|", "")
                shifts = list(map(f90_float_parse, line.split()[-5:-1]))
                for ind, shift in enumerate(shifts):
                    if shift != 0:
                        if atom not in castep["hubbard_u"]:
                            castep["hubbard_u"][atom] = dict()
                        if ind == 0:
                            orbital = "s"
                        elif ind == 1:
                            orbital = "p"
                        elif ind == 2:
                            orbital = "d"
                        elif ind == 3:
                            orbital = "f"
                        castep["hubbard_u"][atom][orbital] = shift
                i += 1
        elif "Pseudopotential Report" in line:
            if "species_pot" not in castep:
                castep["species_pot"] = dict()
            i = 0
            while i + line_no < len(flines) - 3:
                if "Pseudopotential Report" in flines[line_no + i]:
                    i += 2
                    elem = flines[line_no + i].split(":")[1].split()[0]

                elif "core correction" in flines[line_no + i]:
                    i += 2
                    if not pspot_report_dict.get(elem):
                        castep["species_pot"][elem] = (
                            flines[line_no + i].split('"')[1].replace("[]", "")
                        )
                        pspot_report_dict[elem] = True
                i += 1
        elif "species_pot" not in castep and "Files used for pseudopotentials" in line:
            if "species_pot" not in castep:
                castep["species_pot"] = dict()
            i = 1
            while True:
                if not flines[line_no + i].strip():
                    break
                else:
                    elem = flines[line_no + i].split()[0].strip()
                    if not pspot_report_dict.get(elem):
                        castep["species_pot"][elem] = (
                            flines[line_no + i].split()[1].split("/")[-1]
                        )
                        if castep["species_pot"][elem] == "Pseudopotential":
                            castep["species_pot"][elem] = (
                                flines[line_no + i].split()[0].strip()
                            )
                            castep["species_pot"][elem] += "_OTF.usp"
                        pspot_report_dict[elem] = False
                i += 1
        elif "peak_mem_MB" not in castep and "Peak Memory Use" in line:
            try:
                castep["peak_mem_MB"] = int(f90_float_parse(line.split()[-2]) / 1024)
            except (IndexError, ValueError):
                castep["peak_mem_MB"] = -1
        elif (
            "estimated_mem_per_process_MB" not in castep
            and "total storage required per process" in line
        ):
            try:
                castep["estimated_mem_per_process_MB"] = f90_float_parse(
                    line.split()[-5]
                )
            except (IndexError, ValueError):
                castep["estimated_mem_per_process_MB"] = 0
        elif (
            "num_mpi_processes" not in castep
            and "Calculation parallelised over" in line
        ):
            try:
                castep["num_mpi_processes"] = int(f90_float_parse(line.split()[3]))
            except Exception:
                castep["num_mpi_processes"] = 1
        elif "num_mpi_processes" not in castep and "Calculation not parall" in line:
            castep["num_mpi_processes"] = 1
        elif "_xc_beef" not in castep and "BEEF completed" in line:
            castep["_xc_beef"] = True

    # write zero pressure if not found in file
    if "external_pressure" not in castep:
        castep["external_pressure"] = [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    if "spin_polarized" not in castep:
        castep["spin_polarized"] = False
    if "num_mpi_processes" in castep and "estimated_mem_per_process_MB" in castep:
        castep["estimated_mem_MB"] = (
            castep["estimated_mem_per_process_MB"] * castep["num_mpi_processes"]
        )


def _castep_scrape_final_structure(flines, castep, db=True):
    """Scrape final structure from CASTEP file.

    Parameters:
        flines (list): list of lines contained in file
        castep (dict): dictionary to update with data

    Keyword arguments:
        db (bool): whether to enforce database style, e.g. geometry optimisation only

    Returns:
        dict: dictionary updated with scraped data

    """
    if "task" in castep and castep["task"].strip() == "geometryoptimization":
        castep["optimised"] = False
        finish_line, castep["optimised"] = _castep_find_final_structure(flines)
    elif not db:
        finish_line = 0

    final_flines = flines[finish_line + 1 :]
    for line_no, line in enumerate(final_flines):
        if "Real Lattice" in line:
            castep["lattice_cart"] = []
            i = 1
            while i < 4:
                if not final_flines[line_no + i].strip():
                    break
                else:
                    temp_line = final_flines[line_no + i].split()[0:3]
                    castep["lattice_cart"].append(list(map(f90_float_parse, temp_line)))
                i += 1
        if "Lattice parameters" in line:
            castep["lattice_abc"] = []
            i = 1
            castep["lattice_abc"].append(
                list(
                    map(
                        f90_float_parse,
                        [
                            final_flines[line_no + i]
                            .split("=")[1]
                            .strip()
                            .split(" ")[0],
                            final_flines[line_no + i + 1]
                            .split("=")[1]
                            .strip()
                            .split(" ")[0],
                            final_flines[line_no + i + 2]
                            .split("=")[1]
                            .strip()
                            .split(" ")[0],
                        ],
                    )
                )
            )
            castep["lattice_abc"].append(
                list(
                    map(
                        f90_float_parse,
                        [
                            final_flines[line_no + i].split("=")[-1].strip(),
                            final_flines[line_no + i + 1].split("=")[-1].strip(),
                            final_flines[line_no + i + 2].split("=")[-1].strip(),
                        ],
                    )
                )
            )
        elif "Current cell volume" in line:
            castep["cell_volume"] = f90_float_parse(
                line.split("=")[1].split()[0].strip()
            )
        elif "Cell Contents" in line:
            castep["positions_frac"] = []
            i = 1
            atoms = False
            while True:
                if atoms:
                    if "xxxxxxxxx" in final_flines[line_no + i]:
                        atoms = False
                        break
                    else:
                        temp_frac = final_flines[line_no + i].split()[3:6]
                        castep["positions_frac"].append(
                            list(map(f90_float_parse, temp_frac))
                        )
                if "x------" in final_flines[line_no + i]:
                    atoms = True
                i += 1
        # don't check if final_energy exists, as this will update for each GO step
        elif "Final energy =" in line or "Final energy, E" in line:
            castep["total_energy"] = f90_float_parse(line.split("=")[1].split()[0])
            castep["total_energy_per_atom"] = (
                castep["total_energy"] / castep["num_atoms"]
            )
        elif "Final free energy" in line:
            castep["smeared_free_energy"] = f90_float_parse(
                line.split("=")[1].split()[0]
            )
            castep["smeared_free_energy_per_atom"] = (
                castep["smeared_free_energy"] / castep["num_atoms"]
            )
        elif "0K energy" in line:
            castep["0K_energy"] = f90_float_parse(line.split("=")[1].split()[0])
            castep["0K_energy_per_atom"] = castep["0K_energy"] / castep["num_atoms"]
        elif "(SEDC) Total Energy Correction" in line:
            castep["dispersion_correction_energy"] = f90_float_parse(
                line.split(":")[1].split()[0]
            )
        elif "Dispersion corrected final energy" in line:
            castep["dispersion_corrected_energy"] = f90_float_parse(
                line.split("=")[1].split()[0]
            )
            castep["dispersion_corrected_energy_per_atom"] = (
                castep["dispersion_corrected_energy"] / castep["num_atoms"]
            )
        elif "Dispersion corrected final free energy" in line:
            castep["dispersion_corrected_free_energy"] = f90_float_parse(
                line.split("=")[1].split()[0]
            )
            castep["dispersion_corrected_free_energy_per_atom"] = (
                castep["dispersion_corrected_free_energy"] / castep["num_atoms"]
            )
        elif "Dispersion corrected est. 0K energy" in line:
            castep["dispersion_corrected_0K_energy"] = f90_float_parse(
                line.split("=")[1].split()[0]
            )
            castep["dispersion_corrected_0K_energy_per_atom"] = (
                castep["dispersion_corrected_0K_energy"] / castep["num_atoms"]
            )
        elif " Forces **" in line:
            castep["forces"] = []
            i = 1
            forces = False
            while True:
                if forces:
                    if "*" in final_flines[line_no + i].split()[1]:
                        forces = False
                        break
                    else:
                        castep["forces"].append([])
                        for j in range(3):
                            temp = final_flines[line_no + i].replace("(cons'd)", "")
                            castep["forces"][-1].append(
                                f90_float_parse(temp.split()[3 + j])
                            )
                elif "x" in final_flines[line_no + i]:
                    i += 1  # skip next blank line
                    forces = True
                i += 1
            castep["max_force_on_atom"] = np.max(
                np.linalg.norm(castep["forces"], axis=-1)
            )
        elif "Stress Tensor" in line:
            if "Constrained" in line:
                stress_key = "constrained_stress"
            else:
                stress_key = "stress"
            i = 1
            while i < 20:
                if "Cartesian components" in final_flines[line_no + i]:
                    castep[stress_key] = []
                    for j in range(3):
                        castep[stress_key].append(
                            list(
                                map(
                                    f90_float_parse,
                                    (final_flines[line_no + i + j + 4].split()[2:5]),
                                )
                            )
                        )
                elif "Pressure" in final_flines[line_no + i]:
                    try:
                        castep["pressure"] = f90_float_parse(
                            final_flines[line_no + i].split()[-2]
                        )
                    except ValueError:
                        pass
                    break
                i += 1
        elif "Integrated Spin Density" in line:
            castep["integrated_spin_density"] = f90_float_parse(line.split()[-2])
        elif "Integrated |Spin Density|" in line:
            castep["integrated_mod_spin_density"] = f90_float_parse(line.split()[-2])
        elif "Atomic Populations (Mulliken)" in line:
            # population format seems to change every CASTEP version...
            if float(castep.get("castep_version", 0.0)) >= 18:
                if castep["spin_polarized"]:
                    castep["mulliken_spins"] = []
                    castep["mulliken_net_spin"] = 0.0
                    castep["mulliken_abs_spin"] = 0.0
                castep["mulliken_charges"] = []
                castep["mulliken_spins"] = []
                i = 0
                ind = 0
                while ind < len(castep["atom_types"]):
                    if castep["spin_polarized"]:
                        castep["mulliken_charges"].append(
                            f90_float_parse(final_flines[line_no + i + 4].split()[-2])
                        )
                        castep["mulliken_spins"].append(
                            f90_float_parse(final_flines[line_no + i + 4].split()[-1])
                        )
                        castep["mulliken_net_spin"] += castep["mulliken_spins"][-1]
                        castep["mulliken_abs_spin"] += abs(castep["mulliken_spins"][-1])
                        i += 2
                    else:
                        castep["mulliken_charges"].append(
                            f90_float_parse(final_flines[line_no + i + 4].split()[-1])
                        )
                        i += 1
                    ind += 1
        elif "Final Enthalpy" in line:
            castep["enthalpy"] = f90_float_parse(line.split("=")[-1].split()[0])
            castep["enthalpy_per_atom"] = (
                f90_float_parse(line.split("=")[-1].split()[0]) / castep["num_atoms"]
            )
        elif "Final bulk modulus" in line:
            try:
                castep["bulk_modulus"] = f90_float_parse(line.split("=")[-1].split()[0])
            except ValueError:
                # the above will fail if bulk modulus was not printed (i.e. if it was unchanged)
                pass

        elif (
            "Chemical Shielding and Electric Field Gradient Tensors".lower()
            in line.lower()
        ):
            i = 5
            castep["chemical_shifts"] = []
            while True:
                # break when the line containing just '=' is reached
                if len(flines[line_no + i].split()) == 1:
                    break
                castep["chemical_shifts"].append(flines[line_no + i].split()[3])
                i += 1
            if len(castep["chemical_shifts"]) != len(castep["atom_types"]):
                raise RuntimeError(
                    "Found fewer chemical shifts than atoms (or vice versa)!"
                )

    # calculate kpoint spacing if not found
    if (
        "kpoints_mp_grid" in castep
        and "kpoints_mp_spacing" not in castep
        and "lattice_cart" in castep
    ):
        castep["kpoints_mp_spacing"] = calc_mp_spacing(
            castep["lattice_cart"], castep["kpoints_mp_grid"], prec=4
        )
    return castep


def _castep_scrape_metadata(flines, castep):
    """Scrape addtitional timing metadata from CASTEP file.

    Parameters:
        flines (list): list of lines contained in file
        castep (dict): dictionary to update with data

    Returns:
        dict: dictionary updated with scraped data

    """
    # computing metadata, i.e. parallelism, time, memory, version
    if "total_time_secs" not in castep:
        castep["total_time_secs"] = 0
    if "total_time_hrs" not in castep:
        castep["total_time_hrs"] = 0

    for ind, line in enumerate(flines):
        if "Total time" in line and "matrix elements" not in line:
            try:
                time = f90_float_parse(line.split()[-2])
                castep["total_time_secs"] += time
                castep["total_time_hrs"] += time / 3600
            except (IndexError, ValueError):
                castep["final_calculation_time_secs"] = 0
        elif "Calculation only took" in line:
            castep["_time_estimated"] = f90_float_parse(line.split()[4])


def _castep_find_final_structure(flines):
    """Search for info on final structure in .castep file.

    Parameters:
        flines (list): list of lines in file.

    Returns:
        int: line number in file where total energy of final structure is printed.

    """
    optimised = False
    finish_line = 0
    success_string = "Geometry optimization completed successfully"
    failure_string = "Geometry optimization failed to converge after"
    annoying_string = "WARNING - there is nothing to optimise - skipping relaxation"
    # look for final "success/failure" string in file for geometry optimisation
    for line_no, line in enumerate(reversed(flines)):
        if success_string in line:
            finish_line = len(flines) - line_no
            optimised = True
            break
        if annoying_string in line:
            finish_line = len(flines) - line_no
            optimised = True
            break
        if failure_string in line:
            finish_line = len(flines) - line_no
            optimised = False
            break

    # now wind back to get final total energies and non-symmetrised forces
    for count, line in enumerate(reversed(flines[:finish_line])):
        if "Final energy, E" in line or "Final energy =" in line:
            finish_line -= count + 2
            break

    return finish_line, optimised


def _castep_finalize_snapshot(snapshot: dict, intermediates: list) -> None:
    """Add per-atom keys to snapshot and append it to the intermediates list.

    Parameters:
        snapshot: The document containing the current snapshot.
        intermediates: The list of snapshots so far.

    """
    # if positions frac or lattice didn't change (and thus weren't printed, use the last value)
    if "positions_frac" not in snapshot:
        snapshot["positions_frac"] = intermediates[-1]["positions_frac"]
        snapshot["atom_types"] = intermediates[-1]["atom_types"]
        snapshot["num_atoms"] = len(snapshot["positions_frac"])
    if "lattice_cart" not in snapshot:
        snapshot["lattice_cart"] = intermediates[-1]["lattice_cart"]
        snapshot["lattice_abc"] = intermediates[-1]["lattice_abc"]

    snapshot["smeared_free_energy_per_atom"] = (
        snapshot["smeared_free_energy"] / snapshot["num_atoms"]
    )
    snapshot["0K_energy_per_atom"] = snapshot["0K_energy"] / snapshot["num_atoms"]
    snapshot["total_energy_per_atom"] = snapshot["total_energy"] / snapshot["num_atoms"]
    # handle single atom forces edge-case
    if snapshot["num_atoms"] == 1:
        snapshot["forces"] = [[0, 0, 0]]

    intermediates.append(snapshot)


def _castep_scrape_all_snapshots(flines, intermediates=False):
    """Scrape all intermediate structures from a CASTEP file, both
    geometry optimisation snapshots, and repeated SCF calculations.

    Parameters:
        flines (list): list of lines of file.
        intermediates (bool): whether to save the snapshots or just
            count them.

    Returns:
        :obj:`list` of :obj:`dict`: list of dictionaries containing
            intermediate snapshots (will be empty if `intermediates`
            is false).
        int: number of completed geometry optimisation steps.

    """
    intermediates_list = []
    num_opt_steps = 0
    snapshot = dict()
    if intermediates:
        for line_no, line in enumerate(flines):
            try:
                if "Cell Contents" in line or "Unit Cell" in line:
                    if "total_energy" in snapshot:
                        _castep_finalize_snapshot(snapshot, intermediates_list)
                        snapshot = dict()

                if "Cell Contents" in line:
                    snapshot["positions_frac"] = []
                    snapshot["atom_types"] = []
                    i = 1
                    atoms = False
                    while True:
                        if atoms:
                            if "xxxxxxxxx" in flines[line_no + i]:
                                atoms = False
                                break
                            else:
                                temp_frac = flines[line_no + i].split()[3:6]
                                snapshot["positions_frac"].append(
                                    list(map(f90_float_parse, temp_frac))
                                )
                                snapshot["atom_types"].append(
                                    flines[line_no + i].split()[1]
                                )
                        if "x------" in flines[line_no + i]:
                            atoms = True
                        i += 1
                    for ind, pos in enumerate(snapshot["positions_frac"]):
                        for k in range(3):
                            if pos[k] > 1 or pos[k] < 0:
                                snapshot["positions_frac"][ind][k] %= 1

                    snapshot["num_atoms"] = len(snapshot["positions_frac"])
                    snapshot["stoichiometry"] = get_stoich(snapshot["atom_types"])

                elif "Real Lattice" in line:
                    snapshot["lattice_cart"] = []
                    i = 1
                    while True:
                        if not flines[line_no + i].strip():
                            break
                        else:
                            temp_line = flines[line_no + i].split()[0:3]
                            snapshot["lattice_cart"].append(
                                list(map(f90_float_parse, temp_line))
                            )
                        i += 1

                elif "Lattice parameters" in line:
                    snapshot["lattice_abc"] = []
                    i = 1
                    snapshot["lattice_abc"].append(
                        list(
                            map(
                                f90_float_parse,
                                [
                                    flines[line_no + i]
                                    .split("=")[1]
                                    .strip()
                                    .split(" ")[0],
                                    flines[line_no + i + 1]
                                    .split("=")[1]
                                    .strip()
                                    .split(" ")[0],
                                    flines[line_no + i + 2]
                                    .split("=")[1]
                                    .strip()
                                    .split(" ")[0],
                                ],
                            )
                        )
                    )
                    snapshot["lattice_abc"].append(
                        list(
                            map(
                                f90_float_parse,
                                [
                                    flines[line_no + i].split("=")[-1].strip(),
                                    flines[line_no + i + 1].split("=")[-1].strip(),
                                    flines[line_no + i + 2].split("=")[-1].strip(),
                                ],
                            )
                        )
                    )

                elif "Current cell volume" in line:
                    snapshot["cell_volume"] = f90_float_parse(
                        line.split("=")[1].split()[0].strip()
                    )
                elif "Final energy, E" in line or "Final energy =" in line:
                    snapshot["total_energy"] = f90_float_parse(
                        line.split("=")[1].split()[0]
                    )
                elif "Final free energy" in line:
                    snapshot["smeared_free_energy"] = f90_float_parse(
                        line.split("=")[1].split()[0]
                    )
                elif "0K energy" in line:
                    snapshot["0K_energy"] = f90_float_parse(
                        line.split("=")[1].split()[0]
                    )
                elif " Forces **" in line:
                    snapshot["forces"] = []
                    i = 1
                    max_force = 0
                    forces = False
                    while True:
                        if forces:
                            if "*" in flines[line_no + i].split()[1]:
                                forces = False
                                break
                            else:
                                force_on_atom = 0
                                snapshot["forces"].append([])
                                for j in range(3):
                                    temp = flines[line_no + i].replace("(cons'd)", "")
                                    force_on_atom += (
                                        f90_float_parse(temp.split()[3 + j]) ** 2
                                    )
                                    snapshot["forces"][-1].append(
                                        f90_float_parse(temp.split()[3 + j])
                                    )
                                if force_on_atom > max_force:
                                    max_force = force_on_atom
                        elif "x" in flines[line_no + i]:
                            i += 1  # skip next blank line
                            forces = True
                        i += 1
                    snapshot["max_force_on_atom"] = pow(max_force, 0.5)
                elif "Stress Tensor" in line:
                    if "Constrained" in line:
                        stress_key = "constrained_stress"
                    else:
                        stress_key = "stress"
                    i = 1
                    while i < 20:
                        if "Cartesian components" in flines[line_no + i]:
                            snapshot[stress_key] = []
                            for j in range(3):
                                snapshot[stress_key].append(
                                    list(
                                        map(
                                            f90_float_parse,
                                            (flines[line_no + i + j + 4].split()[2:5]),
                                        )
                                    )
                                )
                        elif "Pressure" in flines[line_no + i]:
                            snapshot["pressure"] = f90_float_parse(
                                flines[line_no + i].split()[-2]
                            )
                            break
                        i += 1
                # use only finished iterations for counting number of complete GO steps
                elif ": finished iteration" in line and "with enthalpy" in line:
                    # don't include the "zeroth" step before anything has been moved
                    if "0" not in line.split():
                        num_opt_steps += 1
            except Exception as exc:
                raise RuntimeError from exc
    else:
        for line_no, line in enumerate(flines):
            # use only finished iterations for counting number of complete GO steps
            if ": finished iteration" in line and "with enthalpy" in line:
                # don't include the "zeroth" step before anything has been moved
                if "0" not in line.split():
                    num_opt_steps += 1

    return intermediates_list, num_opt_steps


def _castep_scrape_beef(flines, castep):
    """Scrape the Bayesian error estimate output from CASTEP, storing
    it under the `_beef` key.

    Parameters:
        flines (list): CASTEP output flines to scrape.
        castep (dict): dictionary to update with BEEF output.

    """
    if castep.get("_xc_beef"):
        from matador.utils.chem_utils import HARTREE_TO_EV

        for line_no, line in enumerate(flines):
            if "Bayesian Error Estimate (BEE)" in line.strip():
                beef_start = line_no
                break
        else:
            return
    else:
        return

    castep["_beef"] = {
        "thetas": [],
        "xc_energy": [],
        "total_energy": [],
        "total_energy_per_atom": [],
        "xc_energy_per_atom": [],
    }

    for line_no, line in enumerate(flines[beef_start:]):
        if "<-- BEEF" in line:
            castep["_beef"]["thetas"].append(
                [f90_float_parse(val) for val in line.strip().split()[1:4]]
            )
            castep["_beef"]["xc_energy"].append(
                HARTREE_TO_EV * f90_float_parse(line.strip().split()[4])
            )
            castep["_beef"]["total_energy"].append(
                HARTREE_TO_EV * f90_float_parse(line.strip().split()[6])
            )
            castep["_beef"]["xc_energy_per_atom"].append(
                castep["_beef"]["xc_energy"][-1] / castep["num_atoms"]
            )
            castep["_beef"]["total_energy_per_atom"].append(
                castep["_beef"]["total_energy"][-1] / castep["num_atoms"]
            )
        if "BEEF completed" in line:
            break
    else:
        raise RuntimeError("End of BEEF estimate not found.")

    castep["_beef"]["mean_total_energy"] = np.mean(castep["_beef"]["total_energy"])
    castep["_beef"]["mean_total_energy_per_atom"] = (
        np.mean(castep["_beef"]["total_energy"]) / castep["num_atoms"]
    )
    castep["_beef"]["std_dev_total_energy"] = np.std(castep["_beef"]["total_energy"])
    castep["_beef"]["std_dev_total_energy_per_atom"] = (
        np.std(castep["_beef"]["total_energy"]) / castep["num_atoms"]
    )


def _castep_scrape_devel_code(flines, castep):
    """Scrape the contents of the developer code block and
    extract any information about nanotube encapsulation.
    Searches for then scrapes the last occurence of "Developer
    Code". If nanotube information is found, the `encapsulated`
    flag is set to True.

    Parameters:
        flines (list): CASTEP output flines to scrape.
        castep (dict): dictionary to update with BEEF output.

    """
    last_devel_code = len(flines)
    for line_no, line in enumerate(flines):
        if line.startswith(" *") and "**** Developer Code ****" in line:
            last_devel_code = line_no

    castep["devel_code"] = []
    for line_no, line in enumerate(flines[last_devel_code + 1 :]):
        # look for end of block
        if "*************" in line.strip():
            break

        line = line.strip()
        if line:
            castep["devel_code"].append(line)

    devel_code = "\n".join(castep["devel_code"])
    if "gaussian_cylinder" in devel_code.lower():
        castep["encapsulated"] = True
        for line in castep["devel_code"]:
            if "radius" in line.lower():
                castep["cnt_radius"] = float(line.strip().split()[-1])

    # store as string rather than list
    castep["devel_code"] = devel_code


def _optados_get_projector_labels(header):
    """Get OptaDOS projector labels from a pdos.dat or pdis.dat file,
    returning None-padded labels as `('species', 'ang_mom', 'spin channel')`,
    e.g. ('Li', 's', None) or ('Cr', None, 'up').

    Parameters:
        header (str): the file header containing projector labels.

    Returns:
        list(tuple(str, str, str)): projector labels formatted as above.

    """
    projectors = []
    for ind, line in enumerate(header):
        separators = ["Projector:", "Column:"]
        if any(sep in line for sep in separators):
            is_spin_pdos = "Spin" in header[ind + 1]
            # skip current line and column headings
            j = 2
            elements = []
            ang_mom_channels = []
            spin_channels = []
            # loop over file finding projector header blocks
            while ind + j + 1 < len(header) and not any(
                keyword in header[ind + j + 1] for keyword in separators
            ):
                elements.append(header[ind + j].split()[1])
                ang_mom_channels.append(header[ind + j].split()[3])
                if is_spin_pdos:
                    spin_channels.append(header[ind + j].split()[4].lower())
                j += 1

            projector_label = []

            # check that this projector contains exactly 1 species
            if len(set(elements)) == 1:
                projector_label.append(elements[0])
            else:
                projector_label.append(None)

            # check that this projector has exactly 1 ang mom channel
            if len(set(ang_mom_channels)) == 1:
                projector_label.append(ang_mom_channels[0])
            else:
                projector_label.append(None)

            # check that this projector has exactly 1 spin channel
            if len(set(spin_channels)) == 1:
                projector_label.append(spin_channels[0])
            else:
                projector_label.append(None)

            projector_label = tuple(projector_label)
            projectors.append(projector_label)

    return projectors


def get_seed_metadata(doc, seed):
    """For a given document and seedname, look for
    ICSD CollCode, MaterialsProject IDs and DOIs to add
    to the document.

    Parameters:
        doc (dict): the input document.
        seed (str): the filename that is being scraped (sans file extension).

    """
    errors = []
    if "-CollCode-" in seed:
        try:
            doc["icsd"] = int(
                seed.split("CollCode-")[-1].split("-")[0].split("_")[0].split(".")[0]
            )
        except ValueError as exc:
            errors.append(exc)
    elif "CollCode" in seed:
        try:
            doc["icsd"] = int(
                seed.split("CollCode")[-1].split("-")[0].split("_")[0].split(".")[0]
            )
        except ValueError as exc:
            errors.append(exc)
    elif "-ICSD-" in seed:
        try:
            doc["icsd"] = int(
                seed.split("-ICSD-")[-1].split("-")[0].split("_")[0].split(".")[0]
            )
        except ValueError as exc:
            errors.append(exc)
    if "-OQMD-" in seed:
        try:
            doc["oqmd_id"] = int(
                seed.split("-OQMD-")[-1].split("-")[0].split("_")[0].split(".")[0]
            )
        except ValueError as exc:
            errors.append(exc)
    elif "-OQMD_" in seed:
        try:
            doc["oqmd_id"] = int(
                seed.split("-OQMD_")[-1].split("-")[0].split("_")[0].split(".")[0]
            )
        except ValueError as exc:
            errors.append(exc)
    if "-MP-" in seed:
        try:
            doc["mp_id"] = int(
                seed.split("-MP-")[-1].split("-")[0].split("_")[0].split(".")[0]
            )
        except ValueError as exc:
            errors.append(exc)
    elif "-MP_" in seed:
        try:
            doc["mp_id"] = int(
                seed.split("-MP_")[-1].split("-")[0].split("_")[0].split(".")[0]
            )
        except ValueError as exc:
            errors.append(exc)
    if "-DOI-" in seed:
        try:
            doc["doi"] = seed.split("-DOI-")[-1].split("-")[0].replace("__", "/")
        except ValueError as exc:
            errors.append(exc)

    if errors:
        raise ValueError(errors)
