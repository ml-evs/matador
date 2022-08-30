# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the scraper functions for the CIF
(Crystallographic Information File) format.

"""

import functools
from pathlib import Path
import numpy as np
from matador.scrapers.utils import scraper_function, get_flines_extension_agnostic
from matador.utils.cell_utils import (
    get_spacegroup_spg,
    abc2cart,
    cart2volume,
    frac2cart,
)

EPS = 1e-13


@scraper_function
def cif2dict(fname, **kwargs):
    """Extract available information from  .cif file and store as a
    dictionary. Raw cif data is stored under the `'_cif'` key. Symmetric
    sites are expanded by the symmetry operations and their occupancies
    are tracked.

    Parameters:
        fname (str/list): filename or list of filenames of .cif file(s)
            (with or without extension).

    Returns:
        (dict/str, bool): if successful, a dictionary containing scraped
            data and True, if not, then an error string and False.


    """

    flines, fname = get_flines_extension_agnostic(fname, "cif")

    doc = dict()
    cif_dict = _cif_parse_raw(flines)

    doc["_cif"] = cif_dict
    doc["source"] = [str(Path(fname).resolve())]

    doc["atom_types"] = []
    atom_labels = cif_dict.get("_atom_site_type_symbol", False)
    if not atom_labels:
        atom_labels = cif_dict.get("_atom_site_label", False)
    if not atom_labels:
        raise RuntimeError(f"Unable to find atom types in cif file {fname}.")

    for atom in atom_labels:
        symbol = ""
        for character in atom:
            if not character.isalpha():
                break
            else:
                symbol += character
        doc["atom_types"].append(symbol)

    doc["positions_frac"] = [
        list(map(lambda x: float(x.split("(")[0]), vector))
        for vector in zip(
            cif_dict["_atom_site_fract_x"],
            cif_dict["_atom_site_fract_y"],
            cif_dict["_atom_site_fract_z"],
        )
    ]

    if "_atom_site_occupancy" in cif_dict:
        doc["site_occupancy"] = [
            float(x.split("(")[0]) for x in cif_dict["_atom_site_occupancy"]
        ]
    else:
        doc["site_occupancy"] = [1.0 for _ in doc["positions_frac"]]

    if "_atom_site_symmetry_multiplicity" in cif_dict:
        doc["site_multiplicity"] = [
            float(x.split("(")[0]) for x in cif_dict["_atom_site_symmetry_multiplicity"]
        ]
    else:
        doc["site_multiplicity"] = [1.0 for _ in doc["positions_frac"]]

    doc["lattice_abc"] = [
        list(
            map(
                _cif_parse_float_with_errors,
                [
                    cif_dict["_cell_length_a"],
                    cif_dict["_cell_length_b"],
                    cif_dict["_cell_length_c"],
                ],
            )
        ),
        list(
            map(
                _cif_parse_float_with_errors,
                [
                    cif_dict["_cell_angle_alpha"],
                    cif_dict["_cell_angle_beta"],
                    cif_dict["_cell_angle_gamma"],
                ],
            )
        ),
    ]

    doc["lattice_cart"] = abc2cart(doc["lattice_abc"])
    doc["cell_volume"] = cart2volume(doc["lattice_cart"])
    doc["stoichiometry"] = _cif_disordered_stoichiometry(doc)
    doc["num_atoms"] = len(doc["positions_frac"])

    if (
        "_space_group_symop_operation_xyz" in doc["_cif"]
        and "_symmetry_equiv_pos_as_xyz" not in doc["_cif"]
    ):
        doc["_cif"]["_symmetry_equiv_pos_as_xyz"] = doc["_cif"][
            "_space_group_symop_operation_xyz"
        ]
    if "_symmetry_equiv_pos_as_xyz" in doc["_cif"]:
        _cif_set_unreduced_sites(doc)

    try:
        doc["space_group"] = get_spacegroup_spg(doc, check_occ=False)
    except RuntimeError:
        pass

    return doc, True


def _cif_parse_float_with_errors(x):
    """Strip bracketed errors from end of float."""
    return float(x.split("(")[0])


def _cif_disordered_stoichiometry(doc):
    """Create a matador stoichiometry normalised to the smallest integer
    number of atoms, unless all occupancies are 1/0.

    Parameters:
        doc: dictionary containing `atom_types`, `site_occupancy` and
            `site_multiplicity` keys.

    Returns:
        list of tuples: a standard matador stoichiometry.

    """
    from collections import defaultdict

    stoich = defaultdict(float)
    eps = 1e-8
    disordered = False
    for ind, site in enumerate(doc["atom_types"]):
        stoich[site] += doc["site_occupancy"][ind] * doc["site_multiplicity"][ind]
        if doc["site_multiplicity"][ind] % 1 > 1e-5:
            disordered = True

    if disordered:
        min_int = 1e10
        for atom in stoich:
            if abs(int(stoich[atom]) - stoich[atom]) < eps:
                if int(stoich[atom]) < min_int:
                    min_int = int(stoich[atom])

        if min_int == 1e10:
            min_int = 1
        for atom in stoich:
            stoich[atom] /= min_int

    return sorted([[atom, stoich[atom]] for atom in stoich])


def _cif_parse_raw(flines):
    """Parse raw CIF file data into a dictionary.

    Parameters:
        flines (:obj:`list` of :obj:`str`): contents of .cif file.

    Returns:
        dict: dictionary containing cif data with native fields/ordering.

    """
    ind = 0
    cif_dict = dict()
    cif_dict["loops"] = list()
    while ind < len(flines):
        jnd = 1
        line = flines[ind].strip()
        # parse single (multi-line) tag
        if line.startswith("_"):
            line = line.split()
            key = line[0]
            data = ""
            if len(line) > 1:
                data += " ".join(line[1:])
            while ind + jnd < len(flines) and _cif_line_contains_data(
                flines[ind + jnd].strip()
            ):
                data += flines[ind + jnd].strip().replace(";", "")
                jnd += 1
            cif_dict[key] = data.strip()
        # parse loop block
        elif line.startswith("loop_"):
            # get loop keys
            keys = []
            while flines[ind + jnd].strip().startswith("_"):
                keys.append(flines[ind + jnd].strip())
                jnd += 1
            for key in keys:
                cif_dict[key] = []
            cif_dict["loops"].append(keys)
            while ind + jnd < len(flines) and _cif_line_contains_data(
                flines[ind + jnd].strip()
            ):
                data = ""
                # loop over line and next lines
                while ind + jnd < len(flines) and _cif_line_contains_data(
                    flines[ind + jnd]
                ):
                    data += flines[ind + jnd]
                    jnd += 1

            loop_dict = _cif_parse_loop(keys, data)
            cif_dict.update(loop_dict)

        ind += jnd

    return cif_dict


def _cif_parse_loop(keys, data_block):
    """A hacky way to parse CIF data loops that can be split by quotes
    or spaces. There must be a better way...

    Parameters:
        keys (list of str): list of keys for the loop.
        data_block (str): raw string of the entire data block.

    Returns:
        Dict[str, str]: a dictionary with keys from ``keys``, containing the
            data split by quotes and spaces. All data is left as
            strings for further processing.

    """

    from collections import deque, defaultdict

    dq = deque(data_block)
    data_list = []
    entry = None
    in_quotes = False
    while dq:
        char = dq.popleft()
        if not char.strip() and entry is None:
            continue
        elif (
            (not char.strip() or char in [" ", ";"])
            and entry is not None
            and not in_quotes
        ):
            data_list.append(entry.strip())
            entry = None
        elif not char.strip() and entry is not None and in_quotes:
            entry += " "
        elif char == "'" and entry and entry is not None:
            in_quotes = False
            data_list.append(entry.strip())
            entry = None
        elif char == "'" and entry is None:
            entry = ""
            in_quotes = True
        else:
            if entry is None:
                entry = char
            else:
                entry += char

    loop_dict = defaultdict(list)
    for ind, entry in enumerate(data_list):
        ind = ind % len(keys)
        loop_dict[keys[ind]].append(entry)

    return loop_dict


def _cif_set_unreduced_sites(doc):
    """Expands sites by symmetry operations found under the key
    `symemtry_equiv_pos_as_xyz` in the cif_dict.

    Parameters:
        doc (dict): matador document to modify. Must contain symops
            under doc['_cif']['_symmetry_equiv_pos_as_xyz']. This doc
            is updated with new `positions_frac`, `num_atoms`, `atom_types`
            and `site_occupancy`.

    """
    from matador.utils.cell_utils import wrap_frac_coords
    from matador.utils.cell_utils import calc_pairwise_distances_pbc
    from matador.fingerprints.pdf import PDF

    species_sites = dict()
    species_occ = dict()

    symmetry_ops = []
    symmetry_functions = []

    def _apply_sym_op(x=None, y=None, z=None, symmetry=None):
        """Returns the site after the applied symmetry operation, in string representation."""
        # cannot use a listcomp here due to interplay with functools
        return [eval(symmetry[0]), eval(symmetry[1]), eval(symmetry[2])]

    for symmetry in doc["_cif"]["_symmetry_equiv_pos_as_xyz"]:
        symmetry = tuple(elem.strip() for elem in symmetry.strip("'").split(","))
        # check the element before doing an eval, as it is so unsafe
        allowed_chars = [
            "x",
            "y",
            "z",
            ".",
            "/",
            "+",
            "-",
            "0",
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
        ]
        for element in symmetry:
            for character in element:
                if character not in allowed_chars:
                    raise RuntimeError(
                        "You are trying to do something naughty with the symmetry element {}".format(
                            element
                        )
                    )

        symmetry_ops.append(symmetry)
        symmetry_functions.append(functools.partial(_apply_sym_op, symmetry=symmetry))

    for ind, site in enumerate(doc["positions_frac"]):
        species = doc["atom_types"][ind]
        occupancy = doc["site_occupancy"][ind]
        if doc["atom_types"][ind] not in species_sites:
            species_sites[species] = []
            species_occ[species] = []
        for symmetry in symmetry_functions:
            x, y, z = site
            new_site = symmetry(x=x, y=y, z=z)
            new_site = wrap_frac_coords([new_site])[0]
            species_sites[species].append(new_site)
            species_occ[species].append(occupancy)

    unreduced_sites = []
    unreduced_occupancies = []
    unreduced_species = []

    # this loop assumes that no symmetry operation can map 2 unlike sites upon one another
    for species in species_sites:
        unreduced_sites.extend(species_sites[species])
        unreduced_occupancies.extend(species_occ[species])
        unreduced_species.extend(len(species_sites[species]) * [species])

    # check that the symmetry procedure has not generated overlapping atoms
    # this can happen for certain symmetries/cells if positions are not
    # reported to sufficient precision
    images = PDF._get_image_trans_vectors_auto(
        doc["lattice_cart"],
        0.1,
        0.01,
        max_num_images=1,
    )

    poscarts = frac2cart(doc["lattice_cart"], unreduced_sites)
    distances = calc_pairwise_distances_pbc(
        poscarts,
        images,
        doc["lattice_cart"],
        0.01,
        compress=False,
        filter_zero=False,
        per_image=True,
    )

    dupe_set = set()
    for img in distances:
        try:
            i_s, j_s = np.where(~img.mask)
        except ValueError:
            # ValueError will be raised if there is only one atom as i_s, j_s cannot be unpacked
            continue
        for i, j in zip(i_s, j_s):
            if i == j:
                continue
            else:
                # sites can overlap if they have partial occupancy
                if i not in dupe_set and unreduced_species[i] == unreduced_species[j]:
                    dupe_set.add(j)

    doc["positions_frac"] = unreduced_sites
    doc["site_occupancy"] = unreduced_occupancies
    doc["atom_types"] = unreduced_species

    doc["site_occupancy"] = [
        atom for ind, atom in enumerate(unreduced_occupancies) if ind not in dupe_set
    ]
    doc["atom_types"] = [
        atom for ind, atom in enumerate(unreduced_species) if ind not in dupe_set
    ]
    doc["positions_frac"] = [
        atom for ind, atom in enumerate(unreduced_sites) if ind not in dupe_set
    ]

    _num_atoms = np.sum(doc["site_occupancy"])
    if abs(_num_atoms - round(_num_atoms, 0)) < EPS:
        _num_atoms = int(round(_num_atoms, 0))
    doc["num_atoms"] = _num_atoms

    if len(doc["site_occupancy"]) != len(doc["positions_frac"]):
        raise RuntimeError(
            "Size mismatch between positions and occs, {} vs {}".format(
                len(doc["site_occupancy"]), len(doc["positions_frac"])
            )
        )
    if len(doc["positions_frac"]) != len(doc["atom_types"]):
        raise RuntimeError("Size mismatch between positions and types")


def _cif_line_contains_data(line):
    """Check if string contains cif-style data."""
    return not any(
        [line.startswith("_"), line.startswith("#"), line.startswith("loop_")]
    )


@scraper_function
def _ase_cif2dict(fname):
    """Read cif file into ASE object,
    then convert ASE Atoms into matador document.

    Parameters:
        fname (str): cif filename

    Returns:
        (dict, bool): simple matador document with error status.

    """
    import ase.io
    from matador.utils.ase_utils import ase2dict

    fname = fname.replace(".cif", "")
    atoms = ase.io.read(fname + ".cif")
    doc = ase2dict(atoms)

    return doc, True
