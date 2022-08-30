# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements some useful functions for real/reciprocal
cell manipulation, symmetry checking and sampling (e.g. grids and paths.)

"""

from __future__ import annotations
from typing import Dict, Any, Union, Tuple, List, TYPE_CHECKING

import numpy as np

from matador.data.periodic_table import PERIODIC_TABLE

if TYPE_CHECKING:
    from matador.crystal import Crystal


EPS = 1e-12


def abc2cart(lattice_abc: List[List[float]]) -> List[List[float]]:
    """Converts lattice parameters into Cartesian lattice vectors.

    Parameters:
        lattice_abc: The lattice parameters [[a, b, c], [alpha, beta, gamma]]

    Returns:
        Cartesian lattice vectors.

    """
    assert len(lattice_abc) == 2
    assert len(lattice_abc[0]) == 3
    assert len(lattice_abc[1]) == 3
    a = lattice_abc[0][0]
    b = lattice_abc[0][1]
    c = lattice_abc[0][2]
    deg2rad = np.pi / 180
    alpha = lattice_abc[1][0] * deg2rad
    beta = lattice_abc[1][1] * deg2rad
    gamma = lattice_abc[1][2] * deg2rad
    lattice_cart = []
    lattice_cart.append([a, 0.0, 0.0])
    bx = b * np.cos(gamma)
    by = b * np.sin(gamma)
    tol = 1e-12
    if abs(bx) < tol:
        bx = 0.0
    if abs(by) < tol:
        by = 0.0
    cx = c * np.cos(beta)
    if abs(cx) < tol:
        cx = 0.0
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    if abs(cy) < tol:
        cy = 0.0
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    lattice_cart.append([bx, by, 0.0])
    lattice_cart.append([cx, cy, cz])
    return lattice_cart


def cart2abcstar(lattice_cart: List[List[float]]) -> np.ndarray:
    """Convert Cartesian lattice vectors to the reciprocal of the lattice vectors,
    NOT the reciprocal lattice vectors (useful when writing PDB files).

    Parameters:
        lattice_cart: Cartesian lattice vectors.

    Returns:
        Reciprocal of the lattice vectors.

    """
    return np.asarray(real2recip(lattice_cart)) / (2 * np.pi)


def cart2volume(lattice_cart: List[List[float]]) -> float:
    """Convert lattice_cart to cell volume.

    Parameters:
        lattice_cart: Cartesian lattice vectors.

    Returns:
        Cell volume in the same unit as the lattice vectors (cubed).

    """
    lattice_cart = np.asarray(lattice_cart)
    vol = np.abs(np.dot(np.cross(lattice_cart[0], lattice_cart[1]), lattice_cart[2]))
    return vol


def cart2abc(lattice_cart: List[List[float]]) -> List[List[float]]:
    """Convert Cartesian lattice vectors to lattice parametres.

    Parameters:
        lattice_cart: Cartesian lattice vectors.

    Returns:
        The lattice parameters :math:`[(a, b, c), (\\alpha, \\beta, \\gamma)]`.

    """
    vecs = lattice_cart
    lattice_abc = []
    a, b, c = (np.sqrt(sum([val**2 for val in vec])) for vec in vecs)
    lattice_abc.append([a, b, c])
    # np.cos(alpha) = b.c /|b * c|
    radians2deg = 180.0 / np.pi
    cos_alpha = sum([val_b * val_c for (val_c, val_b) in zip(vecs[2], vecs[1])]) / (
        b * c
    )
    cos_beta = sum([val_c * val_a for (val_c, val_a) in zip(vecs[2], vecs[0])]) / (
        c * a
    )
    cos_gamma = sum([val_a * val_b for (val_a, val_b) in zip(vecs[0], vecs[1])]) / (
        a * b
    )
    alpha = radians2deg * np.arccos(cos_alpha)
    beta = radians2deg * np.arccos(cos_beta)
    gamma = radians2deg * np.arccos(cos_gamma)
    lattice_abc.append([alpha, beta, gamma])
    return lattice_abc


def frac2cart(
    lattice_cart: List[List[float]], positions_frac: List[List[float]]
) -> List[List[float]]:
    """Convert positions_frac block into positions_abs.

    Parameters:
        lattice_cart: Cartesian lattice vectors.
        positions_frac: List of fractional position vectors.

    Returns:
        List of absolute position vectors.

    """
    _positions_frac = np.asarray(positions_frac)
    reshaped = False
    if len(np.shape(_positions_frac)) == 1:
        reshaped = True
        _positions_frac = _positions_frac.reshape((1, 3))
    _lattice_cart = np.asarray(lattice_cart)
    positions_abs = switch_coords(_lattice_cart, _positions_frac)
    if reshaped:
        positions_abs = positions_abs.reshape(-1)
    return positions_abs.tolist()


def wrap_frac_coords(
    positions: Union[List[List[float]], List[float]], remove: bool = False
) -> Union[List[List[float]], List[float]]:
    """Wrap the given fractional coordinates back into the cell.

    Parameters:
        positions: list of fractional position vectors, or
            a single position.

    Keyword arguments:
        remove: if True, removes points exterior to the cell.

    Returns:
        List of wrapped fractional position vectors.

    """
    from copy import deepcopy

    wrapped = deepcopy(positions)
    single = False
    if len(wrapped) == 3 and isinstance(wrapped[0], float):
        wrapped = [wrapped]
        single = True
    if remove:
        to_remove = len(wrapped) * [False]
    for ind, pos in enumerate(wrapped):
        for k in range(3):
            if pos[k] >= 1 or pos[k] < 0:
                if remove:
                    to_remove[ind] = True
                    break
                else:
                    wrapped[ind][k] %= 1
    if remove:
        wrapped = [wrapped[ind] for ind, res in enumerate(to_remove) if not res]
    if single:
        return wrapped[0]

    return wrapped


def switch_coords(
    lattice: np.ndarray, pos: np.ndarray, norm: float = None
) -> np.ndarray:
    """Act on coordinates with the relevant lattice
    vectors to switch from fractional to absolute coordinates.

    Parameters:

        lattice: either lattice_cart or reciprocal lattice_cart (3x3 array)
        pos: input positions to convert (3xN array for N atoms).

    Keyword arguments:
        norm (float): divide final coordinates by normalisation factor, e.g.
            :math:`2 \\pi` when lattice is recip and positions are cartesian.

    Returns:
        3xN array of converted positions.

    """
    new_pos = np.zeros_like(pos)
    for ni, _ in enumerate(pos):
        for j in range(3):
            new_pos[ni] += lattice[j] * pos[ni][j]
    if norm is not None:
        new_pos /= norm
    return new_pos


def cart2frac(
    lattice_cart: List[List[float]], positions_abs: List[List[float]]
) -> List[List[float]]:
    """Convert positions_abs block into positions_frac (and equivalent
    in reciprocal space).

    Parameters:
        lattice_cart: Cartesian lattice vectors.
        positions_abs: list of absolute position vectors.

    Returns:
        List of fractional position vectors with the same shape
            as the input list.

    """
    _positions_abs = np.asarray(positions_abs, dtype=np.float64)
    reshaped = False
    if len(np.shape(_positions_abs)) == 1:
        reshaped = True
        _positions_abs = _positions_abs.reshape((1, 3))
    recip_lat = np.asarray(real2recip(lattice_cart))
    recip_lat = recip_lat.T
    positions_frac = switch_coords(recip_lat, _positions_abs, norm=2 * np.pi)
    if reshaped:
        positions_frac = positions_frac.reshape(-1)
    return positions_frac.tolist()


def real2recip(real_lat: List[List[float]]) -> List[List[float]]:
    """Convert the real lattice in Cartesian basis to
    the reciprocal space lattice.

    Parameters:
        real_lat: Cartesian lattice vectors.

    Returns:
        Cartesian lattice vectors of reciprocal lattice.

    """
    real_lat = np.asarray(real_lat)
    recip_lat = np.zeros((3, 3))
    volume = np.dot(real_lat[0], np.cross(real_lat[1], real_lat[2]))
    recip_lat[0] = (2 * np.pi) * np.cross(real_lat[1], real_lat[2])
    recip_lat[1] = (2 * np.pi) * np.cross(real_lat[2], real_lat[0])
    recip_lat[2] = (2 * np.pi) * np.cross(real_lat[0], real_lat[1])
    recip_lat /= volume
    return recip_lat.tolist()


def calc_mp_grid(
    lattice_cart: List[List[float]], spacing: float
) -> Tuple[int, int, int]:
    """Return the Monkhorst-Pack grid based on lattice
    vectors and desired spacing.

    Parameters:
        lattice_cart: Cartesian lattice vectors.
        spacing: desired maximum grid spacing.

    Returns:
        List of 3 integers defining the MP grid.

    """
    recip_lat = real2recip(lattice_cart)
    recip_len = np.zeros((3))
    recip_len = np.sqrt(np.sum(np.power(recip_lat, 2), axis=1))
    mp_grid = recip_len / (2 * np.pi * spacing)
    return [int(np.ceil(elem)) for elem in mp_grid]


def shift_to_include_gamma(mp_grid: Tuple[int, int, int]) -> Tuple[float, float, float]:
    """Calculate the shift required to include the :math:`\\Gamma`-point
    in the Monkhorst-Pack grid.

    Parameters:
        mp_grid: number of grid points
            in each reciprocal space direction.

    Returns:
        The shift required to include the :math:`\\Gamma`.

    """
    shift = [0, 0, 0]
    for ind, val in enumerate(mp_grid):
        if val % 2 == 0:
            shift[ind] = 1.0 / (val * 2)
    return shift


def shift_to_exclude_gamma(mp_grid: Tuple[int, int, int]) -> Tuple[float, float, float]:
    """Calculate the shift required to exclude the :math:`\\Gamma`-point
    in the Monkhorst-Pack grid. Returns the "minimal shift", i.e. only
    one direction will be shifted.

    Parameters:
        mp_grid: number of grid points
            in each reciprocal space direction.

    Returns:
        The shift required to exclude :math:`\\Gamma`.

    """
    shift = [0, 0, 0]
    if all([val % 2 == 1 for val in mp_grid]):
        for ind, val in enumerate(mp_grid):
            if val % 2 == 1:
                shift[ind] = 1.0 / (val * 2)
                break

    return shift


def get_best_mp_offset_for_cell(doc: Union[Dict[str, Any], Crystal]) -> List[float]:
    """Calculates the "best" kpoint_mp_offset to use for the passed
    cell. If the crystal has a hexagonal space group, then the offset
    returned will shift the grid to include :math:`\\Gamma`-point, and vice
    versa for non-hexagonal cells.

    Parameters:
        doc: matador document/Crystal to consider, containing structural
            information and a "kpoints_mp_spacing" key.

    Returns:
        The desired `kpoint_mp_offset`.

    """

    gamma = False
    if "6" in get_spacegroup_spg(doc, symprec=1e-5):
        gamma = True

    if "lattice_cart" not in doc or "kpoints_mp_spacing" not in doc:
        raise RuntimeError("Unable to calculate offset without lattice or spacing")

    mp_grid = calc_mp_grid(doc["lattice_cart"], doc["kpoints_mp_spacing"])
    if gamma:
        return shift_to_include_gamma(mp_grid)

    return shift_to_exclude_gamma(mp_grid)


def calc_mp_spacing(
    real_lat: List[List[float]], mp_grid: Tuple[int, int, int], prec: int = 3
) -> float:
    """Convert real lattice in Cartesian basis and the CASTEP
    `kpoint_mp_grid` into a CASTEP grid spacing.

    Parameters:
        real_lat: Cartesian lattice vectors.
        mp_grid: 3 integers defining the MP grid.

    Keyword arguments:
        prec: desired decimal precision of output.

    Returns:
        `kpoint_mp_spacing` rounded to `prec` decimal places.

    """
    recip_lat = real2recip(real_lat)
    recip_len = np.zeros((3))
    recip_len = np.sqrt(np.sum(np.power(recip_lat, 2), axis=1))
    spacing = recip_len / (2 * np.pi * np.asarray(mp_grid))
    max_spacing = np.max(spacing)
    exponent = round(np.log10(max_spacing) - prec)
    return round(max_spacing + 0.5 * 10**exponent, prec)


def get_seekpath_kpoint_path(
    doc: Union[Dict[str, Any], Tuple],
    standardize: bool = True,
    explicit: bool = True,
    spacing: float = 0.01,
    threshold: float = 1e-7,
    debug: bool = False,
    symmetry_tol: float = None,
) -> Tuple[Dict[str, Any], List[List[float]], Dict[str, Any]]:
    """Return the conventional kpoint path of the relevant crystal system
    according to the definitions by "HKPOT" in
    Comp. Mat. Sci. 128, 2017:

    http://dx.doi.org/10.1016/j.commatsci.2016.10.015

    Parameters:
        doc: matador doc or spglib tuple to find kpoint path for.

    Keyword arguments:
        spacing: desired kpoint spacing
        threshold: internal seekpath threshold
        symmetry_tol: spglib symmetry tolerance

    Returns:
        dict: standardized version of input doc
        list: list of kpoint positions
        dict: full dictionary of all seekpath results

    """
    try:
        from seekpath import get_explicit_k_path, get_path
    except ImportError:
        raise ImportError(
            "SeeK-Path dependency missing, please install it with `pip install seekpath`."
        )

    if symmetry_tol is None:
        symmetry_tol = 1e-5

    if isinstance(doc, tuple):
        spg_structure = doc
    else:
        if standardize:
            spg_structure = doc2spg(standardize_doc_cell(doc, symprec=symmetry_tol))
        else:
            spg_structure = doc2spg(doc)

    if explicit:
        seekpath_results = get_explicit_k_path(
            spg_structure,
            reference_distance=spacing,
            with_time_reversal=True,
            symprec=symmetry_tol,
            threshold=threshold,
        )

        kpt_path = seekpath_results["explicit_kpoints_rel"]
    else:
        seekpath_results = get_path(spg_structure)
        kpt_path = []

    primitive_doc = dict()
    primitive_doc["lattice_cart"] = seekpath_results["primitive_lattice"]
    primitive_doc["positions_frac"] = seekpath_results["primitive_positions"]
    primitive_doc["atom_types"] = [
        list(PERIODIC_TABLE)[i] for i in seekpath_results["primitive_types"]
    ]
    primitive_doc["num_atoms"] = len(primitive_doc["atom_types"])
    primitive_doc["lattice_abc"] = cart2abc(primitive_doc["lattice_cart"])
    primitive_doc["cell_volume"] = cart2volume(primitive_doc["lattice_cart"])
    if debug:
        print(
            "Found lattice type {}".format(seekpath_results["bravais_lattice_extended"])
        )
        print("Old lattice:\n", np.asarray(doc["lattice_cart"]))
        print("Contained {} atoms".format(doc["num_atoms"]))
        print("New lattice:\n", np.asarray(primitive_doc["lattice_cart"]))
        print("Contains {} atoms".format(primitive_doc["num_atoms"]))
        print("k-point path contains {} points.".format(len(kpt_path)))

    if "site_occupancy" in doc:
        if min(doc["site_occupancy"]) < 1 - EPS:
            print("Ignoring any site occupancy found in this cell.")
        primitive_doc["site_occupancy"] = [1 for atom in primitive_doc["atom_types"]]

    return primitive_doc, kpt_path, seekpath_results


def doc2spg(
    doc: Union[Dict[str, Any], "Crystal"], check_occ: bool = True
) -> Tuple[List[List[float]], List[List[float]], List[int]]:
    """Return an spglib input tuple from a matador doc.

    Parameters:
        doc: matador document or Crystal object.

    Keyword arguments:
        check_occ: check for partial occupancy and raise an error if present.

    Returns:
        tuple: spglib-style tuple of lattice, positions and types.

    """
    from matador.utils.chem_utils import get_atomic_number

    try:
        if "lattice_cart" not in doc:
            doc["lattice_cart"] = abc2cart(doc["lattice_abc"])
    except KeyError:
        raise RuntimeError("doc2spg failed: unable to find lattice in document.")

    required_keys = ("lattice_cart", "positions_frac", "atom_types")
    if all(key in doc for key in required_keys):
        if "lattice_cart" not in doc:
            doc["lattice_cart"] = abc2cart(doc["lattice_abc"])
        cell = (
            doc["lattice_cart"],
            doc["positions_frac"],
            [get_atomic_number(elem) for elem in doc["atom_types"]],
        )
        if check_occ and np.min(doc.get("site_occupancy", [1.0])) < 1.0:
            raise RuntimeError("spglib does not support partial occupancy.")

        return cell

    else:
        raise RuntimeError(
            f"Unable to use doc2spg, one of {required_keys} was missing."
        )


def get_space_group_label_latex(label: str) -> str:
    """Return the LaTeX format of the passed space group label. Takes
    any string, leaves the first character upright, italicses the rest,
    handles subscripts and bars over numbers.

    Parameters:
        label: a given space group in "standard" plain text format,
        e.g. P-63m to convert to '$P\\bar{6}3m$'.

    Returns:
        The best attempt to convert the label to LaTeX format.

    """
    import re

    return "${}$".format(re.sub("-(?P<number>[0-9])", "\\\\bar{\\g<number>}", label))


def standardize_doc_cell(
    doc: Union["Crystal", Dict[str, Any]], primitive: bool = True, symprec: float = 1e-2
) -> Union["Crystal", Dict[str, Any]]:
    """Return standardized cell data from matador doc.

    Parameters:
        doc: matador document or Crystal object to standardize.

    Keyword arguments:
        primitive: whether to reduce cell to primitive.
        symprec: spglib symmetry tolerance.

    Returns:
        A matador document/`Crystal` containing standardized cell.

    """
    import spglib as spg
    from matador.utils.chem_utils import get_atomic_symbol
    from copy import deepcopy
    from matador.crystal import Crystal

    spg_cell = doc2spg(doc)
    spg_standardized = spg.standardize_cell(
        spg_cell, to_primitive=primitive, symprec=symprec
    )
    if not isinstance(doc, Crystal):
        std_doc = deepcopy(doc)
    else:
        std_doc = deepcopy(doc._data)
    std_doc["lattice_cart"] = [list(vec) for vec in spg_standardized[0]]
    std_doc["lattice_abc"] = cart2abc(std_doc["lattice_cart"])
    std_doc["positions_frac"] = [list(atom) for atom in spg_standardized[1]]
    std_doc["atom_types"] = [get_atomic_symbol(atom) for atom in spg_standardized[2]]
    std_doc["site_occupancy"] = len(std_doc["positions_frac"]) * [1]
    std_doc["cell_volume"] = cart2volume(std_doc["lattice_cart"])
    std_doc["space_group"] = get_spacegroup_spg(std_doc, symprec=symprec)
    # if the original document was a crystal, return a new one
    if isinstance(doc, Crystal):
        std_doc = Crystal(std_doc)

    return std_doc


def get_spacegroup_spg(
    doc: Union[Dict[str, Any], Crystal], symprec: float = 0.01, check_occ: bool = True
):
    """Return spglib spacegroup for a cell.

    Parameters:
        doc: matador document or Crystal object.

    Keyword arguments:
        symprec: spglib symmetry tolerance.

    Returns:
        The H-M space group symbol of structure.

    """
    import spglib as spg

    spg_cell = doc2spg(doc, check_occ=check_occ)
    space_group = spg.get_spacegroup(spg_cell, symprec=symprec)
    if space_group is None:
        raise RuntimeError("Spglib was unable to calculate space group.")

    return space_group.split(" ")[0]


def get_compatible_spacegroups(
    doc: Union[Dict[str, Any], Crystal], symprec_range=(-5, 0)
) -> Dict[float, str]:
    """Return the space group of a given crystal for the range of symprecs.

    Parameters:
        doc: The crystal to analyse.
        symprec_range: The range of symprecs to test (log space).

    Returns:
        A mapping from symprec to space group symbol.

    """
    spgs = {}
    for symprec in np.logspace(*symprec_range, num=10):
        spgs[symprec] = get_spacegroup_spg(doc, symprec=symprec)

    return spgs


def add_noise(doc: Dict[str, Any], amplitude: float = 0.1) -> Dict[str, Any]:
    """Add random noise to the positions of structure contained in doc.
    Useful for force convergence tests.

    Parameters:
        doc: dictionary containing matador structure.

    Keyword arguments:
        amplitude: maximum amplitude of noise vector.

    Raises:
        KeyError if (`lattice_cart` and `positions_frac`) or `positions_abs`
            are missing.

    Returns:
        The randomised structure.
    """
    poscart = np.asarray(
        doc.get("positions_abs")
        or frac2cart(doc["lattice_cart"], doc["positions_frac"])
    )
    for atom in poscart:
        noise_vector = np.random.rand(3)
        noise_vector /= np.sqrt(np.sum(noise_vector**2))
        noise_vector *= amplitude * (np.random.rand())
        atom += noise_vector
    doc["positions_abs"] = poscart.tolist()
    doc["positions_frac"] = cart2frac(doc["lattice_cart"], doc["positions_abs"])

    return doc


def calc_pairwise_distances_pbc(
    poscart,
    images,
    lattice,
    rmax,
    poscart_b=None,
    compress=False,
    debug=False,
    filter_zero=False,
    per_image=False,
):
    """Calculate PBC distances with SciPy's cdist, given the
    image cell vectors.

    Parameters:
        poscart (numpy.ndarray): list or array of absolute atomic coordinates.
        images: iterable of lattice vector multiples (e.g. [2, -1, 3])
            required to obtain the translation to desired image cells.
        lattice (:obj:`list` if :obj:`list`): list of lattice vectors of
            the real cell.
        rmax (float): maximum value after which to mask the array.

    Keyword arguments:
        poscart_b (numpy.ndarray): absolute positions of another type of
            atom, where only A-B distances will be calculated.
        debug (bool): print timing data and how many distances were masked.
        compress (bool): whether or not to compressed the output array,
            useful when e.g. creating PDFs but not when atom ID is important.
        filter_zero (bool): whether or not to filter out the "self-interaction"
            zero distances.
        per_image (bool): return a list of distances per image, as opposed to
            one large flat. This preserves atom IDs for use elsewhere.

    Returns:
        distances (numpy.ndarray): pairwise 2-D d_ij masked array with values
            or stripped 1-D array containing just the distances, or a list of
            numpy arrays if per_image is True.

    """
    from scipy.spatial.distance import cdist
    import time

    if debug:
        start = time.time()
    _lattice = np.asarray(lattice)
    _poscart = np.asarray(poscart)
    image_distances = []
    if poscart_b is None:
        _poscart_b = _poscart
    else:
        _poscart_b = np.asarray(poscart_b)

    num_pairs = len(_poscart) * len(_poscart_b)

    if per_image:
        for image_ind, prod in enumerate(images):
            distances = cdist(_poscart, _poscart_b + prod @ _lattice)
            distances = np.ma.masked_where(distances > rmax, distances, copy=False)
            image_distances.append(distances)

        distances = image_distances

    else:
        array_size = len(_poscart) * len(images) * len(_poscart_b)
        distances = np.empty(array_size)
        for image_ind, prod in enumerate(images):
            distances[image_ind * num_pairs : (image_ind + 1) * num_pairs] = cdist(
                _poscart, _poscart_b + prod @ _lattice
            ).flatten()

        distances = np.ma.masked_where(distances > rmax, distances, copy=False)
        if filter_zero:
            distances = np.ma.masked_where(distances < EPS, distances, copy=False)

        if debug:
            print(
                "Calculated: {}, Used: {}, Ignored: {}".format(
                    len(distances),
                    np.ma.count(distances),
                    np.ma.count_masked(distances),
                )
            )
        if compress:
            distances = distances.compressed()

    if debug:
        end = time.time()
        print("Calculated distances in {} s".format(end - start))

    return distances


def create_simple_supercell(
    doc: Union[Dict[str, Any], "Crystal"],
    extension: Tuple[int, int, int],
    standardize: bool = False,
    symmetric: bool = False,
) -> Union[Dict[str, Any], "Crystal"]:
    """Return a document with new supercell, given extension vector.

    Parameters:
        doc: matador doc to construct cell from.
        extension: multiplicity of each lattice vector, e.g. (2,2,1).

    Keyword arguments:
        standardize: whether or not to use spglib to standardize the cell first.
        symmetric: whether or not centre the new cell on the origin.

    Returns:
        The supercell, either as a `Crystal` or as a dictionary,
        depending on the input type.

    """
    from itertools import product
    from copy import deepcopy
    from matador.crystal import Crystal

    if not all([elem >= 1 for elem in extension]) or extension == (1, 1, 1):
        raise RuntimeError("Bad/null supercell {} requested...".format(extension))

    if isinstance(doc, Crystal):
        supercell_doc = deepcopy(doc._data)
    else:
        supercell_doc = deepcopy(doc)

    # copy new doc and delete data that will not be corrected
    keys_to_copy = set(
        ["positions_frac", "lattice_cart", "atom_types", "source", "stoichiometry"]
    )
    supercell_doc = {key: deepcopy(doc[key]) for key in keys_to_copy}

    # standardize cell with spglib
    if standardize:
        supercell_doc = standardize_doc_cell(supercell_doc)

    for i, elem in enumerate(extension):
        supercell_doc["lattice_cart"] = np.asarray(supercell_doc["lattice_cart"])
        for k in range(3):
            supercell_doc["lattice_cart"][i][k] = (
                elem * supercell_doc["lattice_cart"][i][k]
            )
        for ind, atom in enumerate(supercell_doc["positions_frac"]):
            supercell_doc["positions_frac"][ind][i] /= elem

    images = product(*[list(range(elem)) for elem in extension])
    new_positions = []
    new_atoms = []
    for image in images:
        if image != (0, 0, 0):
            for ind, atom in enumerate(supercell_doc["atom_types"]):
                new_pos = deepcopy(supercell_doc["positions_frac"][ind])
                for i, elem in enumerate(image):
                    new_pos[i] += image[i] / extension[i]
                new_positions.append(new_pos)
                new_atoms.append(atom)
    supercell_doc["atom_types"].extend(new_atoms)
    supercell_doc["positions_frac"].extend(new_positions)

    if symmetric:
        supercell_doc["positions_frac"] = np.asarray(supercell_doc["positions_frac"])
        supercell_doc["positions_frac"] -= np.mean(
            supercell_doc["positions_frac"], axis=0
        )
        supercell_doc["positions_frac"] = supercell_doc["positions_frac"].tolist()

    supercell_doc["num_atoms"] = len(supercell_doc["atom_types"])
    supercell_doc["cell_volume"] = cart2volume(supercell_doc["lattice_cart"])
    supercell_doc["lattice_abc"] = cart2abc(supercell_doc["lattice_cart"])
    supercell_doc["source"] = doc["source"]
    supercell_doc["stoichiometry"] = doc["stoichiometry"]

    assert np.isclose(
        supercell_doc["cell_volume"],
        np.prod(extension) * cart2volume(doc["lattice_cart"]),
    )

    if isinstance(doc, Crystal):
        return Crystal(supercell_doc)

    return supercell_doc


def create_supercell_with_minimum_side_length(
    doc: Union[Dict[str, Any], "Crystal"], target: float
):
    """Pad the cell such that the minimum side length is greater than the target length.

    Parameters:
        doc: The crystal structure to pad.
        target: The target minimum side length.

    Returns:
        The supercell.

    """
    extension = [1, 1, 1]
    for ind, lat in enumerate(doc.cell.lengths):
        if lat < target:
            extension[ind] = int(-(target // -lat))

    return create_simple_supercell(doc, extension, standardize=False)
