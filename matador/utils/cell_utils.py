# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements some useful functions for real/reciprocal
cell manipulation, symmetry checking and sampling (e.g. grids and paths.)

"""

import numpy as np
from periodictable import elements
EPS = 1e-12


def abc2cart(lattice_abc):
    """ Converts lattice parameters into Cartesian lattice vectors.

    Parameters:
        lattice_abc (list): [[a, b, c], [alpha, beta, gamma]]

    Returns:
        list: Cartesian lattice vectors.

    """
    assert len(lattice_abc) == 2
    assert len(lattice_abc[0]) == 3
    assert len(lattice_abc[1]) == 3
    a = lattice_abc[0][0]
    b = lattice_abc[0][1]
    c = lattice_abc[0][2]
    deg2rad = np.pi/180
    alpha = lattice_abc[1][0] * deg2rad
    beta = lattice_abc[1][1] * deg2rad
    gamma = lattice_abc[1][2] * deg2rad
    lattice_cart = []
    lattice_cart.append([a, 0.0, 0.0])
    # vec(b) = (b np.cos(gamma), b np.sin(gamma), 0)
    bx = b*np.cos(gamma)
    by = b*np.sin(gamma)
    tol = 1e-12
    if abs(bx) < tol:
        bx = 0.0
    if abs(by) < tol:
        by = 0.0
    cx = c*np.cos(beta)
    if abs(cx) < tol:
        cx = 0.0
    cy = c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
    if abs(cy) < tol:
        cy = 0.0
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    lattice_cart.append([bx, by, 0.0])
    lattice_cart.append([cx, cy, cz])
    return lattice_cart


def cart2abcstar(lattice_cart):
    """ Convert lattice_cart =[[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]]
    to the reciprocal of the lattice vectors, NOT the reciprocal lattice vectors.

    Parameters:
        lattice_cart (list): Cartesian lattice vectors.

    Returns:
        list: lattice parameters [[a,b,c],[alpha,beta,gamma]]

    """
    return np.asarray(real2recip(lattice_cart)) / (2*np.pi)


def cart2volume(lattice_cart):
    """ Convert lattice_cart to cell volume.

    Parameters:
        lattice_cart (list): Cartesian lattice vectors.

    Returns:
        float: cell volume in Angstrom^3.

    """
    lattice_cart = np.asarray(lattice_cart)
    vol = np.abs(np.dot(np.cross(lattice_cart[0], lattice_cart[1]), lattice_cart[2]))
    return vol


def cart2abc(lattice_cart):
    """ Convert Cartesian lattice vectors to lattice parametres.

    Parameters:
        lattice_cart (list): Cartesian lattice vectors.

    Returns:
        list: lattice parameters [[a,b,c],[alpha,beta,gamma]].

    """
    vecs = lattice_cart
    lattice_abc = []
    a, b, c = (np.sqrt(sum([val**2 for val in vec])) for vec in vecs)
    lattice_abc.append([a, b, c])
    # np.cos(alpha) = b.c /|b * c|
    radians2deg = 180.0/np.pi
    cos_alpha = sum([val_b*val_c for (val_c, val_b) in zip(vecs[2], vecs[1])]) / (b*c)
    cos_beta = sum([val_c*val_a for (val_c, val_a) in zip(vecs[2], vecs[0])]) / (c*a)
    cos_gamma = sum([val_a*val_b for (val_a, val_b) in zip(vecs[0], vecs[1])]) / (a*b)
    alpha = radians2deg * np.arccos(cos_alpha)
    beta = radians2deg * np.arccos(cos_beta)
    gamma = radians2deg * np.arccos(cos_gamma)
    lattice_abc.append([alpha, beta, gamma])
    return lattice_abc


def frac2cart(lattice_cart, positions_frac):
    """ Convert positions_frac block into positions_abs.

    Parameters:
        lattice_cart (list): Cartesian lattice vectors.
        positions_frac (list): list of fractional position vectors.

    Returns:
        list: list of absolute position vectors.

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


def wrap_frac_coords(positions, remove=False):
    """ Wrap the given fractional coordinates back into the cell.

    Parameters:
        positions (list): list of fractional position vectors, or
            a single position.

    Keyword arguments:
        remove (bool): if True, removes points exterior to the cell.

    Returns:
        list: list of wrapped fractional position vectors.

    """
    from copy import deepcopy
    wrapped = deepcopy(positions)
    single = False
    if len(wrapped) == 3 and isinstance(wrapped[0], float):
        wrapped = [wrapped]
        single = True
    if remove:
        to_remove = len(wrapped)*[False]
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


def switch_coords(lattice, pos, norm=None):
    """ Act on coordinates with the relevant lattice
    vectors to switch from fractional to absolute coordinates.

    Parameters:

        lattice (np.ndarray(3, 3)): either lattice_cart or reciprocal lattice_cart
        pos (np.ndarray(3, :)): input positions to convert

    Keyword arguments:
        norm (float): divide final coordinates by normalisation factor, e.g.
            2*np.pi when lattice is recip and positions are cartesian.

    Returns:
        np.ndarray(3, :): converted positions

    """
    new_pos = np.zeros_like(pos)
    for ni, _ in enumerate(pos):
        for j in range(3):
            new_pos[ni] += lattice[j] * pos[ni][j]
    if norm is not None:
        new_pos /= norm
    return new_pos


def cart2frac(lattice_cart, positions_abs):
    """ Convert positions_abs block into positions_frac (and equivalent
    in reciprocal space).

    Parameters:
        lattice_cart (list): Cartesian lattice vectors.
        positions_abs (list): list of absolute position vectors.

    Returns:
        list: list of fractional position vectors with the same shape
            as the input list.

    """
    _positions_abs = np.asarray(positions_abs, dtype=np.float64)
    reshaped = False
    if len(np.shape(_positions_abs)) == 1:
        reshaped = True
        _positions_abs = _positions_abs.reshape((1, 3))
    recip_lat = np.asarray(real2recip(lattice_cart))
    recip_lat = recip_lat.T
    positions_frac = switch_coords(recip_lat, _positions_abs, norm=2*np.pi)
    if reshaped:
        positions_frac = positions_frac.reshape(-1)
    return positions_frac.tolist()


def real2recip(real_lat):
    """ Convert the real lattice in Cartesian basis to
    the reciprocal space lattice.

    Parameters:
        real_lat (list): Cartesian lattice vectors.

    Returns:
        list: Cartesian lattice vectors of reciprocal lattice.

    """
    real_lat = np.asarray(real_lat)
    recip_lat = np.zeros((3, 3))
    volume = np.dot(real_lat[0], np.cross(real_lat[1], real_lat[2]))
    recip_lat[0] = (2*np.pi)*np.cross(real_lat[1], real_lat[2])
    recip_lat[1] = (2*np.pi)*np.cross(real_lat[2], real_lat[0])
    recip_lat[2] = (2*np.pi)*np.cross(real_lat[0], real_lat[1])
    recip_lat /= volume
    return recip_lat.tolist()


def calc_mp_grid(lattice_cart, spacing):
    """ Return correct Monkhorst-Pack grid based on lattice
    vectors and desired spacing.

    Parameters:
        lattice_cart (list): Cartesian lattice vectors.
        spacing (float): desired maximum grid spacing.

    Returns:
        list: list of 3 integers defining the MP grid.

    """
    recip_lat = real2recip(lattice_cart)
    recip_len = np.zeros((3))
    recip_len = np.sqrt(np.sum(np.power(recip_lat, 2), axis=1))
    mp_grid = recip_len / (2 * np.pi * spacing)
    return [int(np.ceil(elem)) for elem in mp_grid]


def shift_to_include_gamma(mp_grid):
    """ Calculate the shift required to include $\\Gamma$.
    in the Monkhorst-Pack grid.

    Parameters:
        mp_grid (:obj:`list` of :obj:`int`): number of grid points
            in each reciprocal space direction.

    Returns:
        :obj:`list` of :obj:`float`: shift required to include $\\Gamma$.

    """
    shift = [0, 0, 0]
    for ind, val in enumerate(mp_grid):
        if val % 2 == 0:
            shift[ind] = 1.0/(val*2)
    return shift


def shift_to_exclude_gamma(mp_grid):
    """ Calculate the shift required to exclude $\\Gamma$.
    in the Monkhorst-Pack grid. Returns the "minimal shift", i.e. only
    one direction will be shifted.

    Parameters:
        mp_grid (:obj:`list` of :obj:`int`): number of grid points
            in each reciprocal space direction.

    Returns:
        :obj:`list` of :obj:`float`: shift required to exclude $\\Gamma$.

    """
    shift = [0, 0, 0]
    if all([val % 2 == 1 for val in mp_grid]):
        for ind, val in enumerate(mp_grid):
            if val % 2 == 1:
                shift[ind] = 1.0/(val*2)
                break

    return shift


def get_best_mp_offset_for_cell(doc):
    """ Calculates the "best" kpoint_mp_offset to use for the passed
    cell. If the crystal has a hexagonal space group, then the offset
    returned will shift the grid to include $\\Gamma$ point, and vice
    versa for non-hexagonal cells.

    Parameters:
        doc (dict): matador document to consider, containing structural
            information and a "kpoints_mp_spacing" key.

    Returns:
        :obj:`list` of :obj:`float`: the desired kpoint_mp_offset.

    """

    gamma = False
    if '6' in get_spacegroup_spg(doc, symprec=1e-5):
        gamma = True

    print(doc)

    if 'lattice_cart' not in doc or 'kpoints_mp_spacing' not in doc:
        raise RuntimeError('Unable to calculate offset without lattice or spacing')

    mp_grid = calc_mp_grid(doc['lattice_cart'], doc['kpoints_mp_spacing'])
    if gamma:
        return shift_to_include_gamma(mp_grid)

    return shift_to_exclude_gamma(mp_grid)


def calc_mp_spacing(real_lat, mp_grid, prec=3):
    """ Convert real lattice in Cartesian basis and the
    kpoint_mp_grid into a grid spacing.

    Parameters:
        real_lat (list): Cartesian lattice vectors.
        mp_grid (:obj:`list` of :obj:`int`): 3 integers defining the MP grid.

    Keyword arguments:
        prec (int): desired decimal precision of output.

    Returns:
        float: mp_spacing rounded to `prec`.

    """
    recip_lat = real2recip(real_lat)
    recip_len = np.zeros((3))
    recip_len = np.sqrt(np.sum(np.power(recip_lat, 2), axis=1))
    spacing = recip_len / (2*np.pi*np.asarray(mp_grid))
    max_spacing = np.max(spacing)
    exponent = round(np.log10(max_spacing) - prec)
    return round(max_spacing + 0.5*10**exponent, prec)


def get_seekpath_kpoint_path(
    doc, standardize=True, explicit=True, spacing=0.01, threshold=1e-7, debug=False, symmetry_tol=None
):
    """ Return the conventional kpoint path of the relevant crystal system
    according to the definitions by "HKPOT" in
    Comp. Mat. Sci. 128, 2017:

    http://dx.doi.org/10.1016/j.commatsci.2016.10.015

    Parameters:
        doc (dict/tuple): matador doc or spglib tuple to find kpoint path for.

    Keyword arguments:
        spacing (float): desired kpoint spacing
        threshold (float): internal seekpath threshold
        symmetry_tol (float): spglib symmetry tolerance

    Returns:
        dict: standardized version of input doc
        list: list of kpoint positions
        dict: full dictionary of all seekpath results

    """
    try:
        from seekpath import get_explicit_k_path, get_path
    except ImportError:
        raise ImportError("SeeK-Path dependency missing, please install it with `pip install seekpath`.")

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
        seekpath_results = get_explicit_k_path(spg_structure,
                                               reference_distance=spacing,
                                               with_time_reversal=True,
                                               symprec=symmetry_tol,
                                               threshold=threshold)

        kpt_path = seekpath_results['explicit_kpoints_rel']
    else:
        seekpath_results = get_path(spg_structure)
        kpt_path = []

    primitive_doc = dict()
    primitive_doc['lattice_cart'] = seekpath_results['primitive_lattice']
    primitive_doc['positions_frac'] = seekpath_results['primitive_positions']
    primitive_doc['atom_types'] = [str(elements[i]) for i in seekpath_results['primitive_types']]
    primitive_doc['num_atoms'] = len(primitive_doc['atom_types'])
    primitive_doc['lattice_abc'] = cart2abc(primitive_doc['lattice_cart'])
    primitive_doc['cell_volume'] = cart2volume(primitive_doc['lattice_cart'])
    if debug:
        print('Found lattice type {}'.format(seekpath_results['bravais_lattice_extended']))
        print('Old lattice:\n', np.asarray(doc['lattice_cart']))
        print('Contained {} atoms'.format(doc['num_atoms']))
        print('New lattice:\n', np.asarray(primitive_doc['lattice_cart']))
        print('Contains {} atoms'.format(primitive_doc['num_atoms']))
        print('k-point path contains {} points.'.format(len(kpt_path)))

    if 'site_occupancy' in doc:
        if min(doc['site_occupancy']) < 1 - EPS:
            print('Ignoring any site occupancy found in this cell.')
        primitive_doc['site_occupancy'] = [1 for atom in primitive_doc['atom_types']]

    return primitive_doc, kpt_path, seekpath_results


def doc2spg(doc, check_occ=True):
    """ Return an spglib input tuple from a matador doc.

    Parameters:
        doc (dict or :class:`Crystal`): matador document or Crystal object.

    Keyword arguments:
        check_occ (bool): check for partial occupancy and raise an error if
            present.

    Returns:
        tuple: spglib-style tuple of lattice, positions and types.

    """
    from matador.utils.chem_utils import get_atomic_number
    try:
        if 'lattice_cart' not in doc:
            doc['lattice_cart'] = abc2cart(doc['lattice_abc'])
        cell = (doc['lattice_cart'],
                doc['positions_frac'],
                [get_atomic_number(elem) for elem in doc['atom_types']])
        if check_occ and np.min(doc.get('site_occupancy', [1.0])) < 1.0:
            raise RuntimeError("spglib does not support partial occupancy.")

    except KeyError:
        raise RuntimeError('doc2spg failed, matador document was missing data!')

    return cell


def get_space_group_label_latex(label):
    """ Return the LaTeX format of the passed space group label. Takes
    any string, leaves the first character upright, italicses the rest,
    handles subscripts and bars over numbers.

    Parameters:
        label (str): a given space group in "standard" plain text format,
        e.g. P-63m.

    Returns:
        str: the best attempt to convert the label to LaTeX.

    """
    latex_label = '$'
    if not isinstance(label, str):
        raise RuntimeError("Space group label must be a string, not {}".format(label))

    skip = False
    for ind, char in enumerate(label):
        if skip:
            skip = False
            continue

        if char == '-':
            # add the next char inside the bar
            latex_label += "\\bar{" + label[ind+1] + "}"
            skip = True

        else:
            latex_label += char

    latex_label += "$"

    return latex_label


def standardize_doc_cell(doc, primitive=True, symprec=1e-2):
    """ Return standardized cell data from matador doc.

    Parameters:
        doc (dict or :class:`Crystal`): matador document or Crystal object.

    Keyword arguments:
        primitive (bool): whether to reduce cell to primitive.
        symprec (float): spglib symmetry tolerance.

    Returns:
        dict: matador document containing standardized cell.

    """
    import spglib as spg
    from matador.crystal import Crystal
    from matador.utils.chem_utils import get_atomic_symbol
    from copy import deepcopy

    spg_cell = doc2spg(doc)
    spg_standardized = spg.standardize_cell(spg_cell, to_primitive=primitive, symprec=symprec)
    if not isinstance(doc, Crystal):
        std_doc = deepcopy(doc)
    else:
        std_doc = deepcopy(doc._data)
    std_doc['lattice_cart'] = [list(vec) for vec in spg_standardized[0]]
    std_doc['lattice_abc'] = cart2abc(std_doc['lattice_cart'])
    std_doc['positions_frac'] = [list(atom) for atom in spg_standardized[1]]
    std_doc['atom_types'] = [get_atomic_symbol(atom) for atom in spg_standardized[2]]
    std_doc['site_occupancy'] = len(std_doc['positions_frac']) * [1]
    std_doc['cell_volume'] = cart2volume(std_doc['lattice_cart'])
    std_doc['space_group'] = get_spacegroup_spg(std_doc, symprec=symprec)
    # if the original document was a crystal, return a new one
    if isinstance(doc, Crystal):
        std_doc = Crystal(std_doc)

    return std_doc


def get_spacegroup_spg(doc, symprec=0.01, check_occ=True):
    """ Return spglib spacegroup for a cell.

    Parameters:
        doc (dict or :class:`Crystal`): matador document or Crystal object.

    Keyword arguments:
        symprec (float): spglib symmetry tolerance.

    Returns:
        str: spacegroup symbol of structure.

    """
    import spglib as spg
    spg_cell = doc2spg(doc, check_occ=check_occ)
    space_group = spg.get_spacegroup(spg_cell, symprec=symprec)
    if space_group is None:
        raise RuntimeError('Spglib was unable to calculate space group.')

    return space_group.split(' ')[0]


def add_noise(doc, amplitude=0.1):
    """ Add random noise to the positions of structure contained in doc.
    Useful for force convergence tests.

    Parameters:
        doc (dict): dictionary containing matador structure.

    Keyword arguments:
        amplitude (float): maximum amplitude of noise vector.

    Raises:
        KeyError if (`lattice_cart` and `positions_frac`) or `positions_abs`
            are missing.

    Returns:
        dict: the randomised structure.
    """
    poscart = np.asarray(doc.get('positions_abs') or frac2cart(doc['lattice_cart'], doc['positions_frac']))
    for atom in poscart:
        noise_vector = np.random.rand(3)
        noise_vector /= np.sqrt(np.sum(noise_vector**2))
        noise_vector *= amplitude*(np.random.rand())
        atom += noise_vector
    doc['positions_abs'] = poscart.tolist()
    doc['positions_frac'] = cart2frac(doc['lattice_cart'], doc['positions_abs'])

    return doc


def calc_pairwise_distances_pbc(poscart, images, lattice, rmax,
                                poscart_b=None, compress=False, debug=False,
                                filter_zero=False, per_image=False):
    """ Calculate PBC distances with SciPy's cdist, given the
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
        array_size = (len(_poscart) * len(images) * len(_poscart_b))
        distances = np.empty(array_size)
        for image_ind, prod in enumerate(images):
            distances[image_ind*num_pairs:(image_ind+1)*num_pairs] = cdist(_poscart, _poscart_b + prod @ _lattice).flatten()

        distances = np.ma.masked_where(distances > rmax, distances, copy=False)
        if filter_zero:
            distances = np.ma.masked_where(distances < EPS, distances, copy=False)

        if debug:
            print('Calculated: {}, Used: {}, Ignored: {}'.format(len(distances),
                                                                 np.ma.count(distances),
                                                                 np.ma.count_masked(distances)))
        if compress:
            distances = distances.compressed()

    if debug:
        end = time.time()
        print('Calculated distances in {} s'.format(end - start))

    return distances


def create_simple_supercell(doc, extension, standardize=False, symmetric=False):
    """ Return a document with new supercell, given extension vector.

    Parameters:
        doc (dict): matador doc to construct cell from.
        extension (:obj:`tuple` of :obj:`int`): multiplicity of each lattice vector,
            e.g. (2,2,1).

    Keyword arguments:
        standardize (bool): whether or not to use spglib to standardize the cell first.
        symmetric (bool): whether or not centre the new cell on the origin.

    Returns:
        supercell_doc (dict): matador document containing supercell.

    """
    from itertools import product
    from copy import deepcopy
    if not all([elem >= 1 for elem in extension]) or extension == (1, 1, 1):
        raise RuntimeError('Bad/null supercell {} requested...'.format(extension))

    supercell_doc = deepcopy(doc)
    # standardize cell with spglib
    if standardize:
        supercell_doc = standardize_doc_cell(supercell_doc)

    # copy new doc and delete data that will not be corrected
    keys_to_copy = set(['positions_frac', 'lattice_cart', 'atom_types', 'source', 'stoichiometry'])
    supercell_doc = {key: supercell_doc[key] for key in keys_to_copy}

    for i, elem in enumerate(extension):
        for k in range(3):
            supercell_doc['lattice_cart'][i][k] = elem * supercell_doc['lattice_cart'][i][k]
        for ind, atom in enumerate(supercell_doc['positions_frac']):
            supercell_doc['positions_frac'][ind][i] /= elem

    images = product(*[list(range(elem)) for elem in extension])
    new_positions = []
    new_atoms = []
    for image in images:
        if image != (0, 0, 0):
            for ind, atom in enumerate(supercell_doc['atom_types']):
                new_pos = deepcopy(supercell_doc['positions_frac'][ind])
                for i, elem in enumerate(image):
                    new_pos[i] += image[i] / extension[i]
                new_positions.append(new_pos)
                new_atoms.append(atom)
    supercell_doc['atom_types'].extend(new_atoms)
    supercell_doc['positions_frac'].extend(new_positions)

    if symmetric:
        supercell_doc['positions_frac'] = np.asarray(supercell_doc['positions_frac'])
        supercell_doc['positions_frac'] -= np.mean(supercell_doc['positions_frac'], axis=0)
        supercell_doc['positions_frac'] = supercell_doc['positions_frac'].tolist()

    supercell_doc['num_atoms'] = len(supercell_doc['atom_types'])
    supercell_doc['cell_volume'] = cart2volume(supercell_doc['lattice_cart'])
    supercell_doc['lattice_abc'] = cart2abc(supercell_doc['lattice_cart'])

    assert np.isclose(supercell_doc['cell_volume'],
                      np.prod(extension)*cart2volume(doc['lattice_cart']))

    return supercell_doc
