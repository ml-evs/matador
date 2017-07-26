""" This file implements some useful functions
for cell manipulation.
"""

# external libraries
import numpy as np
# standard library
from math import pi, cos, sin, sqrt, acos, log10
from functools import reduce


def abc2cart(lattice_abc):
    """ Convert lattice_abc=[[a,b,c],[alpha,beta,gamma]]
    (in degrees) to lattice vectors
    lattice_cart=[[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]].
    """
    a = lattice_abc[0][0]
    b = lattice_abc[0][1]
    c = lattice_abc[0][2]
    deg2rad = pi/180
    alpha = lattice_abc[1][0] * deg2rad
    beta = lattice_abc[1][1] * deg2rad
    gamma = lattice_abc[1][2] * deg2rad
    lattice_cart = []
    lattice_cart.append([a, 0.0, 0.0])
    # vec(b) = (b cos(gamma), b sin(gamma), 0)
    bx = b*cos(gamma)
    by = b*sin(gamma)
    tol = 1e-12
    if abs(bx) < tol:
        bx = 0.0
    if abs(by) < tol:
        by = 0.0
    cx = c*cos(beta)
    if abs(cx) < tol:
        cx = 0.0
    cy = c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
    if abs(cy) < tol:
        cy = 0.0
    cz = sqrt(c**2 - cx**2 - cy**2)
    lattice_cart.append([bx, by, 0.0])
    lattice_cart.append([cx, cy, cz])
    return lattice_cart


def cart2abcstar(lattice_cart):
    """ Convert lattice_cart =[[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]]
    to a*, b*, c*.
    """
    lattice_cart = np.asarray(lattice_cart)
    lattice_star = np.zeros_like(lattice_cart)
    lattice_star[0] = np.cross(lattice_cart[1], lattice_cart[2])
    lattice_star[1] = np.cross(lattice_cart[2], lattice_cart[0])
    lattice_star[2] = np.cross(lattice_cart[0], lattice_cart[1])
    vol = np.linalg.det(lattice_star)
    vol = np.dot(np.cross(lattice_cart[0], lattice_cart[1]), lattice_cart[2])
    lattice_star /= vol
    return lattice_star.tolist()


def cart2volume(lattice_cart):
    """ Convert lattice_cart to cell volume. """
    lattice_cart = np.asarray(lattice_cart)
    vol = np.dot(np.cross(lattice_cart[0], lattice_cart[1]), lattice_cart[2])
    assert vol > 0
    return vol


def cart2abc(lattice_cart):
    """ Convert lattice_cart =[[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]]
    to lattice_abc=[[a,b,c],[alpha,beta,gamma]].
    """
    vec_a = lattice_cart[0]
    vec_b = lattice_cart[1]
    vec_c = lattice_cart[2]
    lattice_abc = []
    a = 0
    b = 0
    c = 0
    for i in range(3):
        a += vec_a[i]**2
        b += vec_b[i]**2
        c += vec_c[i]**2
    a = sqrt(a)
    b = sqrt(b)
    c = sqrt(c)
    lattice_abc.append([a, b, c])
    # cos(alpha) = b.c /|b * c|
    cos_alpha = 0
    cos_beta = 0
    cos_gamma = 0
    for i in range(3):
        cos_alpha += vec_b[i] * vec_c[i]
        cos_beta += vec_c[i] * vec_a[i]
        cos_gamma += vec_a[i] * vec_b[i]
    cos_alpha /= b*c
    cos_beta /= c*a
    cos_gamma /= a*b
    alpha = 180.0*acos(cos_alpha)/pi
    beta = 180.0*acos(cos_beta)/pi
    gamma = 180.0*acos(cos_gamma)/pi
    lattice_abc.append([alpha, beta, gamma])
    return lattice_abc


def frac2cart(lattice_cart, positions_frac):
    """ Convert positions_frac block into positions_abs. """
    positions_frac = np.asarray(positions_frac)
    lattice_cart = np.asarray(lattice_cart)
    positions_abs = np.zeros_like(positions_frac)
    assert(len(lattice_cart) == 3)
    for i in range(len(positions_frac)):
        for j in range(3):
            positions_abs[i] += lattice_cart[j]*positions_frac[i][j]
    return positions_abs.tolist()


def real2recip(real_lat):
    """ Convert the real lattice in Cartesian basis to
    the reciprocal space lattice.
    """
    real_lat = np.asarray(real_lat)
    recip_lat = np.zeros((3, 3))
    recip_lat[0] = (2*pi)*np.cross(real_lat[1], real_lat[2]) / \
        (np.dot(real_lat[0], np.cross(real_lat[1], real_lat[2])))
    recip_lat[1] = (2*pi)*np.cross(real_lat[2], real_lat[0]) / \
        (np.dot(real_lat[1], np.cross(real_lat[2], real_lat[0])))
    recip_lat[2] = (2*pi)*np.cross(real_lat[0], real_lat[1]) / \
        (np.dot(real_lat[2], np.cross(real_lat[0], real_lat[1])))
    return recip_lat.tolist()


def calc_mp_spacing(real_lat, mp_grid, prec=2):
    """ Convert real lattice in Cartesian basis and the
    kpoint_mp_grid into a grid spacing.
    """
    recip_lat = real2recip(real_lat)
    recip_len = np.zeros((3))
    recip_len = np.sqrt(np.sum(np.power(recip_lat, 2), axis=1))
    max_spacing = 0
    for j in range(3):
        spacing = recip_len[j] / (2*pi*mp_grid[j])
        max_spacing = (spacing if spacing > max_spacing else max_spacing)
    exponent = round(log10(max_spacing) - 1)
    return round(max_spacing + 0.5*10**exponent, prec)


def doc2spg(doc):
    """ Return an spglib input tuple from a matador doc. """
    from .chem_utils import get_atomic_number
    try:
        if 'lattice_cart' not in doc:
            doc['lattice_cart'] = abc2cart(doc['lattice_abc'])
        cell = (doc['lattice_cart'],
                doc['positions_frac'],
                [get_atomic_number(elem) for elem in doc['atom_types']])
        return cell
    except:
        return False


def standardize_doc_cell(doc):
    """ Inserts spglib e.g. standardized cell data
    into matador doc. """
    import spglib as spg
    from .chem_utils import get_atomic_symbol
    spg_cell = doc2spg(doc)
    spg_standardized = spg.standardize_cell(spg_cell)
    doc['lattice_cart'] = [list(vec) for vec in spg_standardized[0]]
    doc['positions_frac'] = [list(atom) for atom in spg_standardized[1]]
    doc['atom_types'] = [get_atomic_symbol(atom) for atom in spg_standardized[2]]
    return doc


def get_spacegroup_spg(doc, symprec=0.01):
    """ Return spglib spacegroup for a cell. """
    from spglib import get_spacegroup
    spg_cell = doc2spg(doc)
    return get_spacegroup(spg_cell, symprec=symprec).split(' ')[0]


def create_simple_supercell(seed_doc, extension, standardize=False):
    """ Return a document with new supercell, given extension vector.

    Input:

        doc        : dict, matador doc to construct cell from.
        extension  : tuple(int), multiplicity of each lattice vector, e.g. (2,2,1).
        standardize: bool, whether or not to use spglib to standardize the cell first.

    """
    from itertools import product
    from copy import deepcopy
    assert(all([elem >= 1 for elem in extension]))
    assert(all([int(elem) for elem in extension]))
    if extension == (1, 1, 1):
        raise RuntimeError('Redundant supercell {} requested...'.format(extension))

    num_images = reduce(lambda x, y: x*y, extension)
    doc = deepcopy(seed_doc)
    # standardize cell with spglib
    if standardize:
        doc = standardize_doc_cell(doc)

    # copy new doc and delete data that will not be corrected
    supercell_doc = deepcopy(doc)
    if 'positions_abs' in supercell_doc:
        del supercell_doc['positions_abs']
    if 'lattice_abc' in supercell_doc:
        del supercell_doc['lattice_abc']

    for i, elem in enumerate(extension):
        for k in range(3):
            supercell_doc['lattice_cart'][i][k] = elem * doc['lattice_cart'][i][k]
        for ind, atom in enumerate(supercell_doc['positions_frac']):
            supercell_doc['positions_frac'][ind][i] /= elem

    images = product(*[list(range(elem)) for elem in extension])
    _iter = 0
    new_positions = []
    new_atoms = []
    for image in images:
        _iter += 1
        if image == (0, 0, 0):
            continue
        else:
            for ind, atom in enumerate(supercell_doc['atom_types']):
                new_pos = deepcopy(supercell_doc['positions_frac'][ind])
                for i, elem in enumerate(image):
                    new_pos[i] += image[i] / extension[i]
                new_positions.append(new_pos)
                new_atoms.append(atom)
    supercell_doc['atom_types'].extend(new_atoms)
    supercell_doc['positions_frac'].extend(new_positions)
    supercell_doc['num_atoms'] = len(supercell_doc['atom_types'])
    supercell_doc['cell_volume'] = cart2volume(supercell_doc['lattice_cart'])
    supercell_doc['lattice_abc'] = cart2abc(supercell_doc['lattice_cart'])
    assert np.isclose(supercell_doc['cell_volume'], num_images*doc['cell_volume'])
    assert _iter == num_images
    return supercell_doc
