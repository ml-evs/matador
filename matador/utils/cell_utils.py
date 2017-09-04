""" This file implements some useful functions
for real/reciprocal cell manipulation and sampling.
"""

# external libraries
import numpy as np
# standard library
from math import pi, cos, sin, sqrt, acos, log10, ceil
from functools import reduce


def abc2cart(lattice_abc):
    """ Convert lattice_abc=[[a,b,c],[alpha,beta,gamma]]
    (in degrees) to lattice vectors
    lattice_cart=[[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]].
    """
    assert len(lattice_abc) == 2
    assert len(lattice_abc[0]) == 3
    assert len(lattice_abc[1]) == 3
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
    spacing = recip_len / (2*pi*np.asarray(mp_grid))
    max_spacing = np.max(spacing)
    exponent = round(log10(max_spacing) - prec)
    return round(max_spacing + 0.5*10**exponent, prec)


def get_crystal_system(lattice_cart, to_primitive=True, eps=5e-2, symprec=1e-2, debug=False):
    """ Return the name of the crystal system according to the definitions
    by Setyawana & Curtarolo in Comp. Mat. Sci. 49(2), 2010:

    http://dx.doi.org/10.1016/j.commatsci.2010.05.010.

    Input:

        | lattice_cart: list(list(float)), lattice vectors in format
                        [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]].

    Args:

        | to_primitive: bool, whether to use the spglib primitive cell to find
                        the best BZ path. If using this function to plot a bandstructure,
                        this should be False (DEFAULT).
        | eps         : float, numpy tolerance.
        | symprec     : float, spglib tolerance.

    """
    import spglib
    positions = [[0, 0, 0]]
    numbers = [1]
    spg_cell = (lattice_cart, positions, numbers)
    refined_cell = spglib.standardize_cell(spg_cell, to_primitive=to_primitive, symprec=symprec)
    refined_cart = refined_cell[0]
    assert refined_cell is not None, 'Unable to standardize cell'
    refined_abc = cart2abc(refined_cell[0])
    a = refined_abc[0][0]
    b = refined_abc[0][1]
    c = refined_abc[0][2]
    alpha = refined_abc[1][0]
    beta = refined_abc[1][1]
    gamma = refined_abc[1][2]
    if debug:
        print('Standardized cell:')
        print('{l[0]:5.2f} A, {l[1]:5.2f} A, {l[2]:5.2f} A'.format(l=refined_abc[0]))
        print('{l[0]:6.2f}°, {l[1]:6.2f}°, {l[2]:6.2f}°'.format(l=refined_abc[1]))

    def close(values):
        return np.allclose([values[-1]], values, atol=eps, rtol=eps)

    if close([alpha, beta, gamma, 90]):
        if close([a, b, c]):
            if to_primitive:
                if close([refined_cart[i][i] for i in range(3)]+[0]):
                    return 'fcc'
                elif close([refined_cart[i][i] for i in range(3)] + [-refined_cart[0][1]]):
                    return 'bcc'
                else:
                    return 'cubic'
            else:
                return 'cubic'
        elif close([a, b]) or close([b, c]) or close([a, c]):
            return 'tetragonal'
        else:
            return 'orthorhombic'
    elif close([alpha, gamma, 90]):
        return 'monoclinic'
    elif close([alpha, beta, 90]) and close([gamma, 120]) and close([a, b]):
        return 'hexagonal'
    elif close([alpha, beta, gamma]) and not close([alpha, 90]) and close([a, b, c]):
        if alpha > 90:
            return 'rhombohedral type 2'
        else:
            return 'rhombohedral type 1'
    else:
        return 'triclinic'


def get_bs_kpoint_path(lattice_cart, spacing=0.01, debug=False):
    """ Return the conventional kpoint path of the relevant crystal system
    according to the definitions by Setyawana & Curtarolo in
    Comp. Mat. Sci. 49(2), 2010:

    http://dx.doi.org/10.1016/j.commatsci.2010.05.010.

    Input:

        | lattice_cart: list(list(float)), lattice vectors in format
                        [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]].

    Returns:

        | path            : list(list(float)), positions of kpoints.
        | critical_points : list(str), labels of critical points.

    """
    system = get_crystal_system(lattice_cart, debug=debug)
    assert system != 'rhombedral type 2', 'Special points missing for rhombohedral type 2'
    critical_points = SPECIAL_KPOINT_PATHS[system].split(',')
    if debug:
        print('Crystal system detected as:', system)
        print('Constructing path:')
        print('-->'.join(critical_points))
    path = []
    special_points = get_special_kpoints_for_lattice(system, lattice_cart)
    for i in range(len(critical_points)-1):
        point_label = critical_points[i]
        if point_label is '|':
            continue
        else:
            point = np.asarray(special_points[point_label])

        next_point_label = critical_points[i+1]
        if next_point_label is '|':
            path.append(point.tolist())
            continue
        else:
            next_point = np.asarray(special_points[next_point_label])

        diff = next_point - point
        diff_mag = np.sqrt(np.sum(diff**2))
        num_kpts = ceil(diff_mag / spacing)
        diff_vec = diff / num_kpts
        print(diff_vec, diff, num_kpts)
        print(point_label, next_point_label)
        print(point, next_point)
        for j in range(num_kpts):
            path.append((point + j*diff_vec).tolist())
        if i == len(critical_points) - 2:
            path.append(next_point.tolist())

    for point in critical_points:
        if point != '|':
            assert special_points[point] in path, '{} is missing'.format(point)

    return critical_points, path


def get_special_kpoints_for_lattice(crystal_system, lattice_cart):
    """ High-symmetry points in the IBZ for the different crystal systems.

    Taken from ASE,

    https://wiki.fysik.dtu.dk/ase/_modules/ase/dft/kpoints.html

    who themselves took it from Setyawana & Curtarolo, Comp. Mat. Sci. 49(2), 2010:

    http://dx.doi.org/10.1016/j.commatsci.2010.05.010.

    Input:

        | crystal_system : str, one of e.g. 'hexagonal', 'monoclinic'.
        | lattice_cart   : list(list(float), Cartesian lattice vectors.

    """
    if crystal_system is 'monoclinic':
        lattice_abc = cart2abc(lattice_cart)
        b = lattice_abc[0][1]
        c = lattice_abc[0][2]
        alpha = lattice_abc[1][0]
        eta = (1 - b * cos((pi/180)*alpha) / c) / (2 * sin((pi/180)*alpha)**2)
        nu = 0.5 - eta * c * cos((pi/180)*alpha) / b
    elif crystal_system is 'rhombohedral type 1':
        lattice_abc = cart2abc(lattice_cart)
        alpha = lattice_abc[1][0]
        eta = (1 + 4 * cos((pi/180)*alpha)) / (2 + 4 * cos((pi/180)*alpha))
        nu = 0.75 - eta / 2
        print('nu: ', nu, 'eta:', eta)
    else:
        nu = 0
        eta = 0
    return SPECIAL_KPOINTS(crystal_system, nu, eta)


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


def SPECIAL_KPOINTS(crystal_system, nu, eta):
    points = {
        'cubic': {'G': [0, 0, 0],
                  'M': [1 / 2, 1 / 2, 0],
                  'R': [1 / 2, 1 / 2, 1 / 2],
                  'X': [0, 1 / 2, 0]},
        'fcc': {'G': [0, 0, 0],
                'K': [3 / 8, 3 / 8, 3 / 4],
                'L': [1 / 2, 1 / 2, 1 / 2],
                'U': [5 / 8, 1 / 4, 5 / 8],
                'W': [1 / 2, 1 / 4, 3 / 4],
                'X': [1 / 2, 0, 1 / 2]},
        'bcc': {'G': [0, 0, 0],
                'H': [1 / 2, -1 / 2, 1 / 2],
                'P': [1 / 4, 1 / 4, 1 / 4],
                'N': [0, 0, 1 / 2]},
        'tetragonal': {'G': [0, 0, 0],
                       'A': [1 / 2, 1 / 2, 1 / 2],
                       'M': [1 / 2, 1 / 2, 0],
                       'R': [0, 1 / 2, 1 / 2],
                       'X': [0, 1 / 2, 0],
                       'Z': [0, 0, 1 / 2]},
        'orthorhombic': {'G': [0, 0, 0],
                         'R': [1 / 2, 1 / 2, 1 / 2],
                         'S': [1 / 2, 1 / 2, 0],
                         'T': [0, 1 / 2, 1 / 2],
                         'U': [1 / 2, 0, 1 / 2],
                         'X': [1 / 2, 0, 0],
                         'Y': [0, 1 / 2, 0],
                         'Z': [0, 0, 1 / 2]},
        'hexagonal': {'G': [0, 0, 0],
                      'A': [0, 0, 1 / 2],
                      'H': [1 / 3, 1 / 3, 1 / 2],
                      'K': [1 / 3, 1 / 3, 0],
                      'L': [1 / 2, 0, 1 / 2],
                      'M': [1 / 2, 0, 0]},
        'monoclinic': {'G': [0, 0, 0],
                       'A': [1 / 2, 1 / 2, 0],
                       'C': [0, 1 / 2, 1 / 2],
                       'D': [1 / 2, 0, 1 / 2],
                       'D1': [1 / 2, 0, -1 / 2],
                       'E': [1 / 2, 1 / 2, 1 / 2],
                       'H': [0, eta, 1 - nu],
                       'H1': [0, 1 - eta, nu],
                       'H2': [0, eta, -nu],
                       'M': [1 / 2, eta, 1 - nu],
                       'M1': [1 / 2, 1 - eta, nu],
                       'M2': [1 / 2, eta, -nu],
                       'X': [0, 1 / 2, 0],
                       'Y': [0, 0, 1 / 2],
                       'Y1': [0, 0, -1 / 2],
                       'Z': [1 / 2, 0, 0]},
        'rhombohedral type 1': {'G': [0, 0, 0],
                                'B': [eta, 1 / 2, 1 - eta],
                                'B1': [1 / 2, 1 - eta, eta - 1],
                                'F': [1 / 2, 1 / 2, 0],
                                'L': [1 / 2, 0, 0],
                                'L1': [0, 0, - 1 / 2],
                                'P': [eta, nu, nu],
                                'P1': [1 - nu, 1 - nu, 1 - eta],
                                'P2': [nu, nu, eta - 1],
                                'Q': [1 - nu, nu, 0],
                                'X': [nu, 0, -nu],
                                'Z': [0.5, 0.5, 0.5]}}
    return points[crystal_system]


SPECIAL_KPOINT_PATHS = {
    'cubic': 'G,X,M,G,R,X,|,M,R',
    'fcc': 'G,X,W,K,G,L,U,W,L,K,|,U,X',
    'bcc': 'G,H,N,G,P,H,|,P,N',
    'tetragonal': 'G,X,M,G,Z,R,A,Z,X,R,|,M,A',
    'orthorhombic': 'G,X,S,Y,G,Z,U,R,T,Z,|,Y,T,|,U,X,|,S,R',
    'hexagonal': 'G,M,K,G,A,L,H,A,|,L,M,|,K,H',
    'monoclinic': 'G,Y,H,C,E,M1,A,X,H1,|,M,D,Z,|,Y,D',
    'rhombohedral type 1': 'G,L,B1,|,B,Z,G,X,|,Q,F,P1,Z,|,L,P',
    'rhombohedral type 2': 'G,P,Z,Q,G,F,P1,Q1,L,Z'}
