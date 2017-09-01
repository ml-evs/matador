# coding: utf-8
""" This file defines some useful chemistry. """

# external libraries
import numpy as np

# global consts
FARADAY_CONSTANT_Cpermol = 96.485332e3
Cperg_to_mAhperg = 2.778e-1
HARTREE_TO_EV = 27.21139


def get_periodic_table():
    """ Return some periodic table macros. """
    periodic_table = dict()
    periodic_table['I'] = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
    periodic_table['II'] = ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']
    periodic_table['III'] = ['B', 'Al', 'Ga', 'In', 'Tl']
    periodic_table['IV'] = ['C', 'Si', 'Ge', 'Sn', 'Pb']
    periodic_table['V'] = ['N', 'P', 'As', 'Sb', 'Bi']
    periodic_table['VI'] = ['O', 'S', 'Se', 'Te', 'Po']
    periodic_table['VII'] = ['F', 'Cl', 'Br', 'I', 'At']
    periodic_table['Tran'] = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                              'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                              'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
    periodic_table['Lan'] = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                             'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    periodic_table['Act'] = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
                             'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
    periodic_table['X'] = [elem for group in periodic_table.keys() for elem in periodic_table[group]]
    return periodic_table


def get_molar_mass(elem):
    """ Returns molar mass of chosen element. """
    import periodictable
    return periodictable.elements.symbol(elem).mass


def get_atomic_number(elem):
    """ Returns atomic number of chosen element. """
    import periodictable
    return periodictable.elements.symbol(elem).number


def get_atomic_symbol(atomic_number):
    """ Returns elemental symbol from atomic number. """
    import periodictable
    return periodictable.elements[atomic_number].symbol


def get_concentration(doc, elements):
    """ Returns x for A_x B_{1-x}
    or x,y for A_x B_y C_z, (x+y+z=1). """
    concs = [0.0] * (len(elements)-1)
    for ind, elem in enumerate(doc['stoichiometry']):
        if elem[0] in elements[:-1]:
            concs[elements.index(elem[0])] = elem[1]/float(get_atoms_per_fu(doc))
    return concs


def get_num_intercalated(cursor):
    """ Return array of the number of intercalated atoms
    per host atom from a list of structures. """
    from .cursor_utils import get_array_from_cursor
    x = np.zeros((len(cursor)))
    comps = get_array_from_cursor(cursor, 'concentration')
    for idx, comp in enumerate(comps):
        if len(comp) > 1:
            comp = comp[0]
        if 1-comp == 0:
            x[idx] = np.NaN
        else:
            x[idx] = comp/(1-comp)
    return x


def get_binary_grav_capacities(x, m_B):
    """ Returns capacity in mAh/g from x/y in A_x B_y
    and m_B in a.m.u.
    """
    x = np.array(x)
    if m_B != 0:
        return x * FARADAY_CONSTANT_Cpermol * Cperg_to_mAhperg / m_B
    else:
        return float('NaN')


def get_generic_grav_capacity(concs, elements):
    """ Returns gravimetric capacity of
    <elements[0]> in mAh/g of matador doc.
    """
    tmp_concs = np.array(concs, copy=True)
    # if no Li, capacity = 0...
    # tmp_concs /= np.min(concs)
    x = tmp_concs[0]
    if x == 0:
        return 0.0
    masses = dict()
    m_B = 0
    for elem in elements:
        masses[elem] = get_molar_mass(elem)
    for ind, elem in enumerate(elements):
        if ind == 0:
            continue
        else:
            m_B += masses[elem]*tmp_concs[ind]
    Q = get_binary_grav_capacities(x, m_B)
    return Q


def get_atoms_per_fu(doc):
    """ Calculate the number of atoms per formula unit. """
    atoms_per_fu = 0
    for j in range(len(doc['stoichiometry'])):
        atoms_per_fu += doc['stoichiometry'][j][1]
    return atoms_per_fu


def get_formation_energy(chempots, doc):
    """ From given chemical potentials, calculate the simplest
    formation energy per atom of the desired document.
    """
    formation = doc['enthalpy_per_atom']
    num_atoms_per_fu = get_atoms_per_fu(doc)
    for mu in chempots:
        for j in range(len(doc['stoichiometry'])):
            for i in range(len(mu['stoichiometry'])):
                if mu['stoichiometry'][i][0] == doc['stoichiometry'][j][0]:
                    formation -= (mu['enthalpy_per_atom'] * doc['stoichiometry'][j][1] /
                                  num_atoms_per_fu)
    return formation


def get_stoich(atom_types):
    """ Return integer stoichiometry from atom_types list.

    Input:

        atom_types: list of elements, e.g. ['Li', 'P', 'P']

    Returns:

        stoich : [['Li', 1], ['P', 2]]

    """
    from collections import defaultdict
    from math import gcd
    stoich = defaultdict(float)
    for atom in atom_types:
        if atom not in stoich:
            stoich[atom] = 0
        stoich[atom] += 1
    gcd_val = 0
    for atom in atom_types:
        if gcd_val == 0:
            gcd_val = stoich[atom]
        else:
            gcd_val = gcd(stoich[atom], gcd_val)
    # convert stoichiometry to tuple for fryan
    temp_stoich = []
    try:
        for key, value in stoich.items():
            if float(value)/gcd_val % 1 != 0:
                temp_stoich.append([key, float(value)/gcd_val])
            else:
                temp_stoich.append([key, value/gcd_val])
    except AttributeError:
        for key, value in stoich.iteritems():
            if float(value)/gcd_val % 1 != 0:
                temp_stoich.append([key, float(value)/gcd_val])
            else:
                temp_stoich.append([key, value/gcd_val])
    return sorted(temp_stoich)


def get_formula_from_stoich(stoich, elements=None, tex=False):
    """ Get the chemical formula of a structure from
    its matador stoichiometry.

    Input: stoich = [['Li', 3.0], ['P', 2.0]]
    Output: 'Li3P2'

    """
    form = ''
    if not isinstance(stoich, list):
        stoich = stoich.tolist()
    if elements is not None:
        for targ_elem in elements:
            for elem in stoich:
                if elem[0] == targ_elem:
                    if elem[1] == 1:
                        form += elem[0]
                    elif int(elem[1]) != 0:
                        if tex:
                            form += elem[0] + '\\textsubscript{' + str(int(elem[1])) + '}'
                        else:
                            form += elem[0] + str(int(elem[1]))
        assert form != ''
    else:
        for elem in sorted(stoich):
            if elem[1] == 1:
                form += elem[0]
            elif int(elem[1]) != 0:
                if tex:
                    form += elem[0] + '\\textsubscript{' + str(int(elem[1])) + '}'
                else:
                    form += elem[0] + str(int(elem[1]))
    return form
