# coding: utf-8
""" This file defines some useful chemistry. """

# external libraries
import numpy as np

# global consts
FARADAY_CONSTANT_Cpermol = 96.485332e3
Cperg_to_mAhperg = 2.778e-1
C_TO_mAh = Cperg_to_mAhperg
HARTREE_TO_EV = 27.21139
BOHR_TO_ANGSTROM = 0.529177211
RY_TO_EV = 13.605693009
KBAR_TO_GPA = 0.1
AVOGADROS_NUMBER = 6.022141e23
ANGSTROM_CUBED_TO_CENTIMETRE_CUBED = 1e-24
ELECTRON_CHARGE = 1.6021766e-19


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
    periodic_table[
        'X'] = [elem for group in periodic_table.keys() for elem in periodic_table[group]]
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
    concs = [0.0] * (len(elements) - 1)
    for ind, elem in enumerate(doc['stoichiometry']):
        if elem[0] in elements[:-1]:
            concs[elements.index(elem[0])] = elem[1] / float(get_atoms_per_fu(doc))
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
        if 1 - comp == 0:
            x[idx] = np.NaN
        else:
            x[idx] = comp / (1 - comp)
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
            m_B += masses[elem] * tmp_concs[ind]
    Q = get_binary_grav_capacities(x, m_B)
    return Q


def get_binary_volumetric_capacity(initial_doc, final_doc):
    """ For initial (delithiated/sodiated) (single element) structure
    and final (maximally charged) binary structure, calculate the volumetric capacity.

    Input:

        | initial_doc: dict, matador doc of delithiated phase
        | final_doc  : dict, matador doc of maximally lithiated phase

    Returns:

       | volumetric_capacity: float, capacity in mAh/cm^3.

    """

    assert len(initial_doc['stoichiometry']) == 1
    assert len(final_doc['stoichiometry']) == 2

    for _ion in final_doc['stoichiometry']:
        if _ion[0] != initial_doc['stoichiometry'][0][0]:
            ion = _ion[0]

    for species in final_doc['stoichiometry']:
        if species[0] == ion:
            num_ion = species[1]
        else:
            num_B = species[1]

    num_ions_per_initial_fu = num_ion / num_B
    volume_per_fu_cm3 = initial_doc[
        'cell_volume'] * ANGSTROM_CUBED_TO_CENTIMETRE_CUBED / initial_doc['num_fu']
    return ((num_ions_per_initial_fu / volume_per_fu_cm3) * (ELECTRON_CHARGE * Cperg_to_mAhperg))


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
                    formation -= (
                        mu['enthalpy_per_atom'] * doc['stoichiometry'][j][1] / num_atoms_per_fu)
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
            if float(value) / gcd_val % 1 != 0:
                temp_stoich.append([key, float(value) / gcd_val])
            else:
                temp_stoich.append([key, value / gcd_val])
    except AttributeError:
        for key, value in stoich.iteritems():
            if float(value) / gcd_val % 1 != 0:
                temp_stoich.append([key, float(value) / gcd_val])
            else:
                temp_stoich.append([key, value / gcd_val])
    return sorted(temp_stoich)


def get_ratios_from_stoichiometry(stoichiometry):
    """ Get a dictionary of pairwise atomic ratios.

    Parameters:
        stoichiometry (list): matador-style stoichiometry.

    Returns:
        dict: dictionary of pairwise ratios, e.g. for K8SnP4,
            ratio_dict = {'KSn': 8, 'KP': 2, 'SnP': 0.25,
                          'SnK': 0.125, 'PSn': 4, 'PK': 0.5}.

    """
    ratio_dict = dict()
    for i, elem_i in enumerate(stoichiometry):
        for j, elem_j in enumerate(stoichiometry):
            if elem_j != elem_i:
                ratio_dict[stoichiometry[i][0] + stoichiometry[j]
                           [0]] = round(float(stoichiometry[i][1]) / stoichiometry[j][1], 3)
    return ratio_dict


def get_stoich_from_formula(formula: str):
    """ Convert formula string, e.g. Li2TiP4 into a matador-style
    stoichiometry, e.g. [['Li', 2], ['Ti', 1], ['P', 4]].

    Input:

        | formula: str, chemical formula of compound

    Returns:

        | stoich: list, list of lists containing symbol and amount

    """
    from math import gcd
    import re

    parsed_elements = parse_element_string(formula, stoich=True)
    elements = []
    fraction = []
    for i, _ in enumerate(parsed_elements):
        if not bool(re.search(r'\d', parsed_elements[i])):
            elements.append(parsed_elements[i])
            try:
                fraction.append(float(parsed_elements[i + 1]))
            except:
                fraction.append(1.0)
    gcd_val = 0
    for frac in fraction:
        if gcd_val == 0:
            gcd_val = frac
        else:
            gcd_val = gcd(int(frac), int(gcd_val))
    fraction = np.asarray(fraction)
    fraction /= gcd_val
    stoich = [[elements[ind], fraction[ind]] for ind, _ in enumerate(elements)]
    return stoich


def parse_element_string(elements_str, stoich=False):
    """ Parse element query string with macros.

    e.g.
        Input: '[VII][Fe,Ru,Os][I]'
        Returns: ['[VII]', '[Fe,Ru,Os]', '[I]']

    e.g.2
        Input: '[VII]2[Fe,Ru,Os][I]'
        Returns: ['[VII]2', '[Fe,Ru,Os]', '[I]']

    Input:

        | elements_str: str, chemical formula, including macros.

    Args:

        | stoich: bool, parse as a stoichiometry, i.e. check for numbers

    Returns:

        | elements: list, split list of elements contained in input

    """
    import re
    valid = False
    for char in elements_str:
        if char.isupper():
            valid = True
    if not valid:
        exit('Composition must contain at least one upper case character.')
    elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements_str) if elem]
    if stoich:
        tmp_stoich = elements
        for ind, strng in enumerate(elements):
            if not any(char.isdigit() for char in strng):
                tmp_stoich[ind] = [strng]
            else:
                tmp_stoich[ind] = [elem for elem in re.split(r'([0-9]*)', strng) if elem]
        elements = [item for sublist in tmp_stoich for item in sublist]
    # split macros
    while '[' in elements or '][' in elements:
        tmp_stoich = list(elements)
        cleaned = True
        while cleaned:
            for ind, tmp in enumerate(tmp_stoich):
                if tmp == '][':
                    del tmp_stoich[ind]
                    tmp_stoich.insert(ind, '[')
                    tmp_stoich.insert(ind, ']')
                    cleaned = True
                elif ind == len(tmp_stoich) - 1:
                    cleaned = False
        for ind, tmp in enumerate(tmp_stoich):
            if tmp == '[':
                end_bracket = False
                while not end_bracket:
                    if tmp_stoich[ind + 1] == ']':
                        end_bracket = True
                    tmp_stoich[ind] += tmp_stoich[ind + 1]
                    del tmp_stoich[ind + 1]
        try:
            tmp_stoich.remove(']')
        except:
            pass
        try:
            tmp_stoich.remove('')
        except:
            pass
        elements = tmp_stoich
    return elements


def get_root_source(source):
    """ Get the main file source from a doc's source list.

    Parameters:
        source (list): contents of doc['source'].

    Returns:
        str: "root" filename, e.g. if source = ['KP.cell', 'KP.param',
            'KP_specific_structure.res'] then root = 'KP_specific_structure'.

    """
    src_list = set()
    for src in source:
        if src.endswith('.res') or src.endswith('.castep') or src.endswith('.history'):
            src_list.add('.'.join(src.split('/')[-1].split('.')[0:-1]))
        elif 'OQMD' in src.upper():
            src_list.add(src)

    if len(src_list) > 1:
        raise RuntimeError('Ambiguous root source')

    return list(src_list)[0]


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
                            form += elem[0] + '$_\\mathrm{' + str(int(elem[1])) + '}$'
                        else:
                            form += elem[0] + str(int(elem[1]))
        assert form != ''
    else:
        for elem in sorted(stoich):
            if elem[1] == 1:
                form += elem[0]
            elif int(elem[1]) != 0:
                if tex:
                    form += elem[0] + '$_\\mathrm{' + str(int(elem[1])) + '}$'
                else:
                    form += elem[0] + str(int(elem[1]))
    return form
