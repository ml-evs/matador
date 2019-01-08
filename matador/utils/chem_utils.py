# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule defines some useful chemical functions and
constants, with a focus on battery materials.

"""


import copy
import numpy as np

# global consts
FARADAY_CONSTANT_Cpermol = 96.485332e3
Cperg_to_mAhperg = 2.778e-1
C_TO_mAh = Cperg_to_mAhperg
HARTREE_TO_EV = 27.21139
BOHR_TO_ANGSTROM = 0.529177211
RY_TO_EV = 13.605693009
KBAR_TO_GPA = 0.1
eV_PER_ANGSTROM_CUBED_TO_GPa = 160.21776
AVOGADROS_NUMBER = 6.022141e23
ANGSTROM_CUBED_TO_CENTIMETRE_CUBED = 1e-24
ELECTRON_CHARGE = 1.6021766e-19
KELVIN_TO_EV = 8.61733e-5

EPS = 1e-12


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
    periodic_table['Lan'] = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    periodic_table['Act'] = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
    periodic_table['X'] = [elem for group in periodic_table for elem in periodic_table[group]]
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


def get_concentration(doc, elements, include_end=False):
    """ Returns x for A_x B_{1-x}
    or x,y for A_x B_y C_z, (x+y+z=1).

    Parameters:
        doc (list/dict): structure to evaluate OR matador-style stoichiometry.
        elements (list): list of element symbols to enforce ordering.

    Keyword arguments:
        include_end (bool): whether or not to return the final value, i.e.
            [x, y, z] rather than [x, y] in the above.

    Returns:
        list of float: concentrations of elements in given order.

    """
    if isinstance(doc, dict):
        if doc.get('stoichiometry') is None:
            raise RuntimeError('No stoichiometry found.')

        stoich = doc['stoichiometry']
    else:
        stoich = doc

    concs = [0.0] * (len(elements) - bool(not include_end))
    for _, elem in enumerate(stoich):
        if (include_end and elem[0] in elements) or (not include_end and elem[0] in elements[:-1]):
            concs[elements.index(elem[0])] = elem[1] / float(get_atoms_per_fu(doc))
    return concs


def get_num_intercalated(cursor):
    """ Return array of the number of intercalated atoms
    per host atom from a list of structures, of type defined by
    the first entry in the structures' concentration vectors.

    Parameters:
        cursor (list of dict): structures to evaluate.

    Returns:
        ndarray: number of intercalated ions in each structure.

    """
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

    Parameters:
        initial_doc (dict): matador doc of delithiated phase
        final_doc (dict): matador doc of maximally lithiated phase

    Returns:
       volumetric_capacity (float): capacity in mAh/cm^3.

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
    volume_per_fu_cm3 = initial_doc['cell_volume'] * ANGSTROM_CUBED_TO_CENTIMETRE_CUBED / initial_doc['num_fu']
    return ((num_ions_per_initial_fu / volume_per_fu_cm3) * (ELECTRON_CHARGE * Cperg_to_mAhperg))


def get_atoms_per_fu(doc):
    """ Calculate and return the number of atoms per formula unit.

    Parameters:
        doc (list/dict): structure to evaluate OR matador-style stoichiometry.

    """
    if isinstance(doc, dict):
        return sum([elem[1] for elem in doc['stoichiometry']])
    else:
        return sum([elem[1] for elem in doc])


def get_formation_energy(chempots, doc, energy_key='enthalpy_per_atom', temperature=None):
    """ From given chemical potentials, calculate the simplest
    formation energy per atom of the desired document.

    Parameters:
        chempots (list of dict): list of chempot structures.
        doc (dict): structure to evaluate.

    Keyword arguments:
        energy_key (str): name of energy field to use to calculate formation energy.
        temperature (float): if not None, use doc[energy_key][temperature] to calculate formation energy.

    Returns:
        float: formation energy per atom.

    """
    if temperature is not None:
        formation = doc[energy_key][temperature]
    else:
        formation = doc[energy_key]

    num_chempots = get_number_of_chempots(doc, chempots)
    num_atoms_per_fu = get_atoms_per_fu(doc)
    for ind, mu in enumerate(chempots):
        num_atoms_per_mu = get_atoms_per_fu(mu)
        if temperature is not None:
            formation -= mu[energy_key][temperature] * num_chempots[ind] * num_atoms_per_mu / num_atoms_per_fu
        else:
            formation -= mu[energy_key] * num_chempots[ind] * num_atoms_per_mu / num_atoms_per_fu
    return formation


def get_number_of_chempots(stoich, chempot_stoichs):
    """ Return the required number of each (arbitrary) chemical potentials
    to construct one formula unit of the input stoichiometry.

    Parameters:
        stoich (list/dict): matador-style stoichiometry,
            e.g. [['Li', 3], ['P', 1]], or the full document.
        chempot_stoichs (list/dict): list of stoichiometries of the input
            chemical potentials, or the full documents.

    Returns:
        list: number of each chemical potential required to create
            1 formula unit.

    Raises:
        RuntimeError: if the stoichiometry provided cannot be created
            with the given chemical potentials.

    """
    import scipy.linalg

    if isinstance(stoich, dict):
        stoich = stoich['stoichiometry']
    if isinstance(chempot_stoichs[0], dict):
        chempot_stoichs = [mu['stoichiometry'] for mu in chempot_stoichs]

    # find all elements present in the chemical potentials
    elements = set()
    for mu in chempot_stoichs:
        for elem, num in mu:
            elements.add(elem)
    elements = sorted(list(elements))

    chempot_matrix = np.asarray([get_padded_composition(mu, elements) for mu in chempot_stoichs])
    num_extraneous_equations = max(np.shape(chempot_matrix)) - min(np.shape(chempot_matrix))

    try:
        solution = np.asarray(get_padded_composition(stoich, elements))
    except RuntimeError:
        raise RuntimeError('Stoichiometry {} could not be created from chemical potentials {}: missing chempot'
                           .format(stoich, chempot_stoichs))

    if num_extraneous_equations == 0:
        num_chempots = scipy.linalg.solve(chempot_matrix.T, solution)
    else:
        num_chempots = scipy.linalg.solve(chempot_matrix[:, :-num_extraneous_equations].T, solution[:-num_extraneous_equations])
        # check equations are consistent
        for i in range(1, num_extraneous_equations+1):
            verify = sum(num_chempots * chempot_matrix[:, -i])
            if np.abs(verify - solution[-i]) > EPS:
                raise RuntimeError('Stoichiometry {} could not be created from chemical potentials {}: stoichiometry inconsistent with chempots'
                                   .format(stoich, chempot_stoichs))

    return num_chempots.tolist()


def get_stoich(atom_types):
    """ Return integer stoichiometry from atom_types list.

    Parameters:
        atom_types (list): list of element symbols of each atom.

    Returns:
        list: matador-style stoichiometry, e.g. [['Li', 1], ['P', 2]].

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


def get_padded_composition(stoichiometry, elements):
    """ Return a list that contains how many of each species in
    elements exists in the given stoichiometry. e.g. for [['Li', 2], ['O', 1]]
    with elements ['O', 'Li', 'Ba'], this function will return [1, 2, 0].

    Parameters:
        stoichiometry (list): matador-style stoichiometry, as above.
        elements (list): order of element labels to pick out.

    """
    composition = []
    for element in elements:
        if not isinstance(element, str):
            raise RuntimeError('Found invalid element symbol {}'.format(element))
        for species in stoichiometry:
            if not isinstance(species, list):
                raise RuntimeError('Found invalid stoichiometry {}'.format(stoichiometry))
            if species[0] == element:
                composition.append(species[1])
                break
            elif species[0] not in elements:
                raise RuntimeError('Extra element {} in stoichiometry'.format(species[0]))
        else:
            composition.append(0)

    return composition


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
                ratio_dict[stoichiometry[i][0]
                           + stoichiometry[j][0]] = round(float(stoichiometry[i][1]) / stoichiometry[j][1], 3)
    return ratio_dict


def get_stoich_from_formula(formula: str, sort=True):
    """ Convert formula string, e.g. Li2TiP4 into a matador-style
    stoichiometry, e.g. [['Li', 2], ['Ti', 1], ['P', 4]].

    Parameters:
        formula (str): chemical formula of compound

    Returns:
        list: sorted matador-style stoichiometry.

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
            except(ValueError, IndexError):
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
    if sort:
        return sorted(stoich)
    else:
        return stoich


def parse_element_string(elements_str, stoich=False):
    """ Parse element query string with macros. Has to parse braces
    too, and throw an error if brackets are unmatched.

    e.g.
        Parameters: '[VII][Fe,Ru,Os][I]'
        Returns: ['[VII]', '[Fe,Ru,Os]', '[I]']

    e.g.2
        Parameters: '[VII]2[Fe,Ru,Os][I]'
        Returns: ['[VII]2', '[Fe,Ru,Os]', '[I]']

    Parameters:
        elements_str: str, chemical formula, including macros.

    Keyword arguments:
        stoich: bool, parse as a stoichiometry, i.e. check for numbers

    Raises:
        RuntimeError: if the composition contains unmatched brackets.

    Returns:
        list: split list of elements contained in input

    """
    import re
    valid = False
    for char in elements_str:
        if char not in ['[', ']', '{', '}', ',', ':'] and not char.isalnum():
            raise RuntimeError('Illegal character {} detected in query.'.format(char))
    valid = False
    for char in elements_str:
        if char.isupper():
            valid = True
            break
    if not valid:
        raise RuntimeError('Composition must contain at least one upper case character.')
    elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements_str) if elem]
    if stoich:
        tmp_stoich = elements
        for ind, strng in enumerate(elements):
            if not any(char.isdigit() for char in strng):
                tmp_stoich[ind] = [strng]
            else:
                tmp_stoich[ind] = [elem for elem in re.split(r'([0-9]+)', strng) if elem]
        elements = [item for sublist in tmp_stoich for item in sublist]
    # split macros
    while '[' in elements or '{' in elements or '][' in elements or '}{' in elements or ']{' in elements or '}[' in elements:
        tmp_stoich = list(elements)
        cleaned = True
        while cleaned:
            for ind, tmp in enumerate(tmp_stoich):
                if tmp == '][':
                    del tmp_stoich[ind]
                    tmp_stoich.insert(ind, '[')
                    tmp_stoich.insert(ind, ']')
                    cleaned = True
                elif tmp == '}{':
                    del tmp_stoich[ind]
                    tmp_stoich.insert(ind, '{')
                    tmp_stoich.insert(ind, '}')
                    cleaned = True
                elif tmp == ']{':
                    del tmp_stoich[ind]
                    tmp_stoich.insert(ind, '{')
                    tmp_stoich.insert(ind, ']')
                    cleaned = True
                elif tmp == '}[':
                    del tmp_stoich[ind]
                    tmp_stoich.insert(ind, '[')
                    tmp_stoich.insert(ind, '}')
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
        for ind, tmp in enumerate(tmp_stoich):
            if tmp == '{':
                end_bracket = False
                while not end_bracket:
                    if tmp_stoich[ind + 1] == '}':
                        end_bracket = True
                    tmp_stoich[ind] += tmp_stoich[ind + 1]
                    del tmp_stoich[ind + 1]
        if ']' in tmp_stoich:
            tmp_stoich.remove(']')
        if '}' in tmp_stoich:
            tmp_stoich.remove('}')
        if '' in tmp_stoich:
            tmp_stoich.remove('')

        elements = tmp_stoich

    for elem in elements:
        if '}[' in elem or ']{' in elem:
            raise RuntimeError('Unmatched brackets in query string')

    return elements


def get_root_source(source):
    """ Get the main file source from a doc's source list.

    Parameters:
        source (list/dict): contents of doc['source'] or the doc itself.

    Returns:
        str: "root" filename, e.g. if source = ['KP.cell', 'KP.param',
            'KP_specific_structure.res'] then root = 'KP_specific_structure'.

    """
    if isinstance(source, dict):
        source = source['source']
    src_list = set()
    for src in source:
        if src.endswith('.res') or src.endswith('.castep') or src.endswith('.history'):
            src_list.add('.'.join(src.split('/')[-1].split('.')[0:-1]))
        elif 'OQMD' in src.upper():
            src_list.add(src)
        elif src == 'command_line':
            src_list.add('command line')

    if len(src_list) > 1:
        raise RuntimeError('Ambiguous root source')
    if len(src_list) < 1:
        raise RuntimeError('Unable to find root source')

    return list(src_list)[0]


def get_formula_from_stoich(stoich, elements=None, tex=False, latex_sub_style=''):
    """ Get the chemical formula of a structure from
    its matador stoichiometry.

    Parameters:
        stoich (list): matador-style stoichiometry.

    Keyword arguments:
        elements (list): list of element symbols to enforce order.
        tex (bool): whether to print a LaTeX-compatibile string.
        latex_sub_style (str): a string to wrap subscripts in, e.g.
            r"\\mathrm" or r"\\text" (default is blank).

    Returns:
        str: the string representation of the chemical formula.

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
                            if elem[1] % 1 == 0:
                                elem[1] = int(elem[1])
                            form += r'{}$_{}{{{}}}$'.format(elem[0], latex_sub_style, elem[1])
                        else:
                            form += elem[0] + str(int(elem[1]))
        assert form != ''
    else:
        for elem in stoich:
            if elem[1] == 1:
                form += elem[0]
            elif int(elem[1]) != 0:
                if tex:
                    if elem[1] % 1 == 0:
                        elem[1] = int(elem[1])
                    form += '{}$_{}{{{}}}$'.format(elem[0], latex_sub_style, elem[1])
                else:
                    form += elem[0] + str(int(elem[1]))
    return form
