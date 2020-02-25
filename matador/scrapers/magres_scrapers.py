# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements some scraper functions for
NMR-related inputs and outputs, e.g. .magres files.

"""


from collections import defaultdict
from os import stat
from pwd import getpwuid

import numpy as np
from matador.utils.cell_utils import cart2abc, cart2frac
from matador.utils.chem_utils import get_stoich
from matador.scrapers.utils import scraper_function


@scraper_function
def magres2dict(seed, **kwargs):
    """ Extract available information from .magres file. Assumes units of
    Angstrom and ppm for relevant quantities.
    """
    magres = defaultdict(list)
    # read .pwout file into array
    seed = seed.replace('.magres', '')
    with open(seed + '.magres', 'r') as f:
        flines = f.readlines()
    # add .magres to source
    magres['source'].append(seed+'.magres')
    # grab file owner username
    try:
        magres['user'] = getpwuid(stat(seed + '.magres').st_uid).pw_name
    except Exception:
        magres['user'] = 'xxx'
        if kwargs.get('debug'):
            print(seed + '.magres has no owner.')

    for line_no, line in enumerate(flines):
        line = line.lower().strip()
        if line in ['<atoms>', '[atoms]']:
            i = 1
            while flines[line_no+i].strip().lower() not in ['</atoms>', '[/atoms]']:
                if i > 1e5:
                    raise RuntimeError
                if 'units' in flines[line_no+i]:
                    i += 1
                    continue
                elif 'lattice' in flines[line_no+i]:
                    lattice = flines[line_no+i].split()[1:]
                    for j in range(3):
                        magres['lattice_cart'].append([float(elem) for elem in lattice[j*3:(j+1)*3]])
                    magres['lattice_abc'] = cart2abc(magres['lattice_cart'])
                elif 'atom' in flines[line_no+i]:
                    atom = flines[line_no+i].split()
                    magres['atom_types'].append(atom[1])
                    magres['positions_abs'].append([float(elem) for elem in atom[-3:]])
                i += 1
            break

    magres['num_atoms'] = len(magres['atom_types'])
    magres['positions_frac'] = cart2frac(magres['lattice_cart'], magres['positions_abs'])
    magres['stoichiometry'] = get_stoich(magres['atom_types'])

    for line_no, line in enumerate(flines):
        line = line.lower().strip()
        if line in ['<magres>', '[magres]']:
            i = 1
            while flines[line_no+i].strip().lower() not in ['</magres>', '[/magres]']:
                if i > 1e5:
                    raise RuntimeError
                if 'units' in flines[line_no+i]:
                    i += 1
                    continue
                elif 'sus' in flines[line_no+i]:
                    sus = flines[line_no+i].split()[1:]
                    for j in range(3):
                        magres['susceptibility_tensor'].append([float(val) for val in sus[3*j:3*(j+1)]])
                elif 'ms' in flines[line_no+i]:
                    ms = flines[line_no+i].split()[3:]
                    magres['magnetic_shielding_tensor'].append([])
                    for j in range(3):
                        magres['magnetic_shielding_tensor'][-1].append([float(val) for val in ms[3*j:3*(j+1)]])
                    magres['chemical_shifts'].append(0)
                    magres['chemical_shift_anisos'].append(0)
                    magres['chemical_shift_asymmetries'].append(0)
                    for j in range(3):
                        magres['chemical_shifts'][-1] += magres['magnetic_shielding_tensor'][-1][j][j] / 3

                    # find eigenvalues of symmetric part of shielding and order them to calc anisotropy eta
                    symmetric_shielding = (
                        0.5 *
                        (magres['magnetic_shielding_tensor'][-1] + np.asarray(magres['magnetic_shielding_tensor'][-1]).T)
                    )
                    eig_vals, eig_vecs = np.linalg.eig(symmetric_shielding)
                    eig_vals, eig_vecs = zip(*sorted(zip(eig_vals, eig_vecs),
                                                     key=lambda eig: abs(eig[0] - magres['chemical_shifts'][-1])))
                    # Haeberlen convention: |s_zz - s_iso| >= |s_xx - s_iso| >= |s_yy - s_iso|
                    s_yy, s_xx, s_zz = eig_vals
                    s_iso = magres['chemical_shifts'][-1]
                    # convert from reduced anistropy to CSA
                    magres['chemical_shift_anisos'][-1] = s_zz - (s_xx + s_yy)/2.0
                    magres['chemical_shift_asymmetries'][-1] = (s_yy - s_xx) / (s_zz - s_iso)
                i += 1

    for line_no, line in enumerate(flines):
        line = line.lower().strip()
        if line in ['<calculation>', '[calculation]']:
            i = 1
            while flines[line_no+i].strip().lower() not in ['</calculation>', '[/calculation]']:
                if i > 1e5:
                    raise RuntimeError
                # space important as it excludes other calc_code_x variables
                if 'calc_code ' in flines[line_no+i]:
                    magres['calculator'] = flines[line_no+i].split()[1]
                if 'calc_code_version' in flines[line_no+i]:
                    magres['calculator_version'] = flines[line_no+i].split()[1]
                i += 1

    return magres, True
