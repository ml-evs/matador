# coding: utf-8
""" This file implements some scraper functions for
NMR-related inputs and outputs.
"""

# matador modules
from matador.utils.cell_utils import cart2abc, cart2frac
# external libraries
import bson.json_util as json
# standard library
from collections import defaultdict
from os import stat
from pwd import getpwuid
from traceback import print_exc


def magres2dict(seed, **kwargs):
    """ Extract available information from .magres file. Assumes units of
    Angstrom and ppm for relevant quantities.
    """
    magres = defaultdict(list)
    try:
        # read .pwout file into array
        seed = seed.replace('.magres', '')
        with open(seed + '.magres', 'r') as f:
            flines = f.readlines()
        # add .magres to source
        magres['source'].append(seed+'.magres')
        # grab file owner username
        try:
            magres['user'] = getpwuid(stat(seed + '.magres').st_uid).pw_name
        except:
            if kwargs.get('debug'):
                print(seed + '.magres has no owner.')
            magres['user'] == 'xxx'

        for line_no, line in enumerate(flines):
            if '<atoms>' in line.lower():
                i = 1
                while '</atoms>' not in flines[line_no+i].lower():
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

        magres['num_atoms'] = len(magres['atom_types'])
        magres['positions_frac'] = cart2frac(magres['lattice_cart'], magres['positions_abs'])

        for line_no, line in enumerate(flines):
            if '<magres>' in line.lower():
                i = 1
                while '</magres>' not in flines[line_no+i].lower():
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
                        for j in range(3):
                            magres['chemical_shifts'][-1] += magres['magnetic_shielding_tensor'][-1][j][j] / 3
                    i += 1

        for line_no, line in enumerate(flines):
            if '<calculation>' in line.lower():
                i = 1
                while '</calculation>' not in flines[line_no+i].lower():
                    if i > 1e5:
                        raise RuntimeError
                    if 'calc_code' in flines[line_no+i] and 'version' not in flines[line_no+i]:
                        magres['calculator'] = flines[line_no+i].split()[1]
                    i += 1

    except Exception as oops:
        if kwargs.get('verbosity') is not None:
            if kwargs.get('verbosity') > 0:
                print_exc()
                print('Error in .pwout file', seed + '.pwout, skipping...')
        if type(oops) == IOError:
            print_exc()
        return seed+'.magres\t\t' + str(type(oops)) + ' ' + str(oops) + '\n', False

    if kwargs.get('verbosity') is not None:
        if kwargs.get('verbosity') > 4:
            print(json.dumps(magres, indent=2))
    return magres, True
