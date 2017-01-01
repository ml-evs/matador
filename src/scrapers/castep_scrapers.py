#!/usr/bin/python
# coding: utf-8
""" This file implements the scraper functions for CASTEP-related
inputs and outputs.
"""
from __future__ import print_function
# matador modules
from utils.cell_utils import abc2cart, calc_mp_spacing
# external libraries
try:
    import bson.json_util as json
except:
    pass
# standard library
from collections import defaultdict
from os import getcwd, stat
from os.path import isfile
from time import strptime
from math import gcd
from pwd import getpwuid
from traceback import print_exc
import glob
import gzip


def res2dict(seed, db=True, verbosity=0, **kwargs):
    """ Extract available information from .res file; preferably
    used in conjunction with cell or param file.
    """
    # use defaultdict to allow for easy appending
    res = defaultdict(list)
    try:
        # read .res file into array
        if seed.endswith('.res'):
            seed = seed.replace('.res', '')
        with open(seed+'.res', 'r') as f:
            flines = f.readlines()
        # add .res to source
        res['source'].append(seed+'.res')
        # grab file owner username
        try:
            res['user'] = getpwuid(stat(seed+'.res').st_uid).pw_name
        except:
            if kwargs.get('debug'):
                print(seed+'.res has no owner.')
            res['user'] == 'xxx'
        if 'CollCode' in seed:
            res['icsd'] = seed.split('CollCode')[-1]
        # alias special lines in res file
        titl = ''
        cell = ''
        remark = ''
        for line in flines:
            if 'TITL' in line and db:
                # if not db, then don't read title
                titl = line.split()
                if len(titl) != 12:
                    raise RuntimeError('missing some TITL info')
            elif 'CELL' in line:
                cell = line.split()
            elif 'REM' in line:
                remark = line.split()
        if cell == '':
            raise RuntimeError('missing CELL info')
        elif titl == '' and db:
            raise RuntimeError('missing TITL')
        if db:
            res['pressure'] = float(titl[2])
            res['cell_volume'] = float(titl[3])
            res['enthalpy'] = float(titl[4])
            res['num_atoms'] = int(titl[7])
            res['space_group'] = titl[8].strip('()')
            res['enthalpy_per_atom'] = res['enthalpy'] / res['num_atoms']
            res['total_energy'] = res['enthalpy'] - res['pressure']*res['cell_volume']
            res['total_energy_per_atom'] = res['total_energy'] / res['num_atoms']
        res['lattice_abc'] = [list(map(float, cell[2:5])), list(map(float, cell[5:8]))]
        # calculate lattice_cart from abc
        res['lattice_cart'] = abc2cart(res['lattice_abc'])
        for line_no, line in enumerate(flines):
            if 'SFAC' in line:
                i = 1
                while 'END' not in flines[line_no+i] and line_no+i < len(flines):
                    cursor = flines[line_no+i].split()
                    res['atom_types'].append(cursor[0])
                    res['positions_frac'].append(list(map(float, cursor[2:5])))
                    i += 1
        # deal with implicit encapsulation
        if len(remark) > 0:
            if 'NTPROPS' in remark:
                res['cnt_chiral'] = [0, 0]
                res['encapsulated'] = True
                for ind, entry in enumerate(remark):
                    if 'chiralN' in entry:
                        res['cnt_chiral'][0] = int(remark[ind+1].replace(',', ''))
                    if 'chiralM' in entry:
                        res['cnt_chiral'][1] = int(remark[ind+1].replace(',', ''))
                    if entry == '\'r\':':
                        res['cnt_radius'] = float(remark[ind+1].replace(',', ''))
                    if entry == '\'z\':':
                        temp_length = remark[ind+1].replace(',', '')
                        temp_length = temp_length.replace('\n', '')
                        temp_length = temp_length.replace('}', '')
                        res['cnt_length'] = float(temp_length)
        # calculate stoichiometry
        res['stoichiometry'] = defaultdict(float)
        for atom in res['atom_types']:
            if atom not in res['stoichiometry']:
                res['stoichiometry'][atom] = 0
            res['stoichiometry'][atom] += 1
        gcd_val = 0
        for atom in res['atom_types']:
            if gcd_val == 0:
                gcd_val = res['stoichiometry'][atom]
            else:
                gcd_val = gcd(res['stoichiometry'][atom], gcd_val)
        # convert stoichiometry to tuple for fryan
        temp_stoich = []
        for key, value in res['stoichiometry'].iteritems():
            if float(value)/gcd_val % 1 != 0:
                temp_stoich.append([key, float(value)/gcd_val])
            else:
                temp_stoich.append([key, value/gcd_val])
        res['stoichiometry'] = temp_stoich
        atoms_per_fu = 0
        for elem in res['stoichiometry']:
            atoms_per_fu += elem[1]
        res['num_fu'] = len(res['atom_types']) / atoms_per_fu
    except Exception as oops:
        if verbosity > 0:
            print_exc()
            print('Error in .res file', seed + '.res, skipping...')
        if type(oops) == IOError:
            print_exc()
        return seed+'.res\t\t' + str(type(oops)) + ' ' + str(oops) + '\n', False
    if verbosity > 4:
        print(json.dumps(res, indent=2))
    return res, True


def cell2dict(seed, db=True, verbosity=0, **kwargs):
    """ Extract available information from .cell file; probably
    to be merged with another dict from a .param or .res file.
    """
    cell = defaultdict(list)
    try:
        if seed.endswith('.cell'):
            seed = seed.replace('.cell', '')
        with open(seed+'.cell', 'r') as f:
            flines = f.readlines()
        # add cell file to source
        cell['source'].append(seed+'.cell')
        for line_no, line in enumerate(flines):
            if line.startswith(('#', '!')):
                continue
            elif '%block lattice_cart' in line.lower() and not db:
                cell['lattice_cart'] = []
                i = 1
                while 'endblock' not in flines[line_no+i].lower():
                    if not flines[line_no+i].strip()[0].isalpha():
                        cell['lattice_cart'].append(list(map(float, flines[line_no+i].split())))
                    i += 1
                assert(len(cell['lattice_cart']) == 3)
            elif '%block species_pot' in line.lower():
                cell['species_pot'] = dict()
                i = 1
                while 'endblock' not in flines[line_no+i].lower():
                    cell['species_pot'][flines[line_no+i].split()[0]] = \
                        flines[line_no+i].split()[1].split('/')[-1]
                    cell['species_pot'][flines[line_no+i].split()[0]] = \
                        cell['species_pot'][flines[line_no+i].split()[0]].replace('compat7', '')
                    cell['species_pot'][flines[line_no+i].split()[0]] = \
                        cell['species_pot'][flines[line_no+i].split()[0]].replace(',', '')
                    cell['species_pot'][flines[line_no+i].split()[0]] = \
                        cell['species_pot'][flines[line_no+i].split()[0]].replace('()', '')
                    i += 1
            elif '%block cell_constraints' in line.lower():
                cell['cell_constraints'] = list()
                for j in range(2):
                    cell['cell_constraints'].append(list(map(int, flines[line_no+j+1].split())))
            elif '%block hubbard_u' in line.lower():
                cell['hubbard_u'] = defaultdict(list)
                i = 0
                while 'endblock' not in flines[line_no+i].lower():
                    line = flines[line_no+i]
                    if line == 'eV' or len(line.split()) < 3:
                        i += 1
                        continue
                    else:
                        atom = line.split()[0]
                        orbital = line.split()[1].replace(':', '')
                        shift = float(line.split()[-1])
                        atom = line.split()[0]
                        cell['hubbard_u'][atom] = dict()
                        cell['hubbard_u'][atom][orbital] = shift
                        i += 1
            elif '%block external_pressure' in line.lower():
                cell['external_pressure'] = list()
                for j in range(3):
                    flines[line_no+j+1] = flines[line_no+j+1].replace(',', '')
                    cell['external_pressure'].append(list(map(float, flines[line_no+j+1].split())))
            elif 'kpoints_mp_spacing' in line.lower() or 'kpoint_mp_spacing' in line.lower():
                if 'spectral_kpoints_mp_spacing' in line.lower() or 'spectral_kpoint_mp_spacing' in line.lower():
                    cell['spectral_kpoints_mp_spacing'] = float(line.split()[-1])
                else:
                    cell['kpoints_mp_spacing'] = float(line.split()[-1])
            elif 'kpoints_mp_grid' in line.lower() or 'kpoint_mp_grid' in line.lower():
                if 'spectral_kpoints_mp_grid' in line.lower() or 'spectral_kpoint_mp_grid' in line.lower():
                    cell['spectral_kpoints_mp_grid'] = list(map(int, line.split()[-3:]))
                else:
                    cell['kpoints_mp_grid'] = list(map(int, line.split()[-3:]))
            elif 'kpoints_mp_offset' in line.lower() or 'kpoint_mp_offset' in line.lower():
                if 'spectral_kpoints_mp_offset' in line.lower() or 'spectral_kpoint_mp_offset' in line.lower():
                    cell['spectral_kpoints_mp_offset'] = list(map(float, line.split()[-3:]))
                else:
                    cell['kpoints_mp_offset'] = list(map(float, line.split()[-3:]))
            elif not db:
                if '%block positions_frac' in line.lower():
                    atomic_init_spins = defaultdict(list)
                    i = 0
                    while '%endblock positions_frac' not in flines[line_no+i].lower():
                        if 'spin=' in flines[line_no+i].lower():
                            split_line = flines[line_no+i].split()
                            atomic_init_spins[split_line[0]] = \
                                split_line[-1].lower().replace('spin=', '')
                        i += 1
                    if atomic_init_spins:
                        cell['atomic_init_spins'] = atomic_init_spins
                elif 'fix_com' in line.lower():
                    cell['fix_com'] = line.split()[-1]
                elif 'symmetry_generate' in line.lower():
                    cell['symmetry_generate'] = True
                elif 'symmetry_tol' in line.lower():
                    cell['symmetry_tol'] = float(line.split()[-1])
                elif 'snap_to_symmetry' in line.lower():
                    cell['snap_to_symmetry'] = True
                elif 'quantisation_axis' in line.lower():
                    cell['quantisation_axis'] = list(map(int, line.split()[1:]))
        if 'external_pressure' not in cell:
            cell['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0], [0.0]]
    except Exception as oops:
        if verbosity > 0:
            print_exc()
            print('Error in', seed + '.cell, skipping...')
        return seed + '\t\t' + str(type(oops)) + ' ' + str(oops), False
    try:
        for species in cell['species_pot']:
            if 'OTF' in cell['species_pot'][species].upper():
                pspot_seed = ''
                for dir in seed.split('/')[:-1]:
                    pspot_seed += dir + '/'
                # glob for all .usp files with format species_*OTF.usp
                pspot_seed += species + '_*OTF.usp'
                for globbed in glob.glob(pspot_seed):
                    if isfile(globbed):
                        with open(globbed, 'r') as f:
                            flines = f.readlines()
                            for line_no, line in enumerate(flines):
                                if 'Pseudopotential Report' in line:
                                    i = 0
                                    while i+line_no < len(flines)-3:
                                        if 'Pseudopotential Report' in flines[line_no+i]:
                                            i += 2
                                            elem = flines[line_no+i].split(':')[1].split()[0]
                                        elif 'core correction' in flines[line_no+i]:
                                            i += 2
                                            cell['species_pot'][elem] = \
                                                flines[line_no+i].split('"')[1]
                                        i += 1
    except Exception as oops:
        if verbosity > 0:
            print_exc()
            print('Error in', seed + '.cell, skipping...')
        if type(oops) == IOError:
            print_exc()
        return seed + '.cell\t\t' + str(type(oops)) + ' ' + str(oops), False
    if kwargs.get('debug'):
        print(json.dumps(cell, indent=2))
    return cell, True


def param2dict(seed, db=True, verbosity=0, **kwargs):
    """ Extract available information from .param file; probably
    to be merged with other dicts from other files.

    seed  : filename with or without file extension
    db    : if True, only scrape relevant info, otherwise scrape all

    """
    try:
        param = defaultdict(list)
        if seed.endswith('.param'):
            seed = seed.replace('.param', '')
        with open(seed+'.param', 'r') as f:
            flines = f.readlines()
        param['source'].append(seed+'.param')
        # exclude some useless info if importing to db
        scrub_list = ['checkpoint', 'write_bib', 'mix_history_length',
                      'fix_occupancy', 'page_wvfns', 'num_dump_cycles',
                      'backup_interval', 'geom_max_iter', 'fixed_npw',
                      'write_cell_structure', 'bs_write_eigenvalues',
                      'calculate_stress', 'opt_strategy', 'max_scf_cycles']
        false_str = ['False', 'false', '0']
        splitters = [':', '=', '\t', ' ']
        for line_no, line in enumerate(flines):
            line = line.lower()
            # skip blank lines and comments
            if line.startswith(('#', '!')) or len(line.strip()) == 0:
                continue
            else:
                # if scraping to db, ignore "rubbish"
                if db:
                    if [rubbish for rubbish in scrub_list if rubbish in line]:
                        continue
                # read all other parameters in
                for splitter in splitters:
                    if splitter in line:
                        param[line.split(splitter)[0].strip()] = \
                            line.split(splitter)[-1].strip()
                        if 'spin_polarized' in line:
                            if [false for false in false_str
                                    if false in param['spin_polarized']]:
                                param['spin_polarized'] = False
                            else:
                                param['spin_polarized'] = True
                        if 'cut_off_energy' in line and 'mix_cut_off_energy' not in line:
                            temp_cut_off = (param['cut_off_energy'].replace('ev', ''))
                            temp_cut_off = temp_cut_off.strip()
                            param['cut_off_energy'] = float(temp_cut_off)
                        if 'xc_functional' in line:
                            param['xc_functional'] = param['xc_functional'].upper()
                        if 'perc_extra_bands' in line:
                            param['perc_extra_bands'] = float(param['perc_extra_bands'])
                        if 'geom_force_tol' in line:
                            param['geom_force_tol'] = float(param['geom_force_tol'])
                        break
    except Exception as oops:
        if verbosity > 0:
            print_exc()
            print('Error in', seed+'.param, skipping...')
        if type(oops) == IOError:
            print_exc()
        return seed + '.param\t\t' + str(type(oops)) + ' ' + str(oops), False
    if kwargs.get('debug'):
        print(json.dumps(param, indent=2))
    return param, True


def dir2dict(seed, verbosity=0, **kwargs):
    """ Try to extract information from directory name; last hope
    if no param file has been found.
    """
    try:
        dir_dict = defaultdict(list)
        info = False
        if seed == '.':
            seed = getcwd().split('/')[-1]
        dirs_as_list = seed.split('/')
        task_list = ['GO', 'NMR', 'OPTICS']
        phase_list = ['alpha', 'beta', 'gamma', 'theta']
        defect_list = ['vacancy', 'interstitial']
        for dir in dirs_as_list:
            if dir == '.':
                dir = getcwd().split('/')[-1]
            if len(dir.split('-')) > 3:
                if dir[0].isalpha():
                    offset = 1
                else:
                    offset = 0
                try:
                    dir_dict['cut_off_energy'] = float(dir.split('-')[offset])
                    dir_dict['kpoints_mp_spacing'] = float(dir.split('-')[offset+1])
                    dir_dict['xc_functional'] = dir.split('-')[offset+4].upper()
                except:
                    pass
                if dir_dict['xc_functional'][0] == 'S':
                    dir_dict['xc_functional'] = dir_dict['xc_functional'][1:]
                    dir_dict['spin_polarized'] = True
                if [task for task in task_list if task in dir.split('-')[offset+5]]:
                    dir_dict['task'] = dir.split('-')[offset+5]
                info = True
            elif 'GPa' in dir:
                try:
                    dir_dict['external_pressure'].append([float(dir.split('_')[0]), 0.0, 0.0])
                    dir_dict['external_pressure'].append([float(dir.split('_')[0]), 0.0])
                    dir_dict['external_pressure'].append([float(dir.split('_')[0])])
                    info = True
                except:
                    pass
            if [phase for phase in phase_list if phase in dir]:
                for phase in phase_list:
                    if phase in dir:
                        dir_dict['phase'] = phase
                        info = True
            if [defect for defect in defect_list if defect in dir]:
                for defect in defect_list:
                    if defect in dir:
                        dir_dict['defect'] = defect
                        info = True
        if 'external_pressure' not in dir_dict:
            dir_dict['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0], [0.0]]
        if info:
            dir_dict['source'].append(seed)
        else:
            if verbosity > 0:
                print('No information found in dirname', seed)
    except Exception as oops:
        if verbosity > 0:
            print_exc()
            print('Error wrangling dir name', seed)
        return seed + '\t\t' + str(type(oops)) + ' ' + str(oops), False
    if kwargs.get('debug'):
        print(json.dumps(dir_dict, indent=2))
    return dir_dict, True


def castep2dict(seed, db=True, verbosity=0, **kwargs):
    """ From seed filename, create dict of the most relevant
    information about a calculation.
    """
    # use defaultdict to allow for easy appending
    castep = defaultdict(list)
    try:
        # read .castep, .history or .history.gz file
        if '.gz' in seed:
            with gzip.open(seed, 'r') as f:
                flines = f.readlines()
        else:
            with open(seed, 'r') as f:
                flines = f.readlines()
        # set source tag to castep file
        castep['source'].append(seed)
        pspot_report_dict = dict()
        # grab file owner
        castep['user'] = getpwuid(stat(seed).st_uid).pw_name
        if 'CollCode' in seed:
            temp_icsd = seed.split('CollCode')[-1].replace('.castep', '').replace('.history', '')
            castep['icsd'] = temp_icsd
        # wrangle castep file for parameters in 3 passes:
        # once forwards to get number and types of atoms
        # once backwards to get the final parameter set for the calculation
        # once more forwards, from the final step, to get the final structure
        for line_no, line in enumerate(flines):
            if 'atom types' not in castep and 'Cell Contents' in line:
                castep['atom_types'] = list()
                castep['positions_frac'] = list()
                i = 1
                atoms = False
                while True:
                    if atoms:
                        if 'xxxxxxxxx' in flines[line_no+i]:
                            atoms = False
                            break
                        else:
                            castep['atom_types'].append(flines[line_no+i].split()[1])
                            castep['positions_frac'].append(list(map(float,
                                                            (flines[line_no+i].split()[3:6]))))
                    if 'x------' in flines[line_no+i]:
                        atoms = True
                    i += 1
                castep['num_atoms'] = len(castep['atom_types'])
                castep['stoichiometry'] = defaultdict(float)
                for atom in castep['atom_types']:
                    if atom not in castep['stoichiometry']:
                        castep['stoichiometry'][atom] = 0
                    castep['stoichiometry'][atom] += 1
                gcd_val = 0
                for atom in castep['atom_types']:
                    if gcd_val == 0:
                        gcd_val = castep['stoichiometry'][atom]
                    else:
                        gcd_val = gcd(castep['stoichiometry'][atom], gcd_val)
                # convert stoichiometry to tuple for fryan
                temp_stoich = []
                for key in castep['stoichiometry']:
                    value = castep['stoichiometry'][key]
                    if float(value)/gcd_val % 1 != 0:
                        temp_stoich.append([key, float(value)/gcd_val])
                    else:
                        temp_stoich.append([key, value/gcd_val])
                castep['stoichiometry'] = temp_stoich
                atoms_per_fu = 0
                for elem in castep['stoichiometry']:
                    atoms_per_fu += elem[1]
                castep['num_fu'] = castep['num_atoms'] / atoms_per_fu
                break
        for line_no, line in enumerate(reversed(flines)):
            line_no = len(flines) - 1 - line_no
            if 'task' not in castep and 'type of calculation' in line:
                castep['task'] = line.split(':')[-1].strip().replace(" ", "")
            elif 'xc_functional' not in castep and 'functional' in line:
                # convert from .castep file xc_functional to param style
                xc_string = line.split(':')[-1].strip()
                if 'Local Density Approximation' in xc_string:
                    castep['xc_functional'] = 'LDA'
                elif 'Perdew Burke Ernzerhof' in xc_string:
                    castep['xc_functional'] = 'PBE'
                elif 'PBE for solids' in xc_string:
                    castep['xc_functional'] = 'PBESol'
                elif 'hybrid B3LYP' in xc_string:
                    castep['xc_functional'] = 'B3LYP'
                elif 'hybrid HSE03' in xc_string:
                    castep['xc_functional'] = 'HSE03'
                elif 'hybrid HSE06' in xc_string:
                    castep['xc_functional'] = 'HSE06'
            elif 'cut_off_energy' not in castep and 'plane wave basis set' in line:
                castep['cut_off_energy'] = float(line.split(':')[-1].split()[0])
            elif 'finite_basis_corr' not in castep and 'finite basis set correction  ' in line:
                castep['finite_basis_corr'] = line.split(':')[-1].strip()
            elif 'kpoints_mp_grid' not in castep and 'MP grid size for SCF' in line:
                castep['kpoints_mp_grid'] = list(map(int, list(line.split('is')[-1].split())))
            elif 'sedc_apply' not in castep and \
                    'DFT+D: Semi-empirical dispersion correction    : on' in line:
                castep['sedc_apply'] = True
                castep['sedc_scheme'] = flines[line_no+1].split(':')[1].split()[0]
            elif 'space_group' not in castep and 'Space group of crystal' in line:
                castep['space_group'] = line.split(':')[-1].split(',')[0].strip().replace(" ", "")
            elif 'external_pressure' not in castep and 'External pressure/stress' in line:
                try:
                    castep['external_pressure'] = list()
                    castep['external_pressure'].append(list(map(float, flines[line_no+1].split())))
                    castep['external_pressure'].append(list(map(float, flines[line_no+2].split())))
                    castep['external_pressure'].append(list(map(float, flines[line_no+3].split())))
                except:
                    if kwargs.get('dryrun') or verbosity > 0:
                        print_exc()
                    castep['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
            elif 'spin_polarized' not in castep and 'treating system as spin-polarized' in line:
                castep['spin_polarized'] = True
            elif 'hubbard_u' not in castep and 'Hubbard U values are eV' in line:
                castep['hubbard_u'] = defaultdict(list)
                i = 5
                while i < castep['num_atoms']:
                    line = flines[line_no+i].strip()
                    atom = line.split()[0].replace('|', '')
                    shifts = list(map(float, line.split()[-5:-1]))
                    for ind, shift in enumerate(shifts):
                        if shift != 0:
                            if atom not in castep['hubbard_u']:
                                castep['hubbard_u'][atom] = dict()
                            if ind == 0:
                                orbital = 's'
                            elif ind == 1:
                                orbital = 'p'
                            elif ind == 2:
                                orbital = 'd'
                            elif ind == 3:
                                orbital = 'f'
                            castep['hubbard_u'][atom][orbital] = shift
                    i += 1
            elif 'Pseudopotential Report' in line:
                if 'species_pot' not in castep:
                    castep['species_pot'] = dict()
                i = 0
                while i+line_no < len(flines)-3:
                    if 'Pseudopotential Report' in flines[line_no+i]:
                        i += 2
                        elem = flines[line_no+i].split(':')[1].split()[0]

                    elif 'core correction' in flines[line_no+i]:
                        i += 2
                        if not pspot_report_dict[elem]:
                            castep['species_pot'][elem] = flines[line_no+i].split('"')[1]
                            pspot_report_dict[elem] = True
                    i += 1
            elif 'species_pot' not in castep and 'Files used for pseudopotentials' in line:
                if 'species_pot' not in castep:
                    castep['species_pot'] = dict()
                i = 1
                while True:
                    if len(flines[line_no+i].strip()) == 0:
                        break
                    else:
                        elem = flines[line_no+i].split()[0].strip()
                        if not pspot_report_dict.get(elem):
                            castep['species_pot'][elem] = flines[line_no+i].split()[1].split('/')[-1]
                            if castep['species_pot'][elem] == 'Pseudopotential':
                                castep['species_pot'][elem] = flines[line_no+i].split()[0].strip()+'_OTF.usp'
                            pspot_report_dict[elem] = False
                    i += 1
        # write zero pressure if not found in file
        if 'external_pressure' not in castep:
            castep['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0], [0.0]]
        if 'spin_polarized' not in castep:
            castep['spin_polarized'] = False
        # task specific options
        if db and castep['task'] != 'geometryoptimization' and \
                castep['task'] != 'geometry optimization':
            raise RuntimeError('CASTEP file does not contain GO calculation')
        else:
            if castep['task'].strip() == 'geometryoptimization':
                final = False
                finish_line = 0
                castep['optimised'] = False
                success_string = 'Geometry optimization completed successfully'
                failure_string = 'Geometry optimization failed to converge after'
                for line_no, line in enumerate(flines):
                    if any(finished in line for finished in [success_string, failure_string]):
                        for line_next in range(line_no, len(flines)):
                            if any(finished in flines[line_next] for
                                    finished in [success_string, failure_string]):
                                finish_line = line_next
                                if success_string in flines[line_next]:
                                    castep['optimised'] = True
                                elif failure_string in flines[line_next]:
                                    castep['optimised'] = False
                        final = True
            elif not db:
                final = True
                finish_line = 0
            if final:
                final_flines = flines[finish_line+1:]
                for line_no, line in enumerate(final_flines):
                    if 'Real Lattice' in line:
                        castep['lattice_cart'] = list()
                        i = 1
                        while True:
                            if len(final_flines[line_no+i].strip()) == 0:
                                break
                            else:
                                temp_line = final_flines[line_no+i].split()[0:3]
                                castep['lattice_cart'].append(list(map(float, temp_line)))
                            i += 1
                    elif 'Lattice parameters' in line:
                        castep['lattice_abc'] = list()
                        i = 1
                        castep['lattice_abc'].append(
                            list(map(float,
                                     [final_flines[line_no+i].split('=')[1].strip().split(' ')[0],
                                      final_flines[line_no+i+1].split('=')[1].strip().split(' ')[0],
                                      final_flines[line_no+i+2].split('=')[1].strip().split(' ')[0]])))
                        castep['lattice_abc'].append(
                            list(map(float,
                                     [final_flines[line_no+i].split('=')[-1].strip(),
                                      final_flines[line_no+i+1].split('=')[-1].strip(),
                                      final_flines[line_no+i+2].split('=')[-1].strip()])))
                    elif 'Current cell volume' in line:
                        castep['cell_volume'] = float(line.split('=')[1].split()[0].strip())
                    elif 'Cell Contents' in line:
                        castep['positions_frac'] = list()
                        i = 1
                        atoms = False
                        while True:
                            if atoms:
                                if 'xxxxxxxxx' in final_flines[line_no+i]:
                                    atoms = False
                                    break
                                else:
                                    temp_frac = final_flines[line_no+i].split()[3:6]
                                    castep['positions_frac'].append(list(map(float, temp_frac)))
                            if 'x------' in final_flines[line_no+i]:
                                atoms = True
                            i += 1
                    # don't check if final_energy exists, as this will update for each GO step
                    elif 'Final energy, E' in line:
                        castep['total_energy'] = float(line.split('=')[1].split()[0])
                        castep['total_energy_per_atom'] = castep['total_energy'] / castep['num_atoms']
                    elif 'Final free energy' in line:
                        castep['free_energy'] = float(line.split('=')[1].split()[0])
                        castep['free_energy_per_atom'] = castep['free_energy'] / castep['num_atoms']
                    elif '0K energy' in line:
                        castep['0K_energy'] = float(line.split('=')[1].split()[0])
                        castep['0K_energy_per_atom'] = castep['0K_energy'] / castep['num_atoms']
                    elif 'Forces' in line:
                        i = 1
                        max_force = 0
                        forces = False
                        while True:
                            if forces:
                                if '*' in final_flines[line_no+i].split()[1]:
                                    forces = False
                                    break
                                else:
                                    force_on_atom = 0
                                    for j in range(3):
                                        temp = final_flines[line_no+i].replace('(cons\'d)', '')
                                        force_on_atom += float(temp.split()[3+j])**2
                                    if force_on_atom > max_force:
                                        max_force = force_on_atom
                            elif 'x' in final_flines[line_no+i]:
                                i += 1                      # skip next blank line
                                forces = True
                            i += 1
                        castep['max_force_on_atom'] = pow(max_force, 0.5)
                    elif 'Stress Tensor' in line:
                        i = 1
                        while i < 20:
                            if 'Cartesian components' in final_flines[line_no+i]:
                                castep['stress'] = []
                                for j in range(3):
                                    castep['stress'].append(list(map(
                                        float, (final_flines[line_no+i+j+4].split()[2:5]))))
                            elif 'Pressure' in final_flines[line_no+i]:
                                try:
                                    castep['pressure'] = float(final_flines[line_no+i].split()[-2])
                                except:
                                    if kwargs.get('dryrun') or verbosity > 0:
                                        print_exc()
                                    pass
                                break
                            i += 1
                    elif 'Atomic Populations (Mulliken)' in line:
                        if castep['spin_polarized']:
                            castep['mulliken_spins'] = []
                            castep['mulliken_net_spin'] = 0.0
                            castep['mulliken_abs_spin'] = 0.0
                        castep['mulliken_charges'] = []
                        castep['mulliken_spins'] = []
                        i = 0
                        while i < len(castep['atom_types']):
                            if castep['spin_polarized']:
                                castep['mulliken_charges'].append(
                                    float(final_flines[line_no+i+4].split()[-2]))
                                castep['mulliken_spins'].append(
                                    float(final_flines[line_no+i+4].split()[-1]))
                                castep['mulliken_net_spin'] += castep['mulliken_spins'][-1]
                                castep['mulliken_abs_spin'] += abs(castep['mulliken_spins'][-1])
                            else:
                                castep['mulliken_charges'].append(
                                    float(final_flines[line_no+i+4].split()[-1]))
                            i += 1
                    elif 'Bond' and 'Population' in line.split():
                        try:
                            castep['bonds'] = []
                            i = 2
                            while True:
                                split = final_flines[line_no+i].split()
                                split[1] = int(split[1])-1
                                split[4] = int(split[4])-1
                                # convert element-wise atom labels to
                                # appropriate index of atom_types
                                for q in range(len(castep['atom_types'])):
                                    if castep['atom_types'][q] == split[0]:
                                        split[1] += q
                                        break
                                for q in range(len(castep['atom_types'])):
                                    if castep['atom_types'][q] == split[3]:
                                        split[4] += q
                                        break
                                castep['bonds'].append(
                                    [[split[1], split[4]], float(split[5]), float(split[6])])
                                i += 1
                                if '===' in final_flines[line_no+i]:
                                    break
                        except:
                            if kwargs.get('dryrun') or verbosity > 0:
                                print_exc()
                            pass
                    elif 'Final Enthalpy' in line:
                        castep['enthalpy'] = float(line.split('=')[-1].split()[0])
                        castep['enthalpy_per_atom'] = (float(line.split('=')[-1].split()[0]) /
                                                       castep['num_atoms'])
                    elif 'Final bulk modulus' in line:
                        try:
                            castep['bulk_modulus'] = float(line.split('=')[-1].split()[0])
                        except:
                            continue
                # calculate kpoint spacing if not found
                if 'kpoints_mp_grid' in castep and 'kpoints_mp_spacing' not in castep and \
                   'lattice_cart' in castep:
                    if db:
                        prec = 2
                    else:
                        prec = 3
                    castep['kpoints_mp_spacing'] = calc_mp_spacing(castep['lattice_cart'],
                                                                   castep['kpoints_mp_grid'],
                                                                   prec)
        # computing metadata, i.e. parallelism, time, memory, version
        for line in flines:
            if 'Release CASTEP version' in line:
                castep['castep_version'] = line.split()[-2]
            elif 'Run started:' in line:
                year = line.split()[5]
                month = str(strptime(line.split()[4], '%b').tm_mon)
                day = line.split()[3]
                castep['date'] = day+'-'+month+'-'+year
            elif 'Total time' in line:
                castep['total_time_hrs'] = float(line.split()[-2])/3600
            elif 'Peak Memory Use' in line:
                castep['peak_mem_MB'] = int(float(line.split()[-2])/1024)
        # check that any optimized results were saved and raise errors if not
        if castep['optimised'] is False:
            raise DFTError('CASTEP GO failed to converge.')
        if 'positions_frac' not in castep:
            raise RuntimeError('Could not find positions')
        if 'enthalpy' not in castep:
            castep['enthalpy'] = 'xxx'
            castep['enthalpy_per_atom'] = 'xxx'
        if 'pressure' not in castep:
            castep['pressure'] = 'xxx'
        if 'cell_volume' not in castep:
            castep['cell_volume'] = 'xxx'
        if 'space_group' not in castep:
            castep['space_group'] = 'xxx'
    except Exception as oops:
        if kwargs.get('dryrun') or verbosity > 0:
            print_exc()
            print('Error in .castep file', seed, 'skipping...')
        if type(oops) == DFTError:
            if db:
                # if importing to db, skip unconverged structure
                # and report in log file
                return seed + '\t\t' + str(type(oops)) + ' ' + str(oops)+'\n', False
            else:
                # if not importing to db, return unconverged structure
                # but notify that it is unconverged
                return castep, False
        else:
            if type(oops) == IOError:
                print_exc()
            return seed + '\t\t' + str(type(oops)) + ' ' + str(oops)+'\n', False
    if kwargs.get('debug'):
        print(json.dumps(castep, indent=2, ensure_ascii=False))
    return castep, True


class DFTError(Exception):
    """ Quick DFT exception class for unconverged or
    non-useful calculations.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
