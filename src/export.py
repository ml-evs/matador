#!/usr/bin/python
# coding: utf-8
""" This file implements functions that
can create a file from a db document.
"""
from __future__ import print_function
import string
import numpy as np
from os.path import exists, isfile, expanduser
from os import system, makedirs
from traceback import print_exc


def query2files(cursor, *args):
    """ Write either .res or .cell + .param files
    for all documents in a cursor.
    """
    args = args[0]
    cell = args.get('cell')
    param = args.get('param')
    res = args.get('res')
    prefix = (args.get('prefix') + '-') if args.get('prefix') is not None else ''
    pressure = args.get('write_pressure')
    if args['subcmd'] == 'polish' or args['subcmd'] == 'swaps':
        info = False
        hash = False
    else:
        info = True
        hash = False
    try:
        num = len(cursor)
    except:
        num = cursor.count()
    print('Writing', num, 'structures to file...')
    try:
        if len(cursor) > 10000:
            write = raw_input('This operation will write ' + str(len(cursor)) + ' structures' +
                              ' are you sure you want to do this? [y/n] ')
            if write == 'y' or write == 'Y':
                print('Writing them all.')
                write = True
            else:
                write = False
                return
        else:
            write = True
    except:
        write = True
        pass
    # handle new file and directory names
    dirname = args['subcmd'] + '-'
    if args['composition'] is not None:
        for comp in args['composition']:
            dirname += comp
    elif args['formula'] is not None:
        dirname += args['formula'][0]
    if args['db'] is not None:
        dirname += '-' + args['db'][0]
    if args.get('swap') is not None:
        for swap in args['swap']:
            dirname += '-' + swap
        if args.get('hull_cutoff') is not None:
            dirname += '-hull-' + str(args.get('hull_cutoff')) + 'eV'
    if args.get('id') is not None:
        dirname += '-' + args.get('id')[0] + '_' + args.get('id')[1]
    dirname = dirname.replace('--', '-')
    dir = False
    dir_counter = 0
    while not dir:
        if dir_counter != 0:
            directory = dirname + str(dir_counter)
        else:
            directory = dirname
        if not exists(directory):
            makedirs(directory)
            dir = True
        else:
            dir_counter += 1
    for ind, doc in enumerate(cursor):
        name = prefix
        path = directory + '/'
        # write either cell, res or both
        for source in doc['source']:
            source = str(source)
            if args['subcmd'] == 'swaps':
                comp_string = ''
                comp_list = []
                for atom in doc['atom_types']:
                    if atom not in comp_list:
                        comp_list.append(atom)
                        comp_string += atom
                name = comp_string + '-' + name
            if '.res' in source:
                name += source.split('/')[-1].split('.')[0]
            elif '.castep' in source:
                name += source.split('/')[-1].split('.')[0]
            elif '.history' in source:
                name += source.split('/')[-1].split('.')[0]
            elif 'OQMD' in source:
                stoich_string = ''
                # prepend new stoich
                if len(doc['stoichiometry']) == 1:
                    stoich_string += doc['stoichiometry'][0][0]
                else:
                    for atom in doc['stoichiometry']:
                        stoich_string += atom[0]
                        stoich_string += str(atom[1]) if atom[1] != 1 else ''
                # grab OQMD entry_id
                name += stoich_string + '-OQMD_' + source.split(' ')[-1]
        path += name
        # always write param for each doc; also handles dirs
        if param:
            doc2param(doc, path, hash_dupe=hash)
        if cell:
            doc2cell(doc, path, pressure, hash_dupe=hash)
        if res:
            doc2res(doc, path, info=info, hash_dupe=hash)


def doc2param(doc, path, hash_dupe=True, *args):
    """ Write basic .param file from single doc.

    doc       : the document to write out to file
    path      : the desired path to file
    hash_dupe : hash duplicate file names, or skip?

    """
    try:
        if path.endswith('.param'):
            path.replace('.param', '')
        param_set = set(['task', 'cut_off_energy', 'xc_functional', 'write_cell_structure',
                         'finite_basis_corr', 'spin_polarized', 'smearing_width',
                         'write_bib', 'finite_basis_corr', 'calculate_stress',
                         'page_wvfns', 'geom_method', 'geom_max_iter', 'write_checkpoint',
                         'fix_occupancy', 'metals_method', 'max_scf_cycles', 'cut_off_energy',
                         'opt_strategy', 'page_wvfns', 'num_dump_cycles', 'bs_write_eigenvalues',
                         'backup_interval', 'fixed_npw', 'mix_cut_off_energy', 'mix_charge_amp',
                         'mixing_scheme', 'mix_charge_gmax', 'geom_force_tol', 'perc_extra_bands',
                         'elec_energy_tol', 'grid_scale', 'spin', 'continuation', 'mix_spin_amp',
                         'spin_treatment', 'spin_fix', 'geom_spin_fix', 'spin_unit',
                         'fine_grid_scale', 'spectral_task'])
        param_dict = dict()
        for param in [param for param in param_set if param in doc]:
            param_dict[param] = doc[param]
        if isfile(path+'.param'):
            if hash_dupe:
                print('File name already exists, generating hash...')
                path += '-' + generate_hash()
            else:
                print('File name already exists! Skipping!')
                raise RuntimeError('Duplicate file!')
        with open(path+'.param', 'w') as f:
            f.write('# Param file generated by matador (Matthew Evans 2016)\n')
            for param in param_dict:
                if param != 'source':
                    f.write("{0:30}: {1}\n".format(param, param_dict[param]))
    except Exception as oops:
        print('Writing param file failed!')
        print(oops)


def doc2cell(doc, path, pressure=None, hash_dupe=True, copy_pspots=True, *args):
    """ Write .cell file for single doc.

    doc         : the document to write to file
    path        : the path to the file
    pressure    : the pressure to write the file at
    hash_dupe   : hash duplicate file names or skip?
    copy_pspots : try to copy pspots from ~/pspots?

    """
    if path.endswith('.cell'):
        path.replace('.cell', '')
    try:
        if isfile(path+'.cell'):
            if hash_dupe:
                print('File name already exists, generating hash...')
                path += '-' + generate_hash()
            else:
                raise RuntimeError('Duplicate file!')
        with open(path+'.cell', 'w') as f:
            f.write('# Cell file generated by matador (Matthew Evans 2016)\n')
            # f.write('# enthalpy_per_atom = ' +
                    # '{: 10f} eV\n\n'.format(doc['enthalpy_per_atom']))
            f.write('\n%BLOCK LATTICE_CART\n')
            for vec in doc['lattice_cart']:
                for coeff in vec:
                    f.write(str(coeff) + ' ')
                f.write('\n')
            f.write('%ENDBLOCK LATTICE_CART\n')
            f.write('\n%BLOCK POSITIONS_FRAC\n')
            if 'atomic_init_spins' in doc:
                for ind, atom in enumerate(zip(doc['atom_types'], doc['positions_frac'])):
                    if atom[0] in doc['atomic_init_spins']:
                        f.write("{0:8s} {1[0]: 15f} {1[1]: 15f} {1[2]: 15f} SPIN={2:}\n".format(
                                atom[0], atom[1], doc['atomic_init_spins'][atom[0]]))
                    else:
                        f.write("{0:8s} {1[0]: 15f} {1[1]: 15f} {1[2]: 15f}\n".format(
                            atom[0], atom[1]))
            else:
                for ind, atom in enumerate(zip(doc['atom_types'], doc['positions_frac'])):
                    f.write("{0:8s} {1[0]: 15f} {1[1]: 15f} {1[2]: 15f}\n".format(
                            atom[0], atom[1]))
            f.write('%ENDBLOCK POSITIONS_FRAC\n\n')
            if pressure is not None:
                f.write('\n%block external_pressure\n'.upper())
                for pressures in pressure:
                    pressures = str(pressures)
                if len(pressure) == 1:
                    f.write(pressure[0] + ' 0 0\n')
                    f.write(pressure[0] + ' 0\n')
                    f.write(pressure[0] + '\n')
                elif len(pressure) == 6:
                    f.write(pressure[0] + ' ' + pressure[1] + ' ' + pressure[2] + '\n')
                    f.write(pressure[3] + ' ' + pressure[4] + '\n')
                    f.write(pressure[5] + '\n')
                f.write('%endblock external_pressure\n'.upper())
            if 'kpoints_mp_spacing' in doc:
                f.write('kpoints_mp_spacing : ' + str(doc['kpoints_mp_spacing']) + '\n')
            elif 'kpoints_mp_grid' in doc:
                f.write('kpoints_mp_grid : ' +
                        str(doc['kpoints_mp_grid'][0]) + ' ' +
                        str(doc['kpoints_mp_grid'][1]) + ' ' +
                        str(doc['kpoints_mp_grid'][2]) + '\n')
            if 'spectral_kpoints_mp_offset' in doc:
                f.write('SPECTRAL_KPOINTS_MP_OFFSET ' +
                        str(doc['spectral_kpoints_mp_offset'][0]) + ' ' +
                        str(doc['spectral_kpoints_mp_offset'][1]) + ' ' +
                        str(doc['spectral_kpoints_mp_offset'][2]) + '\n')
            if 'spectral_kpoints_mp_spacing' in doc:
                f.write('SPECTRAL_KPOINTS_MP_SPACING ' +
                        str(doc['spectral_kpoints_mp_spacing']) + '\n')
            if 'cell_constraints' in doc:
                f.write('\n%BLOCK CELL_CONSTRAINTS\n')
                f.write((''.join(str(doc['cell_constraints'][0]).strip('[]'))+'\n').replace(',', ''))
                f.write((''.join(str(doc['cell_constraints'][1]).strip('[]'))+'\n').replace(',', ''))
                f.write('%ENDBLOCK CELL_CONSTRAINTS\n')
            if 'fix_com' in doc:
                f.write('FIX_COM : TRUE\n')
            if 'symmetry_generate' in doc:
                f.write('SYMMETRY_GENERATE\n')
            if 'snap_to_symmetry' in doc:
                f.write('SNAP_TO_SYMMETRY\n')
            if 'symmetry_tol' in doc:
                f.write('SYMMETRY_TOL : ' + str(doc['symmetry_tol']) + '\n')
            if 'hubbard_u' in doc:
                f.write('\n%BLOCK HUBBARD_U\n')
                f.write('eV\n')
                for elem in doc['hubbard_u']:
                    for orbital in doc['hubbard_u'][elem]:
                        shift = str(doc['hubbard_u'][elem][orbital])
                        f.write(elem + ' ' + orbital + ' : ' + shift + '\n')
                f.write('%ENDBLOCK HUBBARD_U\n')
            if 'quantisation_axis' in doc:
                f.write('\nQUANTISATION_AXIS : ')
                for integer in doc['quantisation_axis']:
                    f.write(str(integer) + ' ')
                f.write('\n\n')
            if 'species_pot' in doc:
                f.write('\n%BLOCK SPECIES_POT\n')
                for elem in doc['species_pot']:
                    if copy_pspots:
                        # copy across pspots if they exist
                        if not isfile(''.join(path.split('/')[:-1])+'/'+doc['species_pot'][elem]):
                            if isfile(expanduser('~/pspot/' + doc['species_pot'][elem])):
                                system('cp ' + expanduser('~/pspot/') + doc['species_pot'][elem] +
                                       ' ' + ''.join(path.split('/')[:-1]))
                    f.write(elem + '\t' + doc['species_pot'][elem] + '\n')
                f.write('%ENDBLOCK SPECIES_POT')
            else:
                # by default, fill file with 00PBE's in a sensible dir but comment out
                f.write('\n#%BLOCK SPECIES_POT\n')
                for elem in doc['stoichiometry']:
                    f.write('#' + elem[0] + ' ~/pspot/' + elem[0] + '_00PBE.usp\n')
                f.write('#%ENDBLOCK SPECIES_POT\n')

    except Exception as oops:
        if hash_dupe:
            print(type(oops), oops)
        else:
            print_exc()


def doc2res(doc, path, info=True, hash_dupe=True, *args):
    """ Write .res file for single doc. """
    if path.endswith('.res'):
        path.replace('.res', '')
    try:
        if isfile(path+'.res'):
            if hash_dupe:
                print('File name already exists, generating hash...')
                path += '-' + generate_hash()
            else:
                raise RuntimeError('Skipping duplicate structure...')
        with open(path+'.res', 'w') as f:
            if not info:
                f.write('TITL file generated by matador (Matthew Evans 2016)')
            else:
                f.write('TITL ')
                f.write(path.split('/')[-1] + ' ')
                if type(doc['pressure']) == str:
                    f.write('0.00 ')
                else:
                    f.write(str(doc['pressure']) + ' ')
                f.write(str(doc['cell_volume']) + ' ')
                f.write(str(doc['enthalpy']) + ' ')
                f.write('0 0 ')             # spin
                f.write(str(doc['num_atoms']) + ' ')
                try:
                    if 'x' in doc['space_group']:
                        f.write('(P1) ')
                    else:
                        f.write('(' + str(doc['space_group']) + ')' + ' ')
                except:
                    f.write('(P1) ')
                f.write('n - 1')
            f.write('\n')
            f.write('CELL ')
            f.write('1.0 ')
            for vec in doc['lattice_abc']:
                for coeff in vec:
                    f.write(' ' + str(coeff))
            f.write('\n')
            f.write('LATT -1\n')
            f.write('SFAC \t')
            written_atoms = []
            for elem in doc['atom_types']:
                if elem not in written_atoms:
                    f.write(' ' + str(elem))
                    written_atoms.append(str(elem))
            f.write('\n')
            atom_labels = []
            i = 0
            j = 1
            while i < len(doc['atom_types']):
                num = doc['atom_types'].count(doc['atom_types'][i])
                atom_labels.extend(num*[j])
                i += num
                j += 1
            for atom in zip(doc['atom_types'], atom_labels, doc['positions_frac']):
                f.write("{0:8s}{1:3d}{2[0]: 15f} {2[1]: 15f} {2[2]: 15f}   1.0\n".format(
                    atom[0], atom[1], atom[2]))
            f.write('END')
            # very important newline for compatibliy with cryan
            f.write('\n')
    except Exception as oops:
        if hash_dupe:
            print('Writing res file failed for ', doc['text_id'])
            print(type(oops), oops)
        else:
            pass


def generate_hash(hashLen=6):
    """ Quick hash generator, based on implementation in PyAIRSS by J. Wynn. """
    hashChars = hashChars = [str(x) for x in range(0, 10)]+[x for x in string.ascii_lowercase]
    hash = ''
    for i in range(hashLen):
        hash += np.random.choice(hashChars)
    return hash
