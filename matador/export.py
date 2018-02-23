#!/usr/bin/python
# coding: utf-8
""" This file implements functions that
can create a file from a db document.
"""
from __future__ import print_function
# matador internals
from matador.utils.cell_utils import cart2abcstar, frac2cart, cart2abc
from matador.utils.cell_utils import abc2cart, calc_mp_grid
from matador.utils.cursor_utils import display_results
# external libraries
import numpy as np
# standard library
import string
from os.path import exists, isfile, expanduser
from os import system, makedirs, remove
from traceback import print_exc
from sys import exit


def query2files(cursor, *args, **kwargs):
    """ Write either .res or .cell + .param files
    for all documents in a cursor.
    """
    args = args[0]
    cell = args.get('cell')
    top = args.get('top')
    param = args.get('param')
    res = args.get('res')
    pdb = args.get('pdb')
    md = args.get('markdown')
    tex = args.get('latex')
    argstr = kwargs.get('argstr')
    multiple_files = cell or param or res or pdb
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
    if top is not None:
        if top < num:
            cursor = cursor[:top]
            num = top
    if multiple_files:
        print('Intending to write', num, 'structures to file...')
        if len(cursor) > 10000:
            try:
                write = input('This operation will write ' + str(len(cursor)) + ' structures' +
                              ' are you sure you want to do this? [y/n] ')
            except:
                print_exc()
                print('Stop using Python2!')
                print('Going to assume your answer was yes...')
                write = 'y'
            if write.lower() is 'y':
                print('Writing them all.')
                write = True
            else:
                write = False
                return
        else:
            write = True
    dirname = generate_relevant_path(args)
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
            if '.res' in source or '.castep' in source or '.history' in source:
                if args['subcmd'] == 'swaps':
                    comp_string = ''
                    comp_list = []
                    for atom in doc['atom_types']:
                        if atom not in comp_list:
                            comp_list.append(atom)
                            comp_string += atom
                    name = comp_string + '-swap-'
                root_fname = source.split('/')[-1].split('.')[0].replace('-swap-', '')
                name += root_fname
            elif 'OQMD' in source:
                stoich_string = ''
                # prepend old stoich
                if len(doc['stoichiometry']) == 1:
                    stoich_string += doc['stoichiometry'][0][0]
                else:
                    for atom in doc['stoichiometry']:
                        stoich_string += atom[0]
                        stoich_string += str(atom[1]) if atom[1] != 1 else ''
                name = stoich_string + '-OQMD_' + source.split(' ')[-1]
                # if swaps, prepend new composition
                if args['subcmd'] == 'swaps':
                    comp_string = ''
                    comp_list = []
                    for atom in doc['atom_types']:
                        if atom not in comp_list:
                            comp_list.append(atom)
                            comp_string += atom
                    name = comp_string + '-' + name
                # grab OQMD entry_id
                if 'icsd' in doc and 'CollCode' not in source:
                    name += '-CollCode' + doc['icsd']
        path += name
        if param:
            doc2param(doc, path, hash_dupe=hash)
        if cell:
            doc2cell(doc, path, pressure, hash_dupe=hash)
        if res:
            doc2res(doc, path, info=info, hash_dupe=hash)
        if pdb:
            doc2pdb(doc, path, hash_dupe=hash)

    if md:
        md_path = path.split('/')[0] + '/' + path.split('/')[0] + '.md'
        print('Writing markdown file', md_path + '...')
        hull = True if args['subcmd'] in ['hull', 'voltage'] else False
        md_string = display_results(cursor, args, argstr=argstr, markdown=True, hull=hull)
        with open(md_path, 'w') as f:
            f.write(md_string)
    if tex:
        tex_path = path.split('/')[0] + '/' + path.split('/')[0] + '.tex'
        print('Writing LaTeX file', tex_path + '...')
        hull = True if args['subcmd'] in ['hull', 'voltage'] else False
        tex_string = display_results(cursor, args, argstr=argstr, latex=True, hull=hull)
        with open(tex_path, 'w') as f:
            f.write(tex_string)

    print('Done!')


def doc2param(doc, path, hash_dupe=True, overwrite=False, spin=False, *args):
    """ Write basic .param file from single doc.

    Input:

        | doc       : the document to write out to file
        | path      : the desired path to file

    Args:

        | hash_dupe : hash duplicate file names, or skip?
        | overwrite : bool, overwrite if file exists.
        | spin      : bool, enforce spin symmetry broken to magic number of 5

    """
    try:
        if path.endswith('.param'):
            path = path.replace('.param', '')
        from matador.utils.castep_params import CASTEP_PARAMS, CASTEP_VERSION
        param_set = set(CASTEP_PARAMS)
        param_dict = dict()
        for param in [param for param in param_set if param in doc]:
            param_dict[param] = doc[param]
        if isfile(path+'.param'):
            if overwrite:
                remove(path + '.param')
            elif hash_dupe:
                print('File name already exists, generating hash...')
                path += '-' + generate_hash()
            else:
                print('File name already exists! Skipping!')
                raise RuntimeError('Duplicate file!')
        with open(path+'.param', 'w') as f:
            f.write('# Param file generated by matador (Matthew Evans 2016) from CASTEP version {} parameter set\n'.format(CASTEP_VERSION))
            for param in param_dict:
                if param != 'source':
                    f.write("{0:30}: {1}\n".format(param, param_dict[param]))

            if doc.get('spin_polarized'):
                total_spin = None
                if 'atomic_init_spins' in doc:
                    total_spin = sum(doc['atomic_init_spins'])
                elif spin:
                    total_spin = 5

                if total_spin is not None:
                    f.write("{0:30}: {1}\n".format('spin', total_spin))

        if 'encapsulated' in doc and doc['encapsulated']:
            try:
                from implicit_cnts import implicit_cnt_params
            except ImportError:
                exit('Failed to import implicit_cnt_params, please ensure pyairss is on your PYTHONPATH!')
            cnt_params = implicit_cnt_params(doc['cnt_radius'])
            flines = []
            flines.append('%BLOCK DEVEL_CODE\n')
            flines.append('ADD_EXT_LOCPOT: \"gaussian_cylinder\"\n')
            flines.append('1D_STRESS_TENSOR: True\n')
            flines.append('gaussian_cylinder_pot:\n')
            flines.append('V0 = {V0}\nradius = {radius}\nbroadening = {fwhm}\naxial = \"0 0 1\"\ncentre = \"0.5 0.5 0\"\n'.format(**cnt_params))
            flines.append(':endgaussian_cylinder_pot\n')
            flines.append('%ENDBLOCK DEVEL_CODE\n')

            with open(path+'.param', 'a') as f:
                for line in flines:
                    f.write(line)

    except Exception as oops:
        print_exc()
        print('Writing param file failed!')


def doc2cell(doc, path, pressure=None, hash_dupe=True, copy_pspots=True, overwrite=False, spin=False, *args):
    """ Write .cell file for single doc.

    Input:

        | doc         : dict, the document to write to file
        | path        : str, the desired path to the file

    Args:

        | pressure    : float, the pressure to write the file at
        | hash_dupe   : bool, hash duplicate file names or skip?
        | copy_pspots : bool, try to copy pspots from ~/pspots?
        | overwrite   : bool, overwrite if file exists
        | spin        : bool, break spin symmetry with magic number 5 on first atom

    """
    if path.endswith('.cell'):
        path = path.replace('.cell', '')
    try:
        if isfile(path+'.cell'):
            if hash_dupe:
                print('File already exists, generating hash...')
                path += '-' + generate_hash()
            elif overwrite:
                remove(path+'.cell')
            else:
                raise RuntimeError('Skipping duplicate structure...')
        with open(path+'.cell', 'w') as f:
            f.write('# Cell file generated by matador (Matthew Evans 2016)\n')
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
                    if ind == 0 and spin:
                        # if spin is True, break spin symmetry on first atom a la run.pl
                        f.write("{0:8s} {1[0]: 15f} {1[1]: 15f} {1[2]: 15f}\tSPIN=5\n".format(atom[0], atom[1]))
                    else:
                        f.write("{0:8s} {1[0]: 15f} {1[1]: 15f} {1[2]: 15f}\n".format(atom[0], atom[1]))
            f.write('%ENDBLOCK POSITIONS_FRAC\n\n')
            if 'external_pressure' in doc:
                f.write('\n%block external_pressure\n'.upper())
                pressure = doc['external_pressure'][0]
                f.write(str(pressure[0]) + ' ' + str(pressure[1]) + ' ' + str(pressure[2]) + '\n')
                pressure = doc['external_pressure'][1]
                f.write(str(pressure[0]) + ' ' + str(pressure[1]) + '\n')
                pressure = doc['external_pressure'][2]
                f.write(str(pressure[0]) + '\n')
                f.write('%endblock external_pressure\n'.upper())
            if 'kpoints_list' in doc:
                f.write('%BLOCK KPOINTS_LIST' + '\n')
                for kpoint in doc['kpoints_list']:
                    f.write('{p[0]:f} {p[1]:f} {p[2]:f} {p[3]:f}\n'.format(p=kpoint))
                f.write('%ENDBLOCK KPOINTS_LIST' + '\n')
            elif 'kpoints_mp_spacing' in doc:
                f.write('kpoints_mp_spacing : ' + str(doc['kpoints_mp_spacing']) + '\n')
            elif 'kpoints_mp_grid' in doc:
                f.write('kpoints_mp_grid : ' +
                        str(doc['kpoints_mp_grid'][0]) + ' ' +
                        str(doc['kpoints_mp_grid'][1]) + ' ' +
                        str(doc['kpoints_mp_grid'][2]) + '\n')
            if 'kpoints_mp_offset' in doc:
                f.write('kpoints_mp_offset : ' +
                        str(doc['kpoints_mp_offset'][0]) + ' ' +
                        str(doc['kpoints_mp_offset'][1]) + ' ' +
                        str(doc['kpoints_mp_offset'][2]) + '\n')
            if 'spectral_kpoints_mp_offset' in doc:
                f.write('SPECTRAL_KPOINTS_MP_OFFSET ' +
                        str(doc['spectral_kpoints_mp_offset'][0]) + ' ' +
                        str(doc['spectral_kpoints_mp_offset'][1]) + ' ' +
                        str(doc['spectral_kpoints_mp_offset'][2]) + '\n')
            if 'spectral_kpoints_mp_spacing' in doc:
                f.write('SPECTRAL_KPOINTS_MP_SPACING ' +
                        str(doc['spectral_kpoints_mp_spacing']) + '\n')
            elif 'spectral_kpoints_mp_grid' in doc:
                f.write('SPECTRAL_KPOINTS_MP_GRID ' +
                        str(doc['spectral_kpoints_mp_grid'][0]) + ' ' +
                        str(doc['spectral_kpoints_mp_grid'][1]) + ' ' +
                        str(doc['spectral_kpoints_mp_grid'][2]) + '\n')
            elif 'spectral_kpoints_list' in doc:
                f.write('%BLOCK SPECTRAL_KPOINTS_LIST\n')
                for point in doc['spectral_kpoints_list']:
                    f.write('{p[0]} {p[1]} {p[2]}\n'.format(p=point))
                f.write('%ENDBLOCK SPECTRAL_KPOINTS_LIST\n')
            if 'phonon_kpoint_list' in doc:
                f.write('%BLOCK PHONON_KPOINT_LIST\n')
                for point in doc['phonon_kpoint_list']:
                    f.write('{p[0]} {p[1]} {p[2]}\n'.format(p=point))
                f.write('%ENDBLOCK PHONON_KPOINT_LIST\n')
            if 'phonon_fine_kpoint_list' in doc:
                f.write('%BLOCK PHONON_FINE_KPOINT_LIST\n')
                for point in doc['phonon_fine_kpoint_list']:
                    f.write('{p[0]} {p[1]} {p[2]}\n'.format(p=point))
                f.write('%ENDBLOCK PHONON_FINE_KPOINT_LIST\n')
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
            if 'positions_noise' in doc:
                f.write('\nPOSITIONS_NOISE : {}\n'.format(doc['positions_noise']))
            if 'cell_noise' in doc:
                f.write('\nCELL_NOISE : {}\n'.format(doc['cell_noise']))
            if 'species_pot' in doc:
                f.write('\n%BLOCK SPECIES_POT\n')
                for elem in doc['species_pot']:
                    if copy_pspots:
                        # copy across pspots if they exist
                        if not isfile(''.join(path.split('/')[:-1])+'/'+doc['species_pot'][elem]):
                            if isfile(expanduser('~/pspot/' + doc['species_pot'][elem])):
                                system('cp ' + expanduser('~/pspot/') + doc['species_pot'][elem] +
                                       ' ' + ''.join(path.split('/')[:-1]))
                    if elem == 'library':
                        f.write(doc['species_pot']['library'] + '\n')
                    else:
                        f.write(elem + '\t' + doc['species_pot'][elem] + '\n')
                f.write('%ENDBLOCK SPECIES_POT')
            else:
                # by default, fill file with 00PBE's in a sensible dir but comment out
                f.write('\n#%BLOCK SPECIES_POT\n')
                for elem in set(doc['atom_types']):
                    f.write('#' + elem[0] + ' ~/pspot/' + elem[0] + '_00PBE.usp\n')
                f.write('#%ENDBLOCK SPECIES_POT\n')

    except Exception:
        print_exc()
        print('Continuing...')


def doc2pdb(doc, path, info=True, hash_dupe=True, *args):
    """ Write a simple .pdb for single doc. """
    if path.endswith('.pdb'):
        path = path.replace('.pdb', '')
    try:
        if isfile(path+'.pdb'):
            if hash_dupe:
                print('File already exists, generating hash...')
                path += '-' + generate_hash()
            else:
                raise RuntimeError('Skipping duplicate structure...')
        with open(path+'.pdb', 'w') as f:
            try:
                HEADER = 'HEADER    {} {}'.format(doc['text_id'][0], doc['text_id'][1])
            except:
                HEADER = 'HEADER    Generated with matador.'
            try:
                # write res file header if info
                TITLE = 'TITLE     '
                TITLE += path.split('/')[-1] + ' '
                if type(doc['pressure']) == str:
                    TITLE += '0.00 '
                else:
                    TITLE += str(doc['pressure']) + ' '
                TITLE += str(doc['cell_volume']) + ' '
                TITLE += str(doc['enthalpy']) + ' '
                TITLE += '0 0 '             # spin
                TITLE += str(doc['num_atoms']) + ' '
                try:
                    if 'x' in doc['space_group']:
                        TITLE += '(P1) '
                    else:
                        TITLE += '(' + str(doc['space_group']) + ')' + ' '
                except:
                    TITLE += '(P1) '
                TITLE += 'n - 1'
            except:
                if not info:
                    TITLE = 'TITLE\t' + path.split('/')[-1]
                raise RuntimeError('Failed to get info for res file, turn info off.')
            AUTHOR = 'AUTHOR    Generated with matador (Matthew Evans, 2016)'
            f.write(HEADER + '\n')
            f.write(TITLE + '\n')
            f.write(AUTHOR + '\n')
            # use dummy SG for CRYST1, shouldn't matter
            CRYST1 = 'CRYST1 {v[0][0]:9.3f} {v[0][1]:9.3f} {v[0][2]:9.3f} {v[1][0]:7.2f} {v[1][1]:7.2f} {v[1][2]:7.2f} P 1'.format(v=doc['lattice_abc'])
            f.write(CRYST1 + '\n')
            SCALEn = cart2abcstar(doc['lattice_cart'])
            f.write('SCALE1    {v[0][0]:10.6f} {v[0][1]:10.6f} {v[0][2]:10.6f}      {:10.5f}\n'.format(0.0, v=SCALEn))
            f.write('SCALE2    {v[1][0]:10.6f} {v[1][1]:10.6f} {v[1][2]:10.6f}      {:10.5f}\n'.format(0.0, v=SCALEn))
            f.write('SCALE3    {v[2][0]:10.6f} {v[2][1]:10.6f} {v[2][2]:10.6f}      {:10.5f}\n'.format(0.0, v=SCALEn))
            if 'positions_abs' not in doc:
                doc['positions_abs'] = frac2cart(doc['lattice_cart'], doc['positions_frac'])
            for ind, atom in enumerate(doc['atom_types']):
                try:
                    HETATM = 'HETATM '
                    # append 00 to atom type, a la cell2pdb...
                    HETATM += '{:4d} {:.4} NON A   1     '.format(ind+1, atom+'00')
                    HETATM += '{v[0]:7.3f} {v[1]:7.3f} {v[2]:7.3f} {:5.2f} {:5.2f}          {:.2}'.format(1.0, 0.0, atom, v=doc['positions_abs'][ind])
                    f.write(HETATM + '\n')
                except:
                    print_exc()
            TER = 'TER       {}       NON A   1'.format(len(doc['atom_types']))
            f.write(TER + '\n')
            f.write('END')
    except:
        if hash_dupe:
            print_exc()
            print('Writing pdb file failed for ', doc['text_id'])
        else:
            print_exc()
            pass


def doc2pwscf(doc, path, template=None, spacing=None, *args):
    """ Write the structural part of QE input file based
    on the provided matador doc. Will calculate the correct
    kpoint_mp_grid if spacing is provided.

    Input:

        | doc  : dict, matador document,
        | path : str, filename to write out to.

    Args:

        | template : str, filename of template to prepend before structure
                    (with prefix keyword to be replaced),
        | spacing  : float, kpoint_mp_spacing to use when calculating grid.

    """
    if path.endswith('.in'):
        path = path.replace('.in', '')

    if isfile(path + '.in'):
        print('File already exists, not overwriting...')
        return

    if 'lattice_cart' not in doc:
        doc['lattice_cart'] = abc2cart(doc['lattice_abc'])

    if 'kpoints_mp_spacing' in doc or spacing is not None:
        if 'kpoints_mp_spacing' in doc:
            spacing = doc['kpoints_mp_spacing']
        doc['kpoints_mp_grid'] = calc_mp_grid(doc['lattice_cart'], spacing)

    print(doc['kpoints_mp_grid'])

    file_string = ''
    file_string += 'CELL_PARAMETERS angstrom\n'
    for i in range(3):
        file_string += ' {d[0]: 10.10f} {d[1]: 10.10f} {d[2]: 10.10f}\n'.format(d=doc['lattice_cart'][i])

    file_string += '\n ATOMIC_POSITIONS crystal\n'
    for i in range(len(doc['atom_types'])):
        file_string += '{:4} {d[0]: 10.10f} {d[1]: 10.10f} {d[2]: 10.10f}\n'.format(doc['atom_types'][i], d=doc['positions_frac'][i])
    file_string += '\nK_POINTS automatic\n'
    file_string += '{d[0]} {d[1]} {d[2]} 0 0 0'.format(d=doc['kpoints_mp_grid'])

    if template is not None:
        if isfile(template):
            with open(template, 'r') as f:
                template_string = f.readlines()

    with open(path + '.in', 'w') as f:
        for line in template_string:
            if 'prefix' in line:
                line = '  prefix = \'{}\',\n'.format(path)
            elif 'nat' in line:
                line = '  nat = {},\n'.format(len(doc['atom_types']))
            elif 'nat' in line:
                line = '  ntyp = {},\n'.format(len(set(doc['atom_types'])))
            f.write(line)
        f.write(file_string)


def doc2res(doc, path, info=True, hash_dupe=True, spoof_titl=False, overwrite=False, sort_atoms=True, *args):
    """ Write .res file for single doc.

    Input:

        | doc  : dict, matador document containing structure
        | path : str, desired filename for res file

    Args:

        | info       : bool, print info in res file header
        | hash_dupe  : bool, add random hash to colliding filenames
        | spoof_titl : bool, make up fake info for file header (for use with e.g. cryan)
        | overwrite  : bool, overwrite files with conflicting filenames
        | sorted     : bool, if False, atoms are not sorted (this will not be a valid res file)

    """
    if path.endswith('.res'):
        path = path.replace('.res', '')
    if spoof_titl:
        info = False
    try:
        if isfile(path+'.res'):
            if hash_dupe:
                print('File already exists, generating hash...')
                path += '-' + generate_hash()
            elif overwrite:
                remove(path+'.res')
            else:
                raise RuntimeError('Skipping duplicate structure...')
        if not info:
            space_group = 'P1' if 'space_group' not in doc else doc['space_group']
            if spoof_titl:
                titl = ('TITL {} -1 1 -1 0 0 {} ({}) n - 1').format(path.split('/')[-1], doc['num_atoms'], space_group)
            else:
                titl = ('TITL file generated by matador (Matthew Evans 2016)')
        else:
            try:
                titl = 'TITL '
                titl += (path.split('/')[-1] + ' ')
                if 'pressure' not in doc or type(doc['pressure']) == str:
                    titl += '0.00 '
                else:
                    titl += str(doc['pressure']) + ' '
                titl += str(doc['cell_volume']) + ' '
                titl += str(doc['enthalpy']) + ' '
                titl += '0 0 '             # spin
                titl += str(doc['num_atoms']) + ' '
                if 'space_group' not in doc:
                    titl += '(P1) '
                elif 'x' in doc['space_group']:
                    titl += '(P1) '
                else:
                    titl += '(' + str(doc['space_group']) + ')' + ' '
                titl += 'n - 1'
            except:
                raise RuntimeError('Failed to get info for res file, turn info off.')
        if 'encapsulated' in doc and doc['encapsulated']:
            rem = ("REM NTPROPS {{\'chiralN\': {}, \'chiralM\': {}, \'r\': {}, "
                   "\'offset\': [0.5, 0.5, 0.5], \'date\': \'xxx\', \'eformperfu\': 12345, \'z\': {}}}\n"
                   .format(doc['cnt_chiral'][0], doc['cnt_chiral'][1], doc['cnt_radius'], doc['cnt_length']))
        flines = []
        if 'encapsulated' in doc and doc['encapsulated']:
            flines.append(rem)
        flines.append(titl)
        flines.append('\n')
        flines.append('CELL ')
        flines.append('1.0 ')
        if 'lattice_abc' not in doc or len(doc['lattice_abc']) != 2 or len(doc['lattice_abc'][0]) != 3 or len(doc['lattice_abc'][1]) != 3:
            try:
                doc['lattice_abc'] = cart2abc(doc['lattice_cart'])
            except:
                raise RuntimeError('Failed to get lattice, something has gone wrong for {}'.format(path))
        for vec in doc['lattice_abc']:
            for coeff in vec:
                flines.append(' ' + str(round(coeff, 8)))
        flines.append('\n')
        flines.append('LATT -1\n')
        flines.append('SFAC \t')

        # enforce correct order by elements, sorting only the atom_types, not the positions inside them
        if sort_atoms:
            positions_frac, atom_types = zip(*[(pos, types) for (types, pos) in
                                             sorted(zip(doc['atom_types'], doc['positions_frac']),
                                                    key=lambda k: k[0])])
        else:
            positions_frac = doc['positions_frac']
            atom_types = doc['atom_types']

        written_atoms = []
        for elem in atom_types:
            if elem not in written_atoms:
                flines.append(' ' + str(elem))
                written_atoms.append(str(elem))
        flines.append('\n')
        atom_labels = []
        i = 0
        j = 1
        while i < len(atom_types):
            num = atom_types.count(atom_types[i])
            atom_labels.extend(num*[j])
            i += num
            j += 1
        for atom in zip(atom_types, atom_labels, positions_frac):
            flines.append("{0:8s}{1:3d}{2[0]: 15f} {2[1]: 15f} {2[2]: 15f}   1.0\n".format(
                atom[0], atom[1], atom[2]))
        flines.append('END')
        # very important newline for compatibliy with cryan
        flines.append('\n')

        # actually write to file
        with open(path+'.res', 'w') as f:
            for line in flines:
                f.write(line)

    except Exception:
        print_exc()
        if hash_dupe:
            print('Writing res file failed for ', path)
        else:
            print('Writing res file failed for ', path, 'this is not necessarily a problem...')
            pass


def generate_hash(hashLen=6):
    """ Quick hash generator, based on implementation in PyAIRSS by J. Wynn. """
    hashChars = hashChars = [str(x) for x in range(0, 10)]+[x for x in string.ascii_lowercase]
    hash = ''
    for i in range(hashLen):
        hash += np.random.choice(hashChars)
    return hash


def generate_relevant_path(args):
    """ Generates a suitable path name based on query. """
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
    return dirname
