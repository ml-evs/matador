#!/usr/bin/python
# coding: utf-8
from __future__ import print_function
from collections import defaultdict
from os import walk, getcwd, stat
from pwd import getpwuid
from time import strptime
from sys import argv
from fractions import gcd
from math import pi, log10
import gzip
import argparse
import pymongo as pm
import random
import json

class Spatula:
    ''' The Spatula class implements methods to scrape folders 
    and individual files for crystal structures and create a 
    MongoDB document for each. 

    Files types that can be read are:
    
        * CASTEP output
        * CASTEP .param, .cell input
        * airss.pl / pyAIRSS .res output
    '''

    def __init__(self, dryrun=False, debug=False, verbosity=0, tags=None, scratch=False):
        ''' Set up arguments and initialise DB client. '''
        self.init = True
        self.import_count = 0
        # I/O files 
        self.logfile = open('spatula.log', 'w')
        try:
            # wordfile = open('/home/matthew/crysdb-bacon/src/new_words', 'r')
            # nounfile = open('/home/matthew/crysdb-bacon/src/nouns', 'r')
            wordfile = open('/u/fs1/me388/crysdb-bacon/src/new_words', 'r')
            nounfile = open('/u/fs1/me388/crysdb-bacon/src/nouns', 'r')
        except Exception as oopsy:
            exit(oopsy)
        self.wlines = wordfile.readlines()
        self.num_words = len(self.wlines)
        self.nlines = nounfile.readlines()
        self.num_nouns = len(self.nlines)
        wordfile.close()
        nounfile.close()
        self.dryrun = dryrun
        self.debug = debug 
        self.verbosity = verbosity 
        self.scratch = scratch
        self.tag_dict = dict()
        self.tag_dict['tags'] = tags
        if not self.dryrun:
            self.client = pm.MongoClient()
            self.db = self.client.crystals
            if self.scratch:
                self.repo = self.db.scratch
            else:
                self.repo = self.db.repo
        # scan directory on init
        self.file_lists = self.scan_dir()
        # convert to dict and db if required
        self.files2db(self.file_lists)
        if not self.dryrun:
            print('Successfully imported', self.import_count, 'structures!')
            # index by enthalpy for faster/larger queries
            indexed = False
            for index in self.repo.list_indexes():
                if 'enthalpy_per_atom' in index['name']:
                    print('Index found, rebuilding...')
                    self.repo.reindex()
                    indexed = True
            if indexed == False:
                print('Building enthalpy index...')
                self.repo.create_index([('enthalpy_per_atom', pm.ASCENDING)])
        else:
            print('Dryrun complete!')
        self.logfile.close()
        self.logfile = open('spatula.log', 'r')
        errors = sum(1 for line in self.logfile)
        if errors == 1:
            print('There is', errors, 'error to view in spatala.log')
        elif errors == 0:
            print('There were no errors.')
        elif errors > 1:
            print('There are', errors, 'errors to view in spatala.log')
        self.logfile.close()

    def dict2db(self, struct):
        ''' Insert completed Python dictionary into chosen
        database, with generated text_id.
        '''
        plain_text_id = [self.wlines[random.randint(0,self.num_words-1)].strip(),
                         self.nlines[random.randint(0,self.num_nouns-1)].strip()]
        struct['text_id'] = plain_text_id
        if 'tags' in self.tag_dict:
            struct['tags'] = self.tag_dict['tags']
        struct_id = self.repo.insert_one(struct).inserted_id
        if self.debug:
            print('Inserted', struct_id)
        return 1

    def files2db(self, file_lists):
        ''' Take all files found by scan and appropriately create dicts
        holding all available data; optionally push to database.
        '''
        print('\n###### RUNNING IMPORTER ######\n')
        multi = False
        for root in file_lists:
            if root=='.':
                root_str = getcwd().split('/')[-1]
            else:
                root_str = root
            if self.verbosity > 0:
                print('Dictifying', root_str, '...')
            airss, cell, param, dir = 4*[False]
            if file_lists[root]['res_count'] > 0:
                if file_lists[root]['castep_count'] < file_lists[root]['res_count']:
                    if file_lists[root]['cell_count'] <= file_lists[root]['res_count']:
                        airss = True
            if airss:
                if file_lists[root]['param_count'] == 1:
                    param_dict = self.param2dict(root + '/' + file_lists[root]['param'][0])
                    param = True
                elif file_lists[root]['param_count'] > 1:
                    if self.dryrun:
                        print('Multiple param files found!')
                    multi = True
                if file_lists[root]['cell_count'] == 1:
                    cell_dict = self.cell2dict(root + '/' + file_lists[root]['cell'][0])
                    cell = True
                elif file_lists[root]['cell_count'] > 1:
                    multi = True
                    if self.dryrun:
                        print('Multiple cell files found - searching for param file with same name...')
                if multi:
                    for param_name in file_lists[root]['param']:
                        for cell_name in file_lists[root]['cell']:
                            if param_name.split('.')[0] in cell_name:
                                param_dict = self.param2dict(root + '/' + param_name)
                                cell_dict = self.cell2dict(root + '/' + cell_name)
                                cell = True 
                                param = True
                                if self.verbosity > 0:
                                    print('Found matching cell and param files:', param_name)
                                break
                # always try to scrape dir 
                dir_dict = self.dir2dict(root)
                if 'source' in dir_dict:
                    dir = True
                # combine cell and param dicts for folder
                input_dict = dict()
                if dir:
                    input_dict = dir_dict.copy()
                if cell and param:
                    input_dict.update(cell_dict)
                    input_dict.update(param_dict)
                    input_dict['source'] = cell_dict['source'] + param_dict['source']
                    if dir:
                        input_dict['source'] = input_dict['source'] + dir_dict['source']
                else:
                    if dir:
                        input_dict = dir_dict.copy()
                        if cell:
                            input_dict.update(cell_dict)
                            input_dict['source'] = cell_dict['source'] + dir_dict['source']
                        elif param:
                            input_dict.update(param_dict)
                            input_dict['source'] = param_dict['source'] + dir_dict['source']
                # create res dicts and combine them with input_dict
                for ind, file in enumerate(file_lists[root]['res']):
                    res_dict = self.res2dict(root + '/' + file)
                    if res_dict != False:
                        final_struct = input_dict.copy()
                        final_struct.update(res_dict)
                        # calculate kpoint spacing if not found
                        recip_abc = 3*[0]
                        for j in range(3):
                            recip_abc[j] = 2 * pi / float(final_struct['lattice_abc'][0][j])
                            if 'kpoints_mp_spacing' not in final_struct:
                                if 'kpoints_mp_grid' in final_struct:
                                    max_spacing = 0
                                    for j in range(3):
                                        spacing = recip_abc[j]/(2 * pi * final_struct['kpoints_mp_grid'][j])
                                        max_spacing = spacing if spacing > max_spacing else max_spacing
                                        final_struct['kpoints_mp_spacing'] = round(max_spacing + 0.5*10**(round(log10(max_spacing)-1)), 2)
                        try:
                            final_struct['source'] = res_dict['source'] + input_dict['source']
                        except:
                            pass
                        if not self.dryrun:
                            final_struct.update(self.tag_dict)
                            self.import_count += self.dict2db(final_struct)
            else:
                for ind, file in enumerate(file_lists[root]['castep']):
                    castep_dict = self.castep2dict(root + '/' + file)
                    if castep_dict != False:
                        final_struct = castep_dict
                        if not self.dryrun:
                            final_struct.update(self.tag_dict)
                            self.import_count += self.dict2db(final_struct)
        return
     
    def scan_dir(self):
        ''' Scans folder topdir recursively, returning list of 
        CASTEP/AIRSS input/output files.
        '''
        ResCount, CellCount, CastepCount, ParamCount = 4*[0]
        file_lists = dict()
        topdir = '.'
        topdir_string = getcwd().split('/')[-1]
        print('Scanning', topdir_string, 'for CASTEP/AIRSS output files... ',
              end='')
        for root, dirs, files in walk(topdir, followlinks=True, topdown=True):
            file_lists[root] = defaultdict(list)
            file_lists[root]['res_count'] = 0
            file_lists[root]['cell_count'] = 0
            file_lists[root]['param_count'] = 0
            file_lists[root]['castep_count'] = 0
            for file in files:
                if file.endswith('.res'):
                    file_lists[root]['res'].append(file)
                    file_lists[root]['res_count'] += 1
                    ResCount += 1
                elif file.endswith('.castep') or file.endswith('.history') or file.endswith('.history.gz'):
                    file_lists[root]['castep'].append(file)
                    file_lists[root]['castep_count'] += 1
                    CastepCount += 1
                elif file.endswith('.cell'):
                    if file.endswith('-out.cell'):
                        continue
                    else:
                        file_lists[root]['cell'].append(file)
                        file_lists[root]['cell_count'] += 1
                        CellCount += 1
                elif file.endswith('.param'):
                    file_lists[root]['param'].append(file)
                    file_lists[root]['param_count'] += 1
                    ParamCount += 1
        print('done!\n')
        prefix = '\t\t'
        print(prefix, "{:d}".format(ResCount), '\t\t.res files')
        print(prefix, CastepCount, '\t\t.castep, .history or .history.gz files')
        print(prefix, CellCount, '\t\t.cell files')
        print(prefix, ParamCount, '\t\t.param files\n')
        return file_lists

    ######################## FILE SCRAPER FUNCTIONS ########################

    def res2dict(self, seed):
        ''' Extract available information from .res file; preferably
        used in conjunction with cell or param file.
        '''
        # use defaultdict to allow for easy appending
        res = defaultdict(list)
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
            if self.debug:
                print(seed+'.res has no owner.')
            res['user'] == 'xxx'
        if 'CollCode' in seed:
            res['icsd'] = seed.split('CollCode')[-1] 
        # alias special lines in res file
        try:
            for line in flines:
                if 'TITL' in line:
                    titl = line.split()
                    if len(titl) != 12:
                        raise RuntimeError('res file header missing some data')
                elif 'CELL' in line:
                    cell = line.split()
            res['pressure'] = float(titl[2])
            res['cell_volume'] = float(titl[3])
            res['enthalpy'] = float(titl[4])
            res['num_atoms'] = int(titl[7])
            res['space_group'] = titl[8].strip('()')
            res['enthalpy_per_atom'] = res['enthalpy'] / res['num_atoms']
            res['total_energy'] = res['enthalpy'] - res['pressure']*res['cell_volume']
            res['total_energy_per_atom'] = res['total_energy'] / res['num_atoms']
            res['lattice_abc'] = [map(float, cell[2:5]), map(float, cell[5:8])]
            for line_no, line in enumerate(flines):
                if 'SFAC' in line:
                    i = 1
                    while 'END' not in flines[line_no+i]:
                        cursor = flines[line_no+i].split()
                        res['atom_types'].append(cursor[0])
                        res['positions_frac'].append(map(float, cursor[2:5]))
                        i += 1
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
        except Exception as oopsy:
            if self.verbosity > 0:
                print(oopsy)
                print('Error in .res file', seed+ '.res, skipping...')
            self.logfile.write(seed+'.res\t\t'+str(oopsy)+'\n')
            return False
        if self.verbosity > 4:
            print(json.dumps(res, indent=2))
        return res

    def cell2dict(self, seed):
        ''' Extract available information from .cell file; probably
        to be merged with another dict from a .param or .res file.
        '''
        cell = defaultdict(list)
        if seed.endswith('.cell'):
            seed = seed.replace('.cell', '')
        with open(seed+'.cell', 'r') as f:
            flines = f.readlines()
        try:
            # add cell file to source
            cell['source'].append(seed+'.cell')
            for line_no, line in enumerate(flines):
                if line.startswith('#'):
                    continue
                elif '%block species_pot' in line.lower():
                    cell['species_pot'] = dict()
                    i = 1
                    while 'endblock' not in flines[line_no+i].lower():
                        try:
                            cell['species_pot'][flines[line_no+i].split()[0]] = flines[line_no+i].split()[1].split('/')[-1]
                        except:
                            pass
                        i += 1
                elif '%block external_pressure' in line.lower():
                    cell['external_pressure'] = list()
                    for j in range(3):
                       flines[line_no+j+1] = flines[line_no+j+1].replace(',', '')
                       cell['external_pressure'].append(map(float, flines[line_no+j+1].split()))
                elif 'mp_spacing' in line.lower():
                    cell['kpoints_mp_spacing'] = line.split()[-1]
                elif 'mp_grid' in line.lower():
                    cell['kpoints_mp_grid'] = map(int, line.split()[-3:])
            if 'external_pressure' not in cell:
                cell['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0], [0.0]]
        except Exception as oopsy:
            if self.verbosity > 0:
                print(oopsy)
                print('Error in', seed +'.cell, skipping...')
            self.logfile.write(seed + '\t\t' + str(oopsy))
            return False
        if self.debug:
            print(json.dumps(cell,indent=2))
        return cell

    def param2dict(self, seed):
        ''' Extract available information from .param file; probably
        to be merged with other dicts from other files.
        '''
        param = defaultdict(list)
        if seed.endswith('.param'):
            seed = seed.replace('.param', '')
        with open(seed+'.param', 'r') as f:
            flines = f.readlines()
        param['source'].append(seed+'.param')
        try:
            for line_no, line in enumerate(flines):
                line = line.lower()
                # skip blank lines and comments
                if line.startswith(('#', '!')) or len(line.strip())==0:
                    continue
                else:
                    # exclude some useless info
                    scrub_list = ['checkpoint', 'write_bib', 'mix_history_length',
                                  'fix_occupancy', 'page_wvfns', 'num_dump_cycles',
                                  'backup_interval', 'geom_max_iter', 'fixed_npw',
                                  'write_cell_structure', 'bs_write_eigenvalues',
                                  'calculate_stress', 'opt_strategy', 'max_scf_cycles']
                    false_str = ['False', 'false', '0']
                    splitters = [':', '=', ' ']
                    if [rubbish for rubbish in scrub_list if rubbish in line]:
                        continue
                    # read all other parameters in
                    else:
                        if [splitter for splitter in splitters if splitter in line]: 
                            param[line.split(splitter)[0].strip()] = line.split(':')[-1].strip()
                            if 'spin_polarized' in line:
                                if [false for false in false_str if false in param['spin_polarized']]:
                                    param['spin_polarized'] = False
                                else:
                                    param['spin_polarized'] = True
                            if 'cut_off_energy' in line:
                                param['cut_off_energy'] = int(param['cut_off_energy'].replace('ev','').strip())
                            if 'xc_functional' in line:
                                param['xc_functional'] = param['xc_functional'].upper() 
        except Exception as oopsy:
            if self.verbosity > 0:
                print(oopsy)
                print('Error in', seed+'.param, skipping...')
            self.logfile.write(seed + '\t\t' + str(oopsy))
            return False
        if self.debug:
            print(json.dumps(param,indent=2))
        return param

    def dir2dict(self, seed):
        ''' Try to extract information from directory name; last hope
        if no param file has been found. 
        '''
        dir_dict = defaultdict(list)
        info = False
        if seed == '.':
            seed = getcwd().split('/')[-1]
        dirs_as_list = seed.split('/')
        task_list = ['GO', 'NMR', 'OPTICS'] 
        phase_list = ['alpha', 'beta', 'gamma', 'theta']
        defect_list = ['vacancy', 'interstitial']
        try:
            for dir in dirs_as_list:
                if dir=='.':
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
                    else:
                        dir_dict['species_pot'] = dir.split('-')[offset+5]
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
                if self.verbosity > 0:
                    print('No information found in dirname', seed)
        except Exception as oopsy:
            if self.verbosity > 0:
                print(oopsy)
                print('Error wrangling dir name', seed)
            self.logfile.write(seed + '\t\t' + str(oopsy))
        if self.debug:
            print(json.dumps(dir_dict,indent=2))
        return dir_dict

    def castep2dict(self, seed):
        ''' From seed filename, create dict of the most relevant
        information about a calculation.
        '''
        # use defaultdict to allow for easy appending
        castep = defaultdict(list)
        # read .castep, .history or .history.gz file
        if '.gz' in seed:
            with gzip.open(seed, 'r') as f:
                flines = f.readlines()
        else:
            with open(seed, 'r') as f:
                flines = f.readlines()

        # set source tag to castep file 
        castep['source'].append(seed)
        # grab file owner
        castep['user'] = getpwuid(stat(seed).st_uid).pw_name
        if 'CollCode' in seed:
            castep['icsd'] = seed.split('CollCode')[-1].replace('.castep', '').replace('.history', '').replace('.gz', '')
        try:
            # wrangle castep file for basic parameters
            for line_no, line in enumerate(flines):
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
                    castep['kpoints_mp_grid'] = map(int, list(line.split('is')[-1].split()))
                elif 'sedc_apply' not in castep and 'DFT+D: Semi-empirical dispersion correction    : on' in line:
                    castep['sedc_apply'] = True
                    castep['sedc_scheme'] = flines[line_no+1].split(':')[1].split()[0] 
                elif 'kpoints_calculated' not in castep and 'Number of kpoints used' in line:
                    castep['kpoints_calculated'] = int(line.split('=')[-1])
                elif 'space_group' not in castep and 'Space group of crystal' in line:
                    castep['space_group'] = line.split(':')[-1].split(',')[0].strip().replace(" ", "")
                elif 'external_pressure' not in castep and 'External pressure/stress' in line:
                    castep['external_pressure'] = list()
                    castep['external_pressure'].append(map(float, flines[line_no+1].split()))
                    castep['external_pressure'].append(map(float, flines[line_no+2].split()))
                    castep['external_pressure'].append(map(float, flines[line_no+3].split())) 
                elif 'spin_polarized' not in castep and 'Treating system as spin-polarized' in line:
                    castep['spin_polarized'] = True
                elif 'atom types' not in castep and 'Cell Contents' in line: 
                    castep['atom_types'] = list() 
                    castep['stoichiometry'] = defaultdict(float)
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
                                castep['positions_frac'].append((map(float, (flines[line_no+i].split()[3:6]))))
                        if 'x------' in flines[line_no+i]:
                            atoms= True
                        i += 1
                    castep['num_atoms'] = len(castep['atom_types'])
                    # 
                    min = 2000
                    for atom in castep['atom_types']:
                        if atom not in castep['stoichiometry']:
                            castep['stoichiometry'][atom] = 0
                        castep['stoichiometry'][atom] += 1
                    for atom in castep['atom_types']:
                        if castep['stoichiometry'][atom] < min:
                            min = castep['stoichiometry'][atom]
                    # convert stoichiometry to tuple for fryan
                    temp_stoich = []
                    for key, value in castep['stoichiometry'].iteritems():
                        temp_stoich.append([key, value/min])
                    castep['stoichiometry'] = temp_stoich
                elif 'species_pot' not in castep and 'Files used for pseudopotentials' in line:
                    castep['species_pot'] = dict()
                    i = 1
                    while True:
                        if len(flines[line_no+i].strip())==0:
                            break
                        else:
                            castep['species_pot'][flines[line_no+i].split()[0].strip()] = flines[line_no+i].split()[1].split('/')[-1]
                            i += 1
                # don't check if final_energy exists, as this will update for each GO step
                elif 'Final energy, E' in line:
                    castep['total_energy'] = float(line.split('=')[1].split()[0])
                    castep['total_energy_per_atom'] = castep['total_energy'] / castep['num_atoms']
                elif 'Final free energy' in line:
                    castep['free_energy'] = float(line.split('=')[1].split()[0])
                    castep['free_energy_per_atom'] = castep['free_energy'] / castep['num_atoms']
                elif 'Lattice parameters' in line:
                    castep['lattice_abc'] = list()
                    i = 1
                    castep['lattice_abc'].append(map(float, [flines[line_no+i].split('=')[1].strip().split(' ')[0], 
                                                             flines[line_no+i+1].split('=')[1].strip().split(' ')[0],
                                                             flines[line_no+i+2].split('=')[1].strip().split(' ')[0]]))
                    castep['lattice_abc'].append(map(float, [flines[line_no+i].split('=')[-1].strip(),
                                                             flines[line_no+i+1].split('=')[-1].strip(),
                                                             flines[line_no+i+2].split('=')[-1].strip()]))
                    recip_abc = 3*[0]
                    for j in range(3):
                        recip_abc[j] = 2 * pi / float(castep['lattice_abc'][0][j])
                    if 'kpoints_mp_grid' in castep:
                        max_spacing = 0
                        for j in range(3):
                            spacing = recip_abc[j]/(2 * pi * castep['kpoints_mp_grid'][j])
                            max_spacing = spacing if spacing > max_spacing else max_spacing
                        castep['kpoints_mp_spacing'] = round(max_spacing + 0.5*10**(round(log10(max_spacing)-1)), 2)
                
                elif 'Current cell volume' in line:
                    castep['cell_volume'] = float(line.split('=')[1].split()[0].strip())
            # write zero pressure if not found in file
            if 'external_pressure' not in castep:
                castep['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0], [0.0]]
            # task specific options
            if castep['task'] != 'geometryoptimization' and castep['task'] != 'geometry optimization':
                raise RuntimeError('CASTEP file does not contain GO calculation', castep['task'])
            else:
                final = False
                castep['optimised'] = False
                for line_no, line in enumerate(flines):
                    if 'Geometry optimization completed successfully' in line:
                        castep['optimised'] = True
                        final = True
                    if final:
                        if 'Real Lattice' in line:
                            castep['lattice_cart'] = list()
                            i = 1
                            while True:
                                if len(flines[line_no+i].strip()) == 0:
                                    break
                                else:
                                    castep['lattice_cart'].append((map(float, (flines[line_no+i].split()[0:3]))))
                                i += 1
                        elif 'Lattice parameters' in line:
                            castep['lattice_abc'] = list()
                            i = 1
                            castep['lattice_abc'].append(map(float, [flines[line_no+i].split('=')[1].strip().split(' ')[0], 
                                                                     flines[line_no+i+1].split('=')[1].strip().split(' ')[0],
                                                                     flines[line_no+i+2].split('=')[1].strip().split(' ')[0]]))
                            castep['lattice_abc'].append(map(float, [flines[line_no+i].split('=')[-1].strip(),
                                                                     flines[line_no+i+1].split('=')[-1].strip(),
                                                                     flines[line_no+i+2].split('=')[-1].strip()]))
                            recip_abc = 3*[0]
                            for j in range(3):
                                recip_abc[j] = 2 * pi / float(castep['lattice_abc'][0][j])
                            if 'kpoints_mp_grid' in castep:
                                max_spacing = 0
                                for j in range(3):
                                    spacing = recip_abc[j]/(2 * pi * castep['kpoints_mp_grid'][j])
                                    max_spacing = spacing if spacing > max_spacing else max_spacing
                                castep['kpoints_mp_spacing'] = round(max_spacing + 0.5*10**(round(log10(max_spacing)-1)), 2)
                        
                        elif 'Current cell volume' in line:
                            castep['cell_volume'] = float(line.split('=')[1].split()[0].strip())
                        elif 'Cell Contents' in line :
                            castep['positions_frac'] = list()
                            i = 1
                            atoms = False
                            while True:
                                if atoms:
                                    if 'xxxxxxxxx' in flines[line_no+i]:
                                        atoms = False
                                        break
                                    else:
                                        castep['positions_frac'].append((map(float, (flines[line_no+i].split()[3:6]))))
                                if 'x------' in flines[line_no+i]:
                                    atoms = True
                                i += 1
                        elif 'Forces' in line:
                            i = 1
                            max_force = 0
                            forces = False
                            while True:
                                if forces:
                                    if '*' in flines[line_no+i].split()[1]:
                                        forces = False
                                        break
                                    else:
                                        force_on_atom = 0
                                        for j in range(3):
                                            force_on_atom += float(flines[line_no+i].replace('(cons\'d)', '').split()[3+j])**2
                                        if force_on_atom > max_force:
                                            max_force = force_on_atom
                                elif 'x' in flines[line_no+i]:
                                    i += 1                      # skip next blank line
                                    forces = True
                                i += 1
                            castep['max_force_on_atom'] = pow(max_force, 0.5)
                        elif 'Stress Tensor' in line:
                            i = 1
                            while True and i < 20:
                                if 'Pressure' in flines[line_no+i]:
                                    castep['pressure'] = float(flines[line_no+i].split()[-2])
                                i += 1
                        elif 'Final Enthalpy' in line:
                            castep['enthalpy'] = float(line.split('=')[-1].split()[0])
                            castep['enthalpy_per_atom'] = float(line.split('=')[-1].split()[0])/castep['num_atoms']
                        elif 'Final bulk modulus' in line:
                            try:
                                castep['bulk_modulus'] = float(line.split('=')[-1].split()[0])
                            except:
                                continue 
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
                    castep['peak_mem_MB'] = int(float(line.split()[-2])/1000)
            # check that any optimized results were saved and raise errors if not
            if castep['optimised'] == False:
                raise RuntimeError('CASTEP GO failed to converge.') 
            if 'enthalpy' not in castep:
                raise RuntimeError('Could not find enthalpy')
            if 'positions_frac' not in castep:
                raise RuntimeError('Could not find positions')
            if 'pressure' not in castep:
                castep['pressure'] = 'xxx'
            if 'cell_volume' not in castep:
                castep['cell_volume'] = 'xxx'
            if 'space_group' not in castep:
                castep['space_group'] = 'xxx'
        except Exception as oopsy:
            if self.dryrun:
                print(oopsy)
                print('Error in .castep file', seed, 'skipping...')
            self.logfile.write(seed+ '\t\t' + str(oopsy)+'\n')
            return False
        if self.debug:
            print(json.dumps(castep,indent=2))
        return castep


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Import CASTEP/AIRSS results into MongoDB database.',
            epilog='Written by Matthew Evans (2016)')
    parser.add_argument('-d', '--dryrun', action='store_true',
                        help='run the importer without connecting to the database')
    parser.add_argument('-v', '--verbosity', action='count',
                        help='enable verbose output')
    parser.add_argument('-t', '--tags', nargs='+', type=str,
                        help='set user tags, e.g. nanotube, project name')
    parser.add_argument('--debug', action='store_true',
                        help='enable debug output to print every dict')
    parser.add_argument('-s', '--scratch', action='store_true',
                        help='import to junk collection called scratch')
    args = parser.parse_args()
    importer = Spatula(dryrun=args.dryrun, 
                       debug=args.debug,
                       verbosity=args.verbosity,
                       tags=args.tags,
                       scratch=args.scratch)
