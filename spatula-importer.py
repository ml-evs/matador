#!/usr/bin/python
# coding: utf-8
from __future__ import print_function
from collections import defaultdict
from os import walk, getcwd
from time import strptime
from sys import argv
import argparse
import pymongo as pm
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

    def __init__(self):
        self.init = True
        self.import_count = 0

    def dict2db(self, struct):
        ''' Insert completed Python dictionary into chosen
        database.
        '''
        client = pm.MongoClient()
        db = client.crystals
        pressure = db.pressure
        struct_id = pressure.insert_one(struct).inserted_id
        print('Inserted', struct_id)
        return 1

    def files2db(self, file_lists, push2db=False):
        ''' Take all files found by scan and appropriately create dicts
        holding all available data; optionally push to database.
        '''
        print('###### RUNNING IMPORTER ######')
        for root in file_lists:
            if root=='.':
                root_str = getcwd().split('/')[-1]
            else:
                root_str = root
            print('Dictifying', root_str, '...')
            airss, cell, param, dir = 4*[False]
            if file_lists[root]['res_count'] > 0:
                if file_lists[root]['castep_count'] < file_lists[root]['res_count']:
                    airss = True
            if airss:
                if file_lists[root]['cell_count'] == 1:
                    cell_dict = self.cell2dict(root + '/' + file_lists[root]['cell'][0])
                    cell = True
                elif file_lists[root]['cell_count'] > 1:
                    print('Multiple cell files found!')
                if file_lists[root]['param_count'] == 1:
                    param_dict = self.param2dict(root + '/' + file_lists[root]['param'][0])
                    param = True
                elif file_lists[root]['param_count'] > 1:
                    print('Multiple param files found!')
                if(file_lists[root]['cell_count'] == 0 or \
                    file_lists[root]['param_count'] == 0):
                    print('Will try to get data from dir names')
                    dir_dict = self.dir2dict(root)
                    if 'source' not in dir_dict:
                        print('No information found in dirname', root)
                    else:
                        dir = True
                # combine cell and param dicts for folder
                if cell and param:
                    input_dict = cell_dict.copy()
                    input_dict.update(param_dict)
                else:
                    if dir:
                        input_dict = dir_dict.copy()
                        if cell:
                            input_dict.update(cell_dict)
                        elif param:
                            input_dict.update(param_dict)
                # create res dicts and combine them with input_dict
                for ind, file in enumerate(file_list[root]['res']):
                    res_dict = self.res2dict(root + '/' + file)
                    final_struct = input_dict.copy()
                    final_struct.update(res_dict)
            else:
                for ind, file in enumerate(file_list[root]['castep']):
                    final_struct = self.castep2dict(root + '/' + file)
                    print(final_struct['max_force_on_atom'], final_struct['optimised'])
                # print(json.dumps(final_struct,indent=2))

        return
     
    def scan_dir(self, topdir):
        ''' Scans folder topdir recursively, returning list of 
        CASTEP/AIRSS input/output files.
        '''
        ResCount, CellCount, CastepCount, ParamCount = 4*[0]
        file_lists = dict()
        if topdir == '.':
            topdir_string = getcwd().split('/')[-1]
        else: 
            topdir_string = topdir
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
                elif file.endswith('.castep'):
                    file_lists[root]['castep'].append(file)
                    file_lists[root]['castep_count'] += 1
                    CastepCount += 1
                elif file.endswith('.cell'):
                    file_lists[root]['cell'].append(file)
                    file_lists[root]['cell_count'] += 1
                    CellCount += 1
                elif file.endswith('.param'):
                    file_lists[root]['param'].append(file)
                    file_lists[root]['param_count'] += 1
                    ParamCount += 1
        print('done!\n')
        prefix = '\t'
        print(prefix, ResCount, '\t.res files')
        print(prefix, CastepCount, '\t.castep files')
        print(prefix, CellCount, '\t.cell files')
        print(prefix, ParamCount, '\t.param files\n')
        return file_lists

    ######################## FILE SCRAPER FUNCTIONS ########################

    def res2dict(self, seed):
        ''' Extract available information from .res file; preferably
        used in conjunction with cell or param file, else dict will 
        include an 'incomplete' tag.
        '''
        # use defaultdict to allow for easy appending
        res = defaultdict(list)
        # read .res file into array
        if seed.endswith('.res'):
            seed = seed.replace('.res', '')
        with open(seed+'.res', 'r') as f:
            flines = f.readlines()
        # add .res to source 
        print(seed + '.res')
        res['source'].append(seed+'.res')
        # alias special lines in res file
        for line in flines:
            if 'TITL' in line:
                titl = line.split()
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
        # need to calculate res['lattice_cart'] and res['positions_cart'] 
        # tag document as incomplete, hopefully to be added to by other routines
        res['tags'].append('incomplete')
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
        # add cell file to source
        cell['source'].append(seed+'.cell')
        for line_no, line in enumerate(flines):
            if line.startswith('#'):
                continue
            elif '%block species_pot' in line.lower():
                cell['species_pot'] = dict()
                i = 1
                while 'endblock' not in flines[line_no+i].lower():
                    cell['species_pot'][flines[line_no+i].split()[0]] = flines[line_no+i].split()[1].split('/')[-1]
                    i += 1
            elif '%block external_pressure' in line.lower():
                cell['external_pressure'] = list()
                cell['external_pressure'].append(map(float, flines[line_no+1].split()))
                cell['external_pressure'].append(map(float, flines[line_no+2].split()))
                cell['external_pressure'].append(map(float, flines[line_no+3].split()))
            elif 'mp_spacing' in line.lower():
                cell['kpoints_mp_spacing'] = line.split()[-1]
            elif 'mp_grid' in line.lower():
                cell['kpoints_mp_grid'] = map(int, line.split()[-3:])
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
        for line in flines:
            line = line.lower()
            if line.startswith('#'):
                continue
            else:
                # exclude some useless info
                scrub_list = ['checkpoint', 'write_bib', 'mix_history_length',
                              'fix_occupancy', 'page_wvfns', 'num_dump_cycles',
                              'backup_interval', 'geom_max_iter', 'fixed_npw',
                              'write_cell_structure', 'bs_write_eigenvalues',
                              'calculate_stress', 'opt_strategy']
                if [rubbish for rubbish in scrub_list if rubbish in line]:
                    continue
                else:
                    if ':' in line:
                        param[line.split(':')[0].strip()] = line.split(':')[-1].strip()
                    elif '=' in line:
                        param[line.split('=')[0].strip()] = line.split('=')[-1].strip()
                    else:
                        param[line.split()[0].strip()] = line.split('')[-1].strip()
        return param

    def dir2dict(self, seed):
        ''' Try to extract information from directory name; last hope
        if no param file has been found. 
        '''
        dir_dict = defaultdict(list)
        if seed == '.':
            seed = getcwd().split('/')[-1]
        dir_as_list = seed.split('-')
        if len(dir_as_list) != 7:
            return dict()
        dir_dict['source'].append(seed)
        dir_dict['energy_cut_off'] = dir_as_list[1]
        dir_dict['kpoints_mp_spacing'] = dir_as_list[2]
        dir_dict['xc_functional'] = dir_as_list[-2]
        dir_dict['task'] = dir_as_list[-1]
        # print(json.dumps(dir_dict,indent=2))
        return dir_dict

    def castep2dict(self, seed):
        ''' From seed filename, create dict of the most relevant
        information about a calculation.
        '''
        # use defaultdict to allow for easy appending
        castep = defaultdict(list)
        # read .castep file
        if seed.endswith('.castep'):
            seed = seed.replace('.castep', '')
        with open(seed+'.castep', 'r') as f:
            flines = f.readlines()
        # set source tag to castep file 
        castep['source'].append(seed+'.castep')
        if 'CollCode' in seed:
            castep['icsd_ref'] = seed.split('CollCode')[-1] 
        # wrangle castep file for basic parameters
        for line_no, line in enumerate(flines):
            if 'task' not in castep and 'type of calculation' in line:
                castep['task'] = line.split(':')[-1].strip()
            elif 'xc_functional' not in castep and 'functional' in line:
                castep['xc_functional'] = line.split(':')[-1].strip()
            elif 'cut_off_energy' not in castep and 'plane wave basis set' in line:
                castep['cut_off_energy'] = float(line.split(':')[-1].split()[0])
            elif 'finite_basis_corr' not in castep and 'finite basis set correction  ' in line:
                castep['finite_basis_corr'] = line.split(':')[-1].strip()
            elif 'kpoints_mp_grid' not in castep and 'MP grid size for SCF' in line:
                castep['kpoints_mp_grid'] = map(int, list(line.split('is')[-1].split()))
            elif 'kpoints_calculated' not in castep and 'Number of kpoints used' in line:
                castep['kpoints_calculated'] = int(line.split('=')[-1])
            elif 'space_group' not in castep and 'Space group of crystal' in line:
                castep['space_group'] = line.split(':')[-1].split(',')[0].strip()
            elif 'external_pressure' not in castep and 'External pressure/stress' in line:
                castep['external_pressure'] = list()
                castep['external_pressure'].append(map(float, flines[line_no+1].split()))
                castep['external_pressure'].append(map(float, flines[line_no+2].split()))
                castep['external_pressure'].append(map(float, flines[line_no+3].split()))
            elif 'atom types' not in castep and 'Cell Contents' in line:
                castep['atom_types'] = list()
                i = 1
                atoms = False
                while True:
                    if atoms:
                        if 'xxxxxxxxx' in flines[line_no+i]:
                            atoms = False
                            break
                        else:
                            castep['atom_types'].append(flines[line_no+i].split()[1])
                    if 'x------' in flines[line_no+i]:
                        atoms= True
                    i += 1
                castep['num_atoms'] = len(castep['atom_types'])
            elif 'Mass of species in AMU' in line:
                i = 1
                atomic_masses = list()
                while True:
                    if len(flines[line_no+i].strip())==0:
                        break
                    else:
                        atomic_masses.append([float(flines[line_no+i].split()[1]),
                                                    flines[line_no+i].split()[0]])
                    i += 1
                atomic_masses.sort(reverse=True)
            elif 'species_pot' not in castep and 'Files used for pseudopotentials' in line:
                castep['species_pot'] = dict()
                i = 1
                while True:
                    if len(flines[line_no+i].strip())==0:
                        break
                    else:
                        castep['species_pot'][flines[line_no+i].split()[0].strip()] = flines[line_no+i].split()[1].strip()
                        i += 1
            # don't check if final_energy exists, as this will update for each GO step
            elif 'Final energy, E' in line:
                castep['total_energy'] = float(line.split('=')[1].split()[0])
                castep['total_energy_per_atom'] = castep['total_energy'] / castep['num_atoms']
            elif 'Final free energy' in line:
                castep['free_energy'] = float(line.split('=')[1].split()[0])
                castep['free_energy_per_atom'] = castep['free_energy'] / castep['num_atoms']
        # write zero pressure if not found in file
        if 'external_pressure' not in castep:
            castep['external_pressure'] = [[0.0, 0.0, 0.0], [0.0, 0.0], [0.0]]
        # generate atomic mass-ordered chemical composition string
        castep['composition'] = ''
        count = list()
        for ind, species in enumerate(atomic_masses):
            count_tmp = 0
            for atom in castep['atom_types']:
                if species[1]==atom:
                    count_tmp += 1
            count.append(count_tmp)
        reducible = True
        for num in count:
            if num % min(count) != 0:
                reducible = False
        min_count = min(count)
        if reducible:
            count = [num/min_count for num in count]
        for ind, species in enumerate(atomic_masses):
            castep['composition'] += species[1]
            if count[ind] != 1: 
                castep['composition'] += str(count[ind])
        # task specific options
        if castep['task'] == 'geometry optimization':
            final = False
            for line_no, line in enumerate(flines):
                if 'Geometry optimization failed to converge' in line:
                    castep['optimised'] = False
                elif 'Final Configuration' in line:
                    final = True
                    if 'optimised' not in castep:
                        castep['optimised'] = True
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
                                        force_on_atom += float(flines[line_no+i].split()[3+j])**2
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
            
        return castep
    

if __name__ == '__main__':
    filename = argv[1]
    importer = Spatula()
    file_list = importer.scan_dir(filename)
    # print(json.dumps(file_list,indent=2))
    importer.files2db(file_list)
    # importer.cell2dict(filename)
    # importer.param2dict(filename)
    # for root in file_list:
        # for castep in file_list[root]['castep']:
            # struct = importer.castep2dict(root + '/' + castep)
            # importer.import_count += importer.dict2db(struct)
    # print('Successfully imported', importer.import_count, 'structures.')
