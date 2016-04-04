#!/usr/local/bin/python
# coding: utf-8
from __future__ import print_function
from collections import defaultdict

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

    def dict2db(self):
        ''' Insert completed Python dictionary into chosen
        database.
        '''

    def files2db(self, file_lists, push2db=False):
        ''' Take all files found by scan and appropriately create dicts
        holding all available data; optionally push to database.
        '''
        for root in file_lists:
            if airss:
                cell = cell2dict(root)
                param = param2dict(root)
                if cell and param:
                    input = dict_combine(cell, param)
                elif not cell and not param:
                    input = cell if not param or param if not cell
                if not cell or param:
                    res2dict(root, input, dir)
                res_list = res2dict(root, cell, param)
                dir = dir2dict(root)
                res_list = res2dict(root, dir)
                for res in res_list:
                    dict2db(res)
            else:
                for files in root:
                    dict2db(castep2dict(root, tags))

        return
     
    def scan_dir(self, topdir):
        ''' Scans folder topdir recursively, returning list of 
        CASTEP/AIRSS input/output files.
        '''
        from os import walk, getcwd
        ResCount, CellCount, CastepCount, ParamCount, RootCount = 5*[0]
        file_lists = dict()
        if topdir == '.':
            topdir_string = getcwd()
        else: 
            topdir_string = topdir
        print('Scanning', topdir_string, 'for CASTEP/AIRSS output files... ',
              end='')
        for root, dirs, files in walk(topdir, followlinks=True):
            file_lists[root] = defaultdict(list)
            file_lists[root]['res_count'] = 0
            file_lists[root]['cell_count'] = 0
            file_lists[root]['param_count'] = 0
            file_lists[root]['castep_count'] = 0
            RootCount += 1
            for file in files:
                if '.res' in file:
                    file_lists[root]['res'].append(file)
                    file_lists[root]['res_count'] += 1
                    ResCount += 1
                elif '.castep' in file:
                    file_lists[root]['castep'].append(file)
                    file_lists[root]['castep_count'] += 1
                    CastepCount += 1
                elif '.cell' in file:
                    file_lists[root]['cell'].append(file)
                    file_lists[root]['cell_count'] += 1
                    CellCount += 1
                elif '.param' in file:
                    file_lists[root]['param'].append(file)
                    file_lists[root]['param_count'] += 1
                    ParamCount += 1
        print('done!\n')
        prefix = '\t'
        print(prefix, ResCount, '\t.res files')
        print(prefix, CastepCount, '\t.castep files')
        print(prefix, CellCount, '\t.cell files')
        print(prefix, ParamCount, '\t.param files\n')
        print(getcwd().split('/')[-1])
        print(' │')
        for ind, root in enumerate(file_lists):
            if (file_lists[root]['res_count'] == 0 and \
                    file_lists[root]['castep_count'] == 0):
                continue
            elif file_lists[root]['res_count'] > 0: 
                if (file_lists[root]['cell_count'] == 1 and \
                      file_lists[root]['castep_count'] < file_lists[root]['res_count']):
                    print(' ├──', root[2:], '\n │\t\t└────> AIRSS results.', end='')
                elif (file_lists[root]['param_count'] == 0 or \
                        file_lists[root]['cell_count'] == 0):
                    print(' ├──', root[2:], '\n │\t\t└────> AIRSS results missing',
                            '.cell or .param file.', end='')
            if(ind!=RootCount):
                print('\n │')
        self.files2db(file_lists)
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
        with open(seed+'.res', 'r') as f:
            flines = f.readlines()
        # add .res to source 
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
        return cell

    def param2dict(self, seed):
        ''' Extract available information from .param file; probably
        to be merged with other dicts from other files.
        '''
        param = defaultdict(lsit)
        return param

    def dir2dict(self, seed):
        ''' Try to extract information from directory name; last hope
        if no param file has been found. 
        '''
        dir_dict = defaultdict(list)
        return dir_dict

    def castep2dict(self, seed):
        ''' From seed filename, create dict of the most relevant
        information about a calculation.
        '''
        # use defaultdict to allow for easy appending
        castep = defaultdict(list)
        # read .castep file
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
                        castep['lattice_abc'].append(map(float, [flines[line_no+i:line_no+i+2].split('=')[1].strip().split(' ')[0], 
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
                                    for j in range(1):
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
                        castep['bulk_modulus'] = float(line.split('=')[-1].split()[0])
    
        # computing metadata, i.e. parallelism, time, memory, version
        for line in flines:
            if 'Release CASTEP version' in line:
                castep['castep_version'] = line.split()[-2]
            elif 'Run started:' in line:
                from time import strptime
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
    from sys import argv
    import argparse
    filename = argv[1]
    importer = Spatula()
    importer.scan_dir(filename)
    # importer.res2dict(filename)

    # castep2dict(filename)
