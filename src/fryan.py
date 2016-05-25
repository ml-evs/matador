#!/usr/bin/python
# coding: utf-8
from __future__ import print_function
# import related crysdb functionality
from spatula import param2dict
# import external libraries
import pymongo as pm
import numpy as np
import argparse
import string
import sys
from sys import argv
from os import makedirs, system, uname
from os.path import exists, isfile, expanduser
from copy import deepcopy
from bson.son import SON
import bson.json_util as json
import re

class DBQuery:
    """ Class that implements queries to MongoDB
    structure database.
    """
    
    def __init__(self, **kwargs):
        """ Initialise the query with command line
        arguments and return results.
        """
        # read args
        self.args = kwargs
        for arg in self.args:
            if type(self.args[arg]) == str:
                self.args[arg] = self.args[arg].split() 
        # connect to MongoDB
        local = uname()[1]
        if local == 'cluster2':
            remote = 'node1'
        else:
            remote = None
        self.client = pm.MongoClient(remote)
        self.db = self.client.crystals
        if self.args.get('scratch'):
            self.repo = self.client.crystals.scratch
        else:
            self.repo = self.client.crystals.repo
        # print last spatula report
        self.report = self.client.crystals.spatula
        self.print_report()
        if self.args.get('summary'):
            self.top = self.db.command('collstats', self.repo.name)['count']
        else:
            self.top = self.args.get('top') if self.args.get('top') != None else 10
        # grab all args as string for file dumps
        if self.args.get('sysargs'):
            self.sysargs = ''
            for ind, arg in enumerate(self.args.get('sysargs')):
                if '--write_pressure' in arg:
                    ind += 1
                    self.sysargs += self.args.get('sysargs')[ind] + 'GPa'
                else:
                    self.sysargs += arg
                    self.sysargs += '-'
            self.sysargs = 'query-' + self.sysargs
            self.sysargs = self.sysargs.replace('---','-')
            self.sysargs = self.sysargs.replace('--','-')
        """ PERFORM QUERY """
        self.cursor = EmptyCursor()
        # initalize query_dict to '$and' all queries
        self.query_dict = dict()
        self.query_dict['$and'] = []
        # benchmark enthalpy to display (set by calc_match)
        self.gs_enthalpy = 0.0
        if self.args.get('dbstats'):
            self.dbstats()
            exit()
        if self.args.get('id') != None:
            self.cursor = self.repo.find({'text_id': self.args.get('id')})
            if self.cursor.count() < 1:
                exit('Could not find a match.')
            if self.args.get('calc_match'):
                # save special copy of calc_dict for hulls
                self.calc_dict = dict()
                self.calc_dict['$and'] = []
                # to avoid deep recursion, and since this is always called first
                # don't append, just set
                self.query_dict['$and'] = self.query_calc(self.cursor)
                self.calc_dict['$and'] = list(self.query_dict['$and'])
        if self.args.get('stoichiometry') != None:
            self.query_dict['$and'].append(self.query_stoichiometry())
        if self.args.get('pressure') != None:
            self.query_dict['$and'].append(self.query_pressure())
        if self.args.get('composition') != None:
            self.query_dict['$and'].append(self.query_composition())
        if self.args.get('tags') != None:
            self.query_dict['$and'].append(self.query_tags())
        if self.args.get('user') != None:
            self.query_dict['$and'].append(self.query_user())
        if (len(self.query_dict['$and']) == 0 and
            self.args.get('id') == None and
            not self.args.get('dbstats')):
            self.cursor = self.repo.find().sort('enthalpy_per_atom', pm.ASCENDING)
        # clone cursor for further use after printing
        if self.cursor.count() < 1 or self.args.get('calc_match'):
            if self.args.get('details'):
                print(self.query_dict)
            self.cursor = self.repo.find(SON(self.query_dict)).sort('enthalpy_per_atom', pm.ASCENDING)
        cursor = self.cursor.clone()
        # self.cursor = cursor.clone()
        """ QUERY POST-PROCESSING """
        # write query to res or cell with param files
        if self.args.get('cell') or self.args.get('res'):
            if self.cursor.count() >= 1:
                if self.args.get('top') != None:
                    self.query2files(cursor[:self.top], self.args.get('res'), self.args.get('cell'), top=True, pressure=self.args.get('write_pressure'))
                else:
                    self.query2files(cursor, self.args.get('res'), self.args.get('cell'), pressure=self.args.get('write_pressure'))
        # if called as script, always print results
        print(self.cursor.count(), 'results found for query.')
        if self.args.get('main'):
            if self.cursor.count() >= 1:
                if self.cursor.count() > self.top:
                    self.display_results(cursor[:self.top], details=self.args.get('details'))
                else:
                    self.display_results(cursor, details=self.args.get('details'))

    def __del__(self):
        """ Clean up any temporary databases on garbage 
        collection of DBQuery object.
        """
        try:
            self.temp.drop()
        except:
            pass

    def swaps(self, doc, pairs=1, template_param=None):
        """ Take a db document as input and perform atomic swaps. """
        for source in doc['source']:
            if '.castep' or '.res' in source:
                name = source.split('/')[-1].split('.')[0]
        name = name + '-' + str(pairs) +'-pair-swaps/' + name
        swapDoc = deepcopy(doc)
        swapAtoms = swapDoc['atom_types']
        for i in range(pairs):
            valid = False
            while not valid:    
                swap = np.random.randint(0, len(swapAtoms)-1, size=2)
                if swap[0] != swap[1] and swapAtoms[swap[0]] != swapAtoms[swap[1]]:
                        valid = True
            swapAtoms[swap[1]], swapAtoms[swap[0]] = swapAtoms[swap[0]], swapAtoms[swap[1]]
        swapPos = np.asarray(swapDoc['positions_frac'])
        for i in range(len(swapAtoms)):
           swapPos[i] += np.random.rand(3) * (0.1 / 7.9)
        swapDoc['positions_frac'] = swapPos
        hash = self.generate_hash(8)
        self.doc2cell(swapDoc, name+'-'+hash)
        self.doc2param(swapDoc, name+'-'+hash, template_param)

    def generate_hash(self, hashLen=6):
        """ Quick hash generator, based on implementation in 
        PyAIRSS by J. Wynn. 
        """
        hashChars = hashChars = [str(x) for x in range(0,10)]+[x for x in string.ascii_lowercase]
        hash = ''
        for i in range(hashLen):
            hash += np.random.choice(hashChars)
        return hash

    def query2files(self, cursor, res=False, cell=False, top=False, pressure=None):
        """ Write .res or .cell files for all docs in query,
        including a .param file for each. Eventually handle
        swaps from one element to another from CLI.
        """
        if cursor.count() > 1000 and top == False:
            write = raw_input('This operation will write ' + str(cursor.count()) + ' structures,' \
                    + ' are you sure you want to do this? [y/n] ')
            if write == 'y' or write == 'Y':
                print('Writing them all.')
                write = True
            else:
                write = False
                return
        else:
            write = True
        name = self.sysargs
        dir = False
        dir_counter = 0
        while not dir:
            if dir_counter != 0:
                directory = name + str(dir_counter)
            else:
                directory = name
            if not exists(directory):
                makedirs(directory)
                dir = True
            else:
                dir_counter += 1
        for ind, doc in enumerate(cursor):
            path = directory + '/'
            # write either cell, res or both
            for source in doc['source']:
                if '.res' in source:
                    name = source.split('/')[-1].split('.')[0]
                elif '.castep' in source:
                    name = source.split('/')[-1].split('.')[0]
                elif '.history' in source:
                    name = source.split('/')[-1].split('.')[0]
            path += name
            # always write param for each doc; also handles dirs
            self.doc2param(doc, path)
            if cell:
                self.doc2cell(doc, path, pressure)
            if res:
                self.doc2res(doc, path)

    def doc2param(self, doc, path, template=None):
        """ Write basic .param file from single doc. """
        paramList = ['task', 'cut_off_energy', 'xc_functional',
                'finite_basis_corr', 'spin_polarized']
        seedDict = dict()
        paramDict = dict()
        for param in [param for param in paramList if param in doc]:
            seedDict[param] = doc[param]
        if template != None:
            try:
                paramDict, success = param2dict(template)
                if not success:
                    raise RuntimeError('Failed to open template.')
            except Exception as oops:
                print(oops)
        try:
            if isfile(path+'.param'):
                print('File name already exists, generating hash...')
                path += '-' + self.generate_hash()
        except Exception as oops:
            print('Writing param file failed for ', doc['text_id'])
            print(oops)
        try:
            with open(path+'.param', 'w') as f:
                f.write('# Param file generated by fryan.py (Matthew Evans 2016)\n')
                for param in paramDict:
                    if param != 'source' and param != 'spin_polarized':
                        f.write("{0:20}: {1}\n".format(param, paramDict[param]))
        except Exception as oops:
            if not exists(''.join(path.split('/')[:-1])):
                newDir = ''.join(path.split('/')[:-1])
                makedirs(newDir)
                self.doc2param(doc, path, template=template)

    def doc2cell(self, doc, path, pressure=None):
        """ Write .cell file for single doc. """
        try:
            if isfile(path+'.cell'):
                print('File name already exists, generating hash...')
                path += '-' + self.generate_hash()
            with open(path+'.cell', 'w') as f:
                f.write('# Cell file generated by fryan.py (Matthew Evans 2016)\n\n')
                f.write('%BLOCK LATTICE_ABC\n')
                for vec in doc['lattice_abc']:
                    for coeff in vec:
                        f.write(str(coeff) + ' ')
                    f.write('\n')
                f.write('%ENDBLOCK LATTICE_ABC\n\n')
                f.write('%BLOCK POSITIONS_FRAC\n')
                for ind, atom in enumerate(zip(doc['atom_types'], doc['positions_frac'])):
                    f.write("{0:8s} {1[0]: 15f} {1[1]: 15f} {1[2]: 15f}   1.0\n".format(atom[0], \
                            atom[1]))
                f.write('%ENDBLOCK POSITIONS_FRAC\n\n')
                if pressure != None:
                    f.write('%block external_pressure\n'.upper())
                    for pressures in pressure:
                        pressures = str(pressures)
                    if len(pressure) == 1:
                        f.write(pressure[0] + ' 0 0\n')
                        f.write(pressure[0] + ' 0\n')
                        f.write(pressure[0] + '\n')
                    elif len(pressure)==6:
                        f.write(pressure[0] + ' ' + pressure[1] + ' ' + pressure[2] + '\n')
                        f.write(pressure[3] + ' ' + pressure[4] + '\n')
                        f.write(pressure[5] + '\n')
                    f.write('%endblock external_pressure\n'.upper())
                if 'kpoints_mp_spacing' in doc:
                    f.write('kpoints_mp_spacing : ' + str(doc['kpoints_mp_spacing']) + '\n')
                elif 'kpoints_mp_grid' in doc:
                    f.write('kpoints_mp_grid : ' + str(doc['kpoints_mp_grid'][0]) + ' ' + str(doc['kpoints_mp_grid'][1]) 
                            + ' ' + str(doc['kpoints_mp_grid'][2]) + '\n')
                f.write('\n%BLOCK SPECIES_POT\n')
                for elem in doc['species_pot']:
                    if not isfile(''.join(path.split('/')[:-1]) + '/' + doc['species_pot'][elem]):
                            if isfile(expanduser('~/pspot/' + doc['species_pot'][elem])):
                                system('cp ' + expanduser('~/pspot/') + doc['species_pot'][elem] + \
                                        ' ' + ''.join(path.split('/')[:-1]))
                    f.write(elem + '\t' + doc['species_pot'][elem] + '\n')
                f.write('%ENDBLOCK SPECIES_POT')
        except Exception as oops:
            print('Writing cell file failed for ', doc['text_id'])
            print(oops)
    
    def doc2res(self, doc, path):
        """ Write .res file for single doc. """
        try:
            if isfile(path+'.res'):
                print('File name already exists, generating hash...')
                path += '-' + self.generate_hash()
            with open(path+'.res', 'w') as f:
                # f.write('# Res file generated by fryan.py (Matthew Evans 2016)\n\n')
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
                for ind, atom in enumerate(zip(doc['atom_types'], atom_labels, doc['positions_frac'])):
                    f.write("{0:8s}{1:3d}{2[0]: 15f} {2[1]: 15f} {2[2]: 15f}   1.0\n".format(atom[0], \
                            atom[1], atom[2]))
                f.write('END')
                # very important newline for compatibliy with cryan
                f.write('\n')
        except Exception as oops:
            print('Writing res file failed for ', doc['text_id'])
            print(oops)

    def display_results(self, cursor, details=False):
        """ Print query results in a cryan-like fashion. """
        struct_string = []
        detail_string = []
        detail_substring = []
        source_string = []
        summary_string = []
        formula_string = []
        last_formula = ''
        gs_enthalpy = 0
        header_string = "{:^24}".format('ID')
        header_string += "{:^5}".format('!?!')
        header_string += "{:^12}".format('Pressure')
        header_string += "{:^12}".format('Volume/fu') 
        header_string += "{:^18}".format('Enthalpy/atom')
        header_string += "{:^12}".format('Space group')
        header_string += "{:^10}".format('Formula')
        header_string += "{:^8}".format('# fu')
        for ind, doc in enumerate(cursor):
            formula_substring = ''
            if 'phase' in doc:
                if 'alpha' in doc['phase']:
                    formula_substring += 'α-'
                elif 'beta' in doc['phase']:
                    formula_substring += 'β-'
                elif 'gamma' in doc['phase']:
                    formula_substring += 'γ-'
                elif 'theta' in doc['phase']:
                    formula_substring += 'θ-'
            atom_per_fu = 0
            for item in doc['stoichiometry']:
                for item_ind, subitem in enumerate(item):
                    if item_ind == 0:
                        formula_substring += str(subitem)
                    if item_ind == 1:
                        if subitem != 1:
                            formula_substring += str(subitem)
                        atom_per_fu += subitem
            if last_formula != formula_substring:
                self.gs_enthalpy = 0.0
            formula_string.append(formula_substring)
            struct_string.append(
                    "{:^24}".format(doc['text_id'][0]+' '+doc['text_id'][1]))
            try:
                if doc['quality'] == 0:
                    struct_string[-1] += "{:^5}".format('!!!')
                else:
                    struct_string[-1] += "{:^5}".format((5-doc['quality'])*'?')
            except:
                struct_string[-1] += "{:5}".format(' ')
            try:
                struct_string[-1] += "{:^ 12.3f}".format(doc['pressure'])
            except: 
                struct_string[-1] += "{:^12}".format(doc['pressure'])
            try:
                struct_string[-1] += "{:^12.3f}".format(atom_per_fu * doc['cell_volume'] / doc['num_atoms'])
            except:
                struct_string[-1] += "{:^12}".format('xxx')
            try:
                struct_string[-1] += "{:^18.5f}".format(doc['enthalpy_per_atom'] - self.gs_enthalpy)
            except:
                struct_string[-1] += "{:^18}".format('xxx')
            try:
                struct_string[-1] += "{:^12}".format(doc['space_group'])
            except:
                struct_string[-1] += "{:^12}".format('xxx')
            struct_string[-1] += "{:^10}".format(formula_substring)
            struct_string[-1] += "{:^8}".format(doc['num_atoms']/atom_per_fu)
            if last_formula != formula_substring:
                self.gs_enthalpy = doc['enthalpy_per_atom']
            last_formula = formula_substring
            if details:
                detail_string.append(11 * ' ' + u"├╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ ")
                if self.args.get('source'):
                    detail_substring.append(11 * ' ' + u"├╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ ")
                else:
                    detail_substring.append(11 * ' ' + u"└╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ ")
                if 'spin_polarized' in doc:
                    if doc['spin_polarized']:
                        detail_string[-1] += 'S-'
                if 'sedc_scheme' in doc:
                    detail_string[-1] += doc['sedc_scheme'].upper()+'+'
                if 'xc_functional' in doc:
                    detail_string[-1] += doc['xc_functional']
                else:
                    detail_string[-1] += 'functional unknown for' + doc['source'][0]
                if 'cut_off_energy' in doc:
                    detail_string[-1] += ', ' + "{:4.2f}".format(doc['cut_off_energy']) + ' eV'
                else:
                    detail_string[-1] += 'cutoff unknown'
                if 'external_pressure' in doc:
                    detail_string[-1] += ', ' + "{:4.2f}".format(doc['external_pressure'][0][0]) + ' GPa'
                if 'kpoints_mp_spacing' in doc:
                    detail_string[-1] += ', ' + str(doc['kpoints_mp_spacing']) + ' 1/A'
                if 'species_pot' in doc:
                    try:
                        for species in doc['species_pot']:
                            detail_substring[-1] += doc['species_pot'][species] + ', '
                    except:
                        pass
                if 'icsd' in doc:
                    detail_substring[-1] += 'ICSD-CollCode' + doc['icsd'] + ', '
                if 'tags' in doc:
                    try:
                        for tag in doc['tags']:
                            detail_substring[-1] += tag + ', '
                    except:
                        pass
                if 'user' in doc:
                    detail_substring[-1] += doc['user']
                detail_string[-1] += ' ' + (len(header_string)-len(detail_string[-1])-1)*u"╌"
                detail_substring[-1] += ' ' + (len(header_string)-len(detail_substring[-1])-1)*u"╌"
            if self.args.get('source'):
                source_string.append(11*' ' + u"└───────────────┬──")
                for num, file in enumerate(doc['source']):
                    if num == len(doc['source'])-1:
                        source_string[-1] += (len(u"└────────────── ")+11)*' ' + u'└──'
                    elif num != 0:
                        source_string[-1] += (len(u"└────────────── ")+11)*' ' + u'├──'
                    # elif num == 0:
                    source_string[-1] += ' ' + file[2:] 
                    if num != len(doc['source'])-1:
                        source_string[-1] += '\n'
        print(len(header_string)*'─')
        print(header_string)
        print(len(header_string)*'─')
        if self.args.get('summary'):
            current_formula = ''
            count = 0
            for ind, string in enumerate(formula_string):
                if count > self.top:
                    break
                if string != current_formula and string not in formula_string[:ind]:
                    count += 1
                    print(struct_string[ind])
                    if details:
                        print(detail_string[ind])
                        print(detail_substring[ind])
                    if self.args.get('source'):
                        print(source_string[ind])
                    current_formula = string
        else:
            for ind, string in enumerate(struct_string):
                print(string)
                if details:
                    print(detail_string[ind])
                    print(detail_substring[ind])
                if self.args.get('source'):
                    print(source_string[ind])
                if details or self.args.get('source'):
                    print(len(header_string)*'─')
        
    def query_stoichiometry(self):
        """ Query DB for particular stoichiometry. """
        # alias stoichiometry
        stoich = self.args.get('stoichiometry')
        # if there's only one string, try split it by caps
        if len(stoich) == 1:
            stoich = [elem for elem in re.split(r'([A-Z][a-z]*)', stoich[0]) if elem]
        elements = []
        fraction = []
        for i in range(0, len(stoich), 1):
            if not bool(re.search(r'\d', stoich[i])):
                elements.append(stoich[i])
                try:
                    fraction.append(float(stoich[i+1]))
                except:
                    fraction.append(1.0)
        fraction = np.asarray(fraction)
        fraction /= np.min(fraction)
        query_dict = dict()
        query_dict['$and'] = []
        for ind, elem in enumerate(elements):
            stoich_dict = dict()
            stoich_dict['stoichiometry'] = dict()
            stoich_dict['stoichiometry']['$in'] = [[elem, fraction[ind]]]
            query_dict['$and'].append(stoich_dict)
        if not self.args.get('partial_formula'):
            size_dict = dict()
            size_dict['stoichiometry'] = dict()
            size_dict['stoichiometry']['$size'] = len(elements)
            query_dict['$and'].append(size_dict)
        return query_dict
    
    def query_composition(self, custom_elem=None):
        """ Query DB for all structures containing 
        all the elements taken as input.
        """
        if custom_elem == None:
            elements = self.args.get('composition')
        else:
            elements = custom_elem
        # if there's only one string, try split it by caps
        numeracy = False
        if len(elements) == 1:
            elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements[0]) if elem]
            if elements[0].isdigit():
                numeracy = True
        try:
            if numeracy == False:
                for elem in elements:
                    if bool(re.search(r'\d', elem)):
                        raise RuntimeError('Composition string must be a list of elements or a single number.')
            elif numeracy:
                if self.args.get('partial_formula'):
                    raise RuntimeError('Number of elements not compatible with partial formula.')
        except Exception as oops:
            print(oops)
            return EmptyCursor()
        query_dict = dict()
        query_dict['$and'] = []
        if not numeracy:
            for ind, elem in enumerate(elements):
                types_dict = dict()
                types_dict['atom_types'] = dict()
                types_dict['atom_types']['$in'] = [elem]
                query_dict['$and'].append(types_dict)
        if not self.args.get('partial_formula'):
            size_dict = dict()
            size_dict['stoichiometry'] = dict()
            if numeracy:
                num = int(elements[0])
            else:
                num = len(elements)
            size_dict['stoichiometry']['$size'] = num
            query_dict['$and'].append(size_dict)
        return query_dict
 
    def query_tags(self):
        """ Find all structures matching given tags. """
        query_dict = dict()
        query_dict['$and'] = []
        for tag in self.args.get('tags'):
            temp_dict = dict()
            temp_dict['tags'] = dict()
            temp_dict['tags']['$in'] = [tag]
            query_dict['$and'].append(temp_dict)
        return query_dict

    def query_user(self):
        """ Find all structures matching given tags. """
        query_dict = dict()
        query_dict['$and'] = []
        query_dict['$and'].append(dict())
        query_dict['$and'][-1]['$or'] = []
        for user in self.args.get('user'):
            query_dict['$and'][-1]['$or'].append([dict('user', user)])
        return query_dict

    def query_pressure(self):
        """ Query pressure, either by an exact match on external_pressure
        or an approximate match on the pressure on the cell, with tolerance
        of either 0.05 GPa or 10%.
        """
        input_pressure = self.args.get('pressure')
        print(input_pressure, 'GPa')
        approx_pressure = [0.9*input_pressure-0.05, 1.1*input_pressure+0.05]
        query_dict = dict()
        query_dict['$and'] = []
        temp_dict = dict()
        temp_dict['external_pressure'] = dict()
        temp_dict['external_pressure']['$in'] = [[input_pressure]]
        query_dict['$and'].append(temp_dict)
        temp_dict = dict()
        temp_dict['pressure'] = dict()
        temp_dict['pressure']['$lt'] = approx_pressure[1]
        query_dict['$and'].append(temp_dict)
        temp_dict['pressure'] = dict()
        temp_dict['pressure']['$gt'] = approx_pressure[0]
        query_dict['$and'].append(temp_dict)
        return query_dict
 
    def query_calc(self, cursor):
        """ Find all structures with matching
        accuracy to specified structure.
        """
        doc = cursor[0]
        self.gs_enthalpy = doc['enthalpy_per_atom']
        query_dict = []
        # query_dict['$and'] = []
        temp_dict = dict()
        temp_dict['xc_functional'] = doc['xc_functional']
        query_dict.append(temp_dict)
        temp_dict = dict()
        temp_dict['cut_off_energy'] = doc['cut_off_energy']
        # temp_dict['cut_off_energy'] = dict()
        # temp_dict['cut_off_energy']['$gqe'] = doc['cut_off_energy']
        # query_dict[-1]['cut_off_energy'] = temp_dict
        query_dict.append(temp_dict)
        temp_dict = dict()
        temp_dict['species_pot'] = doc['species_pot']
        query_dict.append(temp_dict)
        # for species in doc['species_pot']:
            # temp_dict = dict()
            # temp_dict['species_pot'] = dict()
            # temp_dict['$in'] = doc['species_pot'][species]
        temp_dict = dict()
        temp_dict['kpoints_mp_spacing'] = doc['kpoints_mp_spacing']
        query_dict.append(temp_dict)
        # query_dict[-1]['kpoints_mp_spacing'] = temp_dict
        return query_dict

    def print_report(self):
        """ Print spatula report on current database. """
        try:
            report = self.report.find_one()
            print('Database last modified on', report['last_modified'], 'with spatula', report['version'], 
                  'changeset (' + report['git_hash'] + ').')
        except:
            print('Failed to print database report: spatula is probably running!')
    
    def dbstats(self):
        """ Print some useful stats about the database. """ 
        db_stats_dict = self.db.command('collstats', self.repo.name)
        print('Database collection', self.db.name + '.' + self.repo.name, 'contains', db_stats_dict['count'],
              'structures at', "{:.1f}".format(db_stats_dict['avgObjSize']/1024), 'kB each, totalling', 
              "{:.1f}".format(db_stats_dict['storageSize']/(1024**2)), 'MB when padding is included.')
        cursor = self.repo.find()
        comp_list = dict()
        for doc in cursor:
            temp = ''
            for ind, elem in enumerate(doc['stoichiometry']):
                temp += str(elem[0])
                if ind != len(doc['stoichiometry'])-1:
                    temp += '+'
            if temp not in comp_list:
                comp_list[temp] = 0
            comp_list[temp] += 1
        keys = list(comp_list.keys())
        vals = list(comp_list.values())
        comp_list = zip(keys, vals)
        comp_list.sort(key = lambda t:t[1], reverse=True)
        small_list = []
        small_count = 0
        first_ind = 1000
        cutoff = 200
        for ind, comp in enumerate(comp_list):
            if comp[1] < cutoff:
                if ind < first_ind:
                    first_ind = ind
                small_list.append(comp[0])
                small_count += comp[1]
        comp_list = comp_list[:first_ind] 
        comp_list.append(['others < ' + str(cutoff), small_count])
        comp_list.sort(key = lambda t:t[1], reverse=True)
        from ascii_graph import Pyasciigraph
        from ascii_graph.colors import Gre, Blu, Yel, Red
        from ascii_graph.colordata import hcolor
        graph = Pyasciigraph(line_length=80, multivalue=False)
        db_stats_dict['count']
        thresholds = {int(db_stats_dict['count']/40): Gre, int(db_stats_dict['count']/10): Blu, int(db_stats_dict['count']/4): Red,}
        data = hcolor(comp_list, thresholds)
        for line in graph.graph(label=None, data=data):
           print(line) 
        print('where others', end=': ')
        for small in small_list:
            print(small, end=', ')
        print('\n') 

    def temp_collection(self, cursor):
        """ Create temporary collection
        for successive filtering. 
        """
        # check temp doesn't already exist; drop if it does
        try:
            self.client.crystals.temp.drop()
        except:
            pass
        self.temp = self.client.crystals.temp
        if cursor.count() != 0:
            self.temp.insert(cursor)
        else:
            self.temp.drop()
            exit('No structures found.')
        return self.temp

class EmptyCursor:
    """ Empty cursor class for failures. """
    def count(self):
        return 0 
    def clone(self):
        return EmptyCursor()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Query MongoDB structure database.',
            epilog='Written by Matthew Evans (2016). Based on the cryan concept by Chris Pickard.')
    group = parser.add_argument_group()
    group.add_argument('-f', '--formula', nargs='+', type=str,
        help='choose a stoichiometry, e.g. Ge 1 Te 1 Si 3, or GeTeSi3')
    group.add_argument('-c', '--composition', nargs='+', type=str,
        help=('find all structures containing the given elements, e.g. GeTeSi, or find' +
        'the number of structures with n elements, e.g. 1, 2, 3'))
    group.add_argument('-i', '--id', type=str, nargs='+',
            help='specify a particular structure by its text_id')
    parser.add_argument('-s', '--summary', action='store_true',
            help='show only the ground state for each formula (i.e. phase+stoichiometry)')
    parser.add_argument('-t', '--top', type=int,
            help='number of structures to show (DEFAULT: 10)')
    parser.add_argument('-d', '--details', action='store_true',
            help='show as much detail about calculation as possible')
    parser.add_argument('-p', '--pressure', type=float,
            help='specify an isotropic external pressure to search for, e.g. 10 (GPa)')
    parser.add_argument('--source', action='store_true',
            help='print filenames from which structures were wrangled')
    parser.add_argument('-ac', '--calc-match', action='store_true',
            help='display calculations of the same accuracy as specified id')
    parser.add_argument('-pf', '--partial-formula', action='store_true',
            help=('stoichiometry/composition queries will include other unspecified species,' +
            'e.g. search for Ge2Te3 will also return Ge2Te3Si or Ge4Te6Fe2, and search for Li will' +
            ' include any structure containing Li, not just pure Li.'))
    parser.add_argument('--dbstats', action='store_true',
            help=('print some stats about the database that is being queried'))
    parser.add_argument('--user', nargs=1, type=str,
            help=('search for a CRSid (or other user tag) to find a particular user\'s structures'))
    parser.add_argument('--tags', nargs='+', type=str,
            help=('search for up to 3 manual tags at once'))
    parser.add_argument('--hull', action='store_true',
            help=('create a convex hull for 2 elements (to be extended to 3, 4 soon)'))
    parser.add_argument('--dis', action='store_true',
            help='smear hull with local stoichiometry')
    parser.add_argument('--scratch', action='store_true',
            help=('query local scratch database'))
    parser.add_argument('--cell', action='store_true',
            help='export query to .cell files in folder name from query string')
    parser.add_argument('--res', action='store_true',
            help='export query to .res files in folder name from query string')
    parser.add_argument('--write_pressure', nargs='+', type=str,
            help='pressure to add to new cell file, either one float' +  
            'for isotropic or 6 floats for anisotropic.')
    args = parser.parse_args()
    if args.calc_match and args.id == None:
        exit('--calc-match requires -i or --id')
    if args.hull == True and args.composition == None:
        exit('--hull requires --composition')
    if args.dis == True and args.hull != True:
        exit('--dis requires --hull')
    if args.write_pressure and args.cell != True:
        exit('--write_pressure requires cell')
    query = DBQuery(stoichiometry=args.formula, composition=args.composition,
                    summary=args.summary, id=args.id, top=args.top, details=args.details,
                    pressure=args.pressure, source=args.source, calc_match=args.calc_match,
                    partial_formula=args.partial_formula, dbstats=args.dbstats, scratch=args.scratch, 
                    tags=args.tags, user=args.user, hull=args.hull, dis=args.dis, res=args.res,
                    cell=args.cell, write_pressure=args.write_pressure, main=True, sysargs=argv[1:])
    # generate hull outside query object
    if args.hull:
        from hull import FryanConvexHull
        FryanConvexHull(query)
