#!/usr/bin/python
# coding: utf-8
""" This file implements all queries to the database,
including parsing user inputs, displaying results
and calling other functionality. """
from __future__ import print_function
# import related crysdb functionality
from spatula import param2dict
# import external libraries
import pymongo as pm
import numpy as np
import string
from os import makedirs, system
from os.path import exists, isfile, expanduser
from copy import deepcopy
from bson.son import SON
import re


class DBQuery:
    """ Class that implements queries to MongoDB
    structure database.
    """
    def __init__(self, client, collections, *args):
        """ Initialise the query with command line
        arguments and return results.
        """
        # read args
        self.args = args[0]
        self.client = client
        self.db = client.crystals
        self.collections = collections
        # print(self.args)
        # for arg in self.args:
            # if type(self.args[arg]) == str:
                # self.args[arg] = self.args[arg].split()
        # if self.args.get('scratch'):
            # self.repo = self.client.crystals.scratch
        # elif self.args.get('oqmd'):
            # self.repo = self.client.crystals.oqmd
        # else:
            # self.repo = self.client.crystals.repo
        # print last spatula report
        if self.args.get('summary'):
            self.top = -1
        else:
            self.top = self.args.get('top') if self.args.get('top') is not None else 10
        # grab all args as string for file dumps; TO-DO: make this better
        """ PERFORM QUERY """
        self.cursor = EmptyCursor()
        # initalize query_dict to '$and' all queries
        self.query_dict = dict()
        self.query_dict['$and'] = []
        empty_query = True
        # benchmark enthalpy to display (set by calc_match)
        self.gs_enthalpy = 0.0
        if self.args.get('id') is not None:
            empty_query = False
            self.cursor = []
            for collection in self.collections:
                temp_cursor = self.collections[collection].find({'text_id': self.args.get('id')})
                for doc in temp_cursor:
                    self.cursor.append(doc)
            if len(self.cursor) < 1:
                exit('Could not find a match, try widening your search.')
            elif self.args.get('calc_match'):
                # save special copy of calc_dict for hulls
                self.calc_dict = dict()
                self.calc_dict['$and'] = []
                # to avoid deep recursion, and since this is always called first
                # don't append, just set
                self.query_dict['$and'] = self.query_calc(self.cursor[0])
                self.calc_dict['$and'] = list(self.query_dict['$and'])
                empty_query = False
            else:
                self.display_results(self.cursor)
        if self.args.get('stoichiometry') is not None:
            self.query_dict['$and'].append(self.query_stoichiometry())
            empty_query = False
        if self.args.get('composition') is not None:
            self.query_dict['$and'].append(self.query_composition())
            empty_query = False
        if self.args.get('pressure') is not None:
            self.query_dict['$and'].append(self.query_pressure())
            empty_query = False
        if self.args.get('space_group') is not None:
            self.query_dict['$and'].append(self.query_space_group())
            empty_query = False
        if self.args.get('num_fu') is not None:
            self.query_dict['$and'].append(self.query_num_fu())
            empty_query = False
        if self.args.get('encapsulated') is True:
            self.query_dict['$and'].append(self.query_encap())
            empty_query = False
        if self.args.get('tags') is not None:
            self.query_dict['$and'].append(self.query_tags())
            empty_query = False
        # only query quality when making a hull
        if not self.args.get('ignore_warnings'):
            self.query_dict['$and'].append(self.query_quality())
        # if no query submitted, find all
        if empty_query:
            for collection in self.collections:
                self.repo = self.collections[collection]
                if self.args.get('debug'):
                    print('Empty query, showing all...')
                self.cursor = self.repo.find().sort('enthalpy_per_atom', pm.ASCENDING)
                if self.top == -1:
                    self.top = self.cursor.count()
                self.display_results(self.cursor[:self.top], details=self.args.get('details'))
        # if no special query has been made already, begin executing the query
        if not empty_query or self.args.get('calc_match'):
            self.cursors = []
            for collection in self.collections:
                self.repo = self.collections[collection]
                if self.args.get('details'):
                    print(self.query_dict)
                # execute query
                self.cursor = self.repo.find(SON(self.query_dict)).sort('enthalpy_per_atom',
                                                                        pm.ASCENDING)
                self.cursors.append(self.cursor.clone())
                cursor_count = self.cursor.count()
                """ QUERY POST-PROCESSING """
                # write query to res or cell with param files
                if self.args.get('cell') or self.args.get('res'):
                    if cursor_count >= 1:
                        if self.args.get('top') is not None:
                            if self.top == -1:
                                self.top = cursor_count
                            self.query2files(self.cursor[:self.top],
                                             self.args.get('res'),
                                             self.args.get('cell'),
                                             top=True,
                                             pressure=self.args.get('write_pressure'))
                        else:
                            self.query2files(self.cursor,
                                             self.args.get('res'),
                                             self.args.get('cell'),
                                             pressure=self.args.get('write_pressure'))
                # if called as script, always print results
                if self.args.get('id') is None:
                    print(cursor_count, 'results found for query in', collection+'.')
                if self.args.get('subcmd') != 'hull':
                    if cursor_count >= 1:
                        if self.top == -1:
                            self.top = cursor_count
                        if cursor_count > self.top:
                            self.display_results(self.cursor.clone()[:self.top],
                                                 details=self.args.get('details'))
                        else:
                            self.display_results(self.cursor.clone(),
                                                 details=self.args.get('details'))
            # building hull from just comp, find best structure to calc_match
            if self.args.get('subcmd') == 'hull' or self.args.get('subcmd') == 'voltage':
                if 'repo' in self.collections:
                    self.repo = self.collections['repo']
                else:
                    exit('Hulls and voltage curves require AJM repo structures, exiting...')
                print('Creating hull from AJM db structures.')
                self.args['summary'] = True
                print('\nFinding biggest calculation set for hull...\n')
                test_cursor = []
                test_cursor_count = []
                test_query_dict = []
                sample = 10
                rerun = False
                i = 0
                count = self.cursor.count()
                while i < sample:
                    # start with sample/2 lowest enthalpy structures
                    if i < int(sample/2):
                        ind = i
                    # then do some random samples
                    else:
                        ind = np.random.randint(5, count-1)
                    id_cursor = self.repo.find({'text_id': self.cursor[ind]['text_id']})
                    self.query_dict['$and'] = self.query_calc(id_cursor[0])
                    self.calc_dict = dict()
                    self.calc_dict['$and'] = list(self.query_dict['$and'])
                    self.query_dict['$and'].append(self.query_composition())
                    test_query_dict.append(self.query_dict)
                    test_cursor.append(
                        self.repo.find(SON(test_query_dict[-1])).sort('enthalpy_per_atom',
                                                                      pm.ASCENDING))
                    test_cursor_count.append(test_cursor[-1].count())
                    print("{:^24}".format(self.cursor[ind]['text_id'][0] + ' ' +
                                          self.cursor[ind]['text_id'][1]) +
                          ': matched' + str(test_cursor_count[-1]), 'structures.')
                    # if we have at least 2/3 of the structures, just plot
                    if test_cursor_count[-1] > 2*int(count/3):
                        print('Matched at least 2/3 of total number, composing hull...')
                        break
                    if i == (sample-1) and not rerun:
                        # if less than half the total structures are to be plotted, rand 5 more
                        if np.max(np.asarray(test_cursor_count)) < int(count/2):
                            i -= 5
                            rerun = True
                    i += 1
                self.cursor = test_cursor[np.argmax(np.asarray(test_cursor_count))]
                # if including oqmd, connect to oqmd collection and generate new query
                if self.args.get('include_oqmd'):
                    self.oqmd_repo = self.client.crystals.oqmd
                    self.oqmd_query = dict()
                    self.oqmd_query['$and'] = []
                    # query only oqmd, assume all calculations are over-converged:
                    # this is bad!
                    self.oqmd_query['$and'].append(self.query_composition())
                    self.oqmd_cursor = self.oqmd_repo.find(SON(self.oqmd_query))
                    self.oqmd_cursor.sort('enthalpy_per_atom', pm.ASCENDING)

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
        name = name + '-' + str(pairs) + '-pair-swaps/' + name
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
        """ Quick hash generator, based on implementation in PyAIRSS by J. Wynn. """
        hashChars = hashChars = [str(x) for x in range(0, 10)]+[x for x in string.ascii_lowercase]
        hash = ''
        for i in range(hashLen):
            hash += np.random.choice(hashChars)
        return hash

    def query2files(self, cursor, res=False, cell=False, top=False, pressure=None):
        """ Write .res or .cell files for all docs in query,
        including a .param file for each. Eventually handle
        swaps from one element to another from CLI.
        """
        if cursor.count() > 1000 and top is False:
            write = raw_input('This operation will write ' + str(cursor.count()) + ' structures,' +
                              ' are you sure you want to do this? [y/n] ')
            if write == 'y' or write == 'Y':
                print('Writing them all.')
                write = True
            else:
                write = False
                return
        else:
            write = True
        name = 'query-'
        if self.args['composition'] is not None:
            for comp in self.args['composition']:
                name += comp
        elif self.args['formula'] is not None:
            name += self.args['formula']
        name += '-' + self.args['db'][0]
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
        if template is not None:
            try:
                paramDict, success = param2dict(template)
                if not success:
                    raise RuntimeError('Failed to open template.')
            except Exception as oops:
                print(type(oops), oops)
        try:
            if isfile(path+'.param'):
                print('File name already exists, generating hash...')
                path += '-' + self.generate_hash()
        except Exception as oops:
            print('Writing param file failed for ', doc['text_id'])
            print(type(oops), oops)
        try:
            with open(path+'.param', 'w') as f:
                f.write('# Param file generated by matador (Matthew Evans 2016)\n')
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
                f.write('# Cell file generated by matador (Matthew Evans 2016)\n\n')
                f.write('# enthalpy_per_atom = ' +
                        '{: 10f} eV\n\n'.format(doc['enthalpy_per_atom']))
                f.write('%BLOCK LATTICE_CART\n')
                for vec in doc['lattice_cart']:
                    for coeff in vec:
                        f.write(str(coeff) + ' ')
                    f.write('\n')
                f.write('%ENDBLOCK LATTICE_CART\n\n')
                f.write('%BLOCK POSITIONS_FRAC\n')
                for ind, atom in enumerate(zip(doc['atom_types'], doc['positions_frac'])):
                    f.write("{0:8s} {1[0]: 15f} {1[1]: 15f} {1[2]: 15f}   1.0\n".format(atom[0],
                            atom[1]))
                f.write('%ENDBLOCK POSITIONS_FRAC\n\n')
                if pressure is not None:
                    f.write('%block external_pressure\n'.upper())
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
                    f.write('kpoints_mp_grid : ' + str(doc['kpoints_mp_grid'][0]) + ' ' +
                            str(doc['kpoints_mp_grid'][1]) + ' ' +
                            str(doc['kpoints_mp_grid'][2]) + '\n')
                if 'species_pot' in doc:
                    f.write('\n%BLOCK SPECIES_POT\n')
                    for elem in doc['species_pot']:
                        if not isfile(''.join(path.split('/')[:-1])+'/'+doc['species_pot'][elem]):
                            if isfile(expanduser('~/pspot/' + doc['species_pot'][elem])):
                                system('cp ' + expanduser('~/pspot/') + doc['species_pot'][elem] +
                                       ' ' + ''.join(path.split('/')[:-1]))
                        f.write(elem + '\t' + doc['species_pot'][elem] + '\n')
                    f.write('%ENDBLOCK SPECIES_POT')
        except Exception as oops:
            print('Writing cell file failed for ', doc['text_id'])
            print(type(oops), oops)

    def doc2res(self, doc, path):
        """ Write .res file for single doc. """
        try:
            if isfile(path+'.res'):
                print('File name already exists, generating hash...')
                path += '-' + self.generate_hash()
            with open(path+'.res', 'w') as f:
                # f.write('# Res file generated by matador (Matthew Evans 2016)\n\n')
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
            print('Writing res file failed for ', doc['text_id'])
            print(type(oops), oops)

    def display_results(self, cursor, details=False):
        """ Print query results in a cryan-like fashion. """
        struct_string = []
        detail_string = []
        detail_substring = []
        source_string = []
        formula_string = []
        last_formula = ''
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
            if 'encapsulated' in doc:
                formula_substring += '+CNT'
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
                struct_string[-1] += "{:^12.3f}".format(atom_per_fu *
                                                        doc['cell_volume'] / doc['num_atoms'])
            except:
                struct_string[-1] += "{:^12}".format('xxx')
            try:
                struct_string[-1] += "{:^18.5f}".format(doc['enthalpy_per_atom'] -
                                                        self.gs_enthalpy)
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
                    detail_string[-1] += (', ' +
                                          "{:4.2f}".format(doc['external_pressure'][0][0]) +
                                          ' GPa')
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
                if 'encapsulated' in doc:
                    try:
                        detail_string[-1] += (', (n,m)=(' + str(doc['cnt_chiral'][0]) +
                                              ',' + str(doc['cnt_chiral'][1]) + ')')
                        detail_string[-1] += ', r=' + "{:4.2f}".format(doc['cnt_radius']) + ' A'
                        detail_string[-1] += ', z=' + "{:4.2f}".format(doc['cnt_length']) + ' A'
                    except:
                        pass
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
            for ind, substring in enumerate(formula_string):
                if count > self.top:
                    break
                if substring != current_formula and substring not in formula_string[:ind]:
                    count += 1
                    print(struct_string[ind])
                    if details:
                        print(detail_string[ind])
                        print(detail_substring[ind])
                    if self.args.get('source'):
                        print(source_string[ind])
                    current_formula = substring
        else:
            for ind, substring in enumerate(struct_string):
                print(substring)
                if details:
                    print(detail_string[ind])
                    print(detail_substring[ind])
                if self.args.get('source'):
                    print(source_string[ind])
                if details or self.args.get('source'):
                    print(len(header_string) * '─')

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
        if custom_elem is None:
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
            if numeracy is False:
                for elem in elements:
                    if bool(re.search(r'\d', elem)):
                        raise RuntimeError('Composition string must be a ' +
                                           'list of elements or a single number.')
            elif numeracy:
                if self.args.get('partial_formula'):
                    raise RuntimeError('Number of elements not compatible with partial formula.')
        except Exception as oops:
            print(type(oops), oops)
            return EmptyCursor()
        query_dict = dict()
        query_dict['$and'] = []
        if not numeracy:
            for ind, elem in enumerate(elements):
                # prototype for chemically motivated searches, e.g. transition metals
                if elem == 'T':
                    types_dict = dict()
                    types_dict['$or'] = list()
                    transition_metals = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                                         'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                                         'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
                    for metal in transition_metals:
                        types_dict['$or'].append(dict())
                        types_dict['$or'][-1]['atom_types'] = metal
                else:
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

    def query_space_group(self):
        """ Query DB for all structures with given
        space group.
        """
        query_dict = dict()
        query_dict['space_group'] = str(self.args.get('space_group'))
        return query_dict

    def query_num_fu(self):
        """ Query DB for all structures with more than a
        given number of formula units in the simulation.
        """
        query_dict = dict()
        query_dict['num_fu'] = dict()
        query_dict['num_fu']['$gte'] = self.args.get('num_fu')
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

    def query_quality(self):
        """ Find all structures with non-zero or
        non-existent (e.g. OQMD) quality. """
        query_dict = dict()
        query_dict['$or'] = []
        query_dict['$or'].append(dict())
        query_dict['$or'][-1]['quality'] = dict()
        query_dict['$or'][-1]['quality']['$gt'] = 0
        query_dict['$or'].append(dict())
        query_dict['$or'][-1]['quality'] = dict()
        query_dict['$or'][-1]['quality']['$exists'] = False
        return query_dict

    def query_pressure(self):
        """ Query pressure, either by an exact match on external_pressure
        or an approximate match on the pressure on the cell, with tolerance
        of either 0.05 GPa or 10%.
        """
        input_pressure = self.args.get('pressure')
        print(input_pressure, 'GPa')
        if input_pressure < 0:
            approx_pressure = [1.1*input_pressure-0.05, 0.9*input_pressure+0.05]
        else:
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

    def query_encap(self):
        """ Query only CNT encapsulated structures. """
        query_dict = dict()
        query_dict['encapsulated'] = dict()
        query_dict['encapsulated']['$exists'] = True
        return query_dict

    def query_calc(self, doc):
        """ Find all structures with matching
        accuracy to specified structure.
        """
        self.gs_enthalpy = doc['enthalpy_per_atom']
        query_dict = []
        temp_dict = dict()
        temp_dict['xc_functional'] = doc['xc_functional']
        query_dict.append(temp_dict)
        temp_dict = dict()
        if 'spin_polarized' in doc:
            temp_dict['spin_polarized'] = doc['spin_polarized']
            query_dict.append(temp_dict)
        if self.args.get('loose'):
            temp_dict = dict()
            query_dict.append(dict())
            temp_dict['$gte'] = float(doc['cut_off_energy'])
            query_dict[-1]['cut_off_energy'] = temp_dict
            # temp_dict = dict()
            # query_dict.append(dict())
            # temp_dict['$lte'] = float(doc['kpoints_mp_spacing'])
            # query_dict[-1]['kpoints_mp_spacing'] = temp_dict
        else:
            # temp_dict = dict()
            # temp_dict['kpoints_mp_spacing'] = doc['kpoints_mp_spacing']
            # query_dict.append(temp_dict)
            query_dict.append(dict())
            query_dict[-1]['cut_off_energy'] = doc['cut_off_energy']
        for species in doc['species_pot']:
            temp_dict = dict()
            temp_dict['$or'] = []
            temp_dict['$or'].append(dict())
            temp_dict['$or'][-1]['species_pot.'+species] = dict()
            temp_dict['$or'][-1]['species_pot.'+species]['$exists'] = False
            temp_dict['$or'].append(dict())
            temp_dict['$or'][-1]['species_pot.'+species] = doc['species_pot'][species]
            query_dict.append(temp_dict)
        return query_dict

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
