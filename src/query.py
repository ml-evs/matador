#!/usr/bin/python
# coding: utf-8
""" This file implements all queries to the database,
including parsing user inputs, displaying results
and calling other functionality. """
from __future__ import print_function
# import related crysdb functionality
from export import query2files
# import external libraries
import pymongo as pm
import numpy as np
from bson.son import SON
# import standard library
import re
from os import uname


class DBQuery:
    """ Class that implements queries to MongoDB
    structure database.
    """
    def __init__(self, client=False, collections=False, *args, **kwargs):
        """ Parse arguments from matador or API call
        before calling query.
        """
        # read args
        try:
            self.args = args[0]
            self.client = client
            self.db = client.crystals
            self.collections = collections
        except:
            pass
        # if empty collections, assume called from API and read kwargs, 
        # also need to connect to db
        if not collections or not client:
            local = uname()[1]
            if local == 'cluster2':
                remote = 'node1'
            else:
                remote = None
            self.client = pm.MongoClient(remote)
            self.db = self.client.crystals
            self.args = dict()
            self.collections = dict()
            if kwargs['db'] is not None:
                for database in kwargs['db']:
                    if database == 'all':
                        self.collections['ajm'] = self.db['repo']
                        self.collections['oqmd'] = self.db['oqmd']
                    elif database == 'ajm':
                        database = 'repo'
                        self.collections['ajm'] = self.db['repo']
                    else:
                        self.collections[database] = self.db[database]
            else:
                self.collections['ajm'] = self.db['repo']
            self.args['tags'] = kwargs['tags']
        if self.args.get('summary'):
            self.top = -1
        else:
            self.top = self.args.get('top') if self.args.get('top') is not None else 10
        # define some periodic table macros
        self.periodic_table = dict()
        self.periodic_table['I'] = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
        self.periodic_table['II'] = ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']
        self.periodic_table['III'] = ['B', 'Al', 'Ga', 'In', 'Tl']
        self.periodic_table['IV'] = ['C', 'Si', 'Ge', 'Sn', 'Pb']
        self.periodic_table['V'] = ['N', 'P', 'As', 'Sb', 'Bi']
        self.periodic_table['VI'] = ['O', 'S', 'Se', 'Te', 'Po']
        self.periodic_table['VII'] = ['F', 'Cl', 'Br', 'I', 'At']
        self.periodic_table['Tran'] = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                                       'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                                       'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
        self.periodic_table['Lan'] = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                                      'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
        self.periodic_table['Act'] = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
                                      'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

        self.perform_query()

    def perform_query(self):
        """ Set up query dict and perform query depending on 
        command-line / API arguments. 
        """
        self.cursor = EmptyCursor()
        # initalize query_dict to '$and' all queries
        self.query_dict = dict()
        self.query_dict['$and'] = []
        empty_query = True
        # benchmark enthalpy to display (set by calc_match)
        self.gs_enthalpy = 0.0
        if self.args.get('id') is not None:
            self.cursor = []
            for collection in self.collections:
                temp_cursor = self.collections[collection].find({'text_id': self.args.get('id')})
                for doc in temp_cursor:
                    self.cursor.append(doc)
            if len(self.cursor) < 1:
                exit('Could not find a match, try widening your search.')
            elif len(self.cursor) == 1:
                if self.args.get('cell') or self.args.get('res'):
                    query2files(self.cursor, self.args)
                self.display_results(self.cursor)
            if self.args.get('calc_match'):
                # save special copy of calc_dict for hulls
                self.calc_dict = dict()
                self.calc_dict['$and'] = []
                # to avoid deep recursion, and since this is always called first
                # don't append, just set
                self.query_dict['$and'] = self.query_calc(self.cursor[0])
                self.calc_dict['$and'] = list(self.query_dict['$and'])
                empty_query = False
        # create alias for formula for backwards-compatibility
        self.args['stoichiometry'] = self.args.get('formula')
        if self.args.get('stoichiometry') is not None:
            self.query_dict['$and'].append(self.query_stoichiometry())
            empty_query = False
        if self.args.get('composition') is not None:
            self.query_dict['$and'].append(self.query_composition())
            empty_query = False
        if self.args.get('num_species') is not None:
            self.query_dict['$and'].append(self.query_num_species())
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
            if self.args.get('id') is None:
                for collection in self.collections:
                    self.repo = self.collections[collection]
                    if self.args.get('debug'):
                        print('Empty query, showing all...')
                    self.cursor = self.repo.find().sort('enthalpy_per_atom', pm.ASCENDING)
                    if self.top == -1:
                        self.top = self.cursor.count()
                    self.display_results(self.cursor[:self.top])
        # if no special query has been made already, begin executing the query
        if not empty_query:
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
                        query2files(self.cursor, self.args)
                # if called as script, always print results
                if self.args.get('id') is None:
                    print(cursor_count, 'results found for query in', collection+'.')
                if self.args.get('subcmd') != 'hull':
                    if cursor_count >= 1:
                        if self.top == -1:
                            self.top = cursor_count
                        if cursor_count > self.top:
                            self.display_results(self.cursor.clone()[:self.top])
                        else:
                            self.display_results(self.cursor.clone())
            # building hull from just comp, find best structure to calc_match
            if self.args.get('subcmd') == 'hull' or self.args.get('subcmd') == 'voltage' or self.args.get('hull_cutoff') is not None:
                if 'oqmd' in self.collections:
                    exit('Use --include_oqmd instead of --db, exiting...')
                if len(self.collections.keys()) == 1:
                    self.repo = self.collections[self.collections.keys()[0]]
                else:
                    exit('Hulls and voltage curves require just one source or --include_oqmd, exiting...')
                print('Creating hull from AJM db structures.')
                self.args['summary'] = True
                self.args['top'] = -1
                if self.args.get('biggest'):
                    print('\nFinding biggest calculation set for hull...\n')
                else:
                    print('\nFinding the best calculation set for hull...')
                test_cursors = []
                test_cursor_count = []
                test_query_dict = []
                calc_dicts = []
                cutoff = []
                sample = 5
                rand_sample = 5 if self.args.get('biggest') else 2
                i = 0
                count = self.cursor.count()
                if count <= 0:
                    exit('No structures found for hull.')
                while i < sample+rand_sample:
                    # start with sample/2 lowest enthalpy structures
                    if i < int(sample):
                        ind = i
                    # then do some random samples
                    else:
                        ind = np.random.randint(rand_sample, count-1)
                    id_cursor = self.repo.find({'text_id': self.cursor[ind]['text_id']})
                    self.query_dict = dict()
                    self.query_dict['$and'] = self.query_calc(id_cursor[0])
                    cutoff.append(id_cursor[0]['cut_off_energy'])
                    calc_dicts.append(dict())
                    calc_dicts[-1]['$and'] = list(self.query_dict['$and'])
                    self.query_dict['$and'].append(self.query_composition())
                    test_query_dict.append(self.query_dict)
                    test_cursors.append(
                        self.repo.find(SON(test_query_dict[-1])).sort('enthalpy_per_atom',
                                                                      pm.ASCENDING))
                    test_cursor_count.append(test_cursors[-1].count())
                    print("{:^24}".format(self.cursor[ind]['text_id'][0] + ' ' +
                                          self.cursor[ind]['text_id'][1]) +
                          ': matched ' + str(test_cursor_count[-1]), 'structures.')
                    try:
                        print(12*' ' + '└╌╌╌╌╌╌╌╌╌╌╌╌╌ ', self.cursor[ind]['xc_functional'] + ',',
                              self.cursor[ind]['cut_off_energy'], 'eV,', self.cursor[ind]['kpoints_mp_spacing'], '1/A')
                    except:
                        pass
                    if self.args.get('biggest'):
                        if test_cursor_count[-1] > 2*int(count/3):
                            print('Matched at least 2/3 of total number, composing hull...')
                            break
                    i += 1
                if self.args.get('biggest'):
                    choice = np.argmax(np.asarray(test_cursor_count))
                else:
                    # by default, find highest cutoff hull as first proxy for quality
                    choice = np.argmax(np.asarray(cutoff))
                self.cursor = test_cursors[choice]
                self.calc_dict = calc_dicts[choice]
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

    def display_results(self, cursor):
        """ Print query results in a cryan-like fashion. """
        details = self.args.get('details')
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
                struct_string[-1] += "{:^12.3f}".format(doc['cell_volume']/doc['num_fu'])
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
            try:
                struct_string[-1] += "{:^8}".format(doc['num_fu'])
            except:
                struct_string[-1] += "{:^8}".format('xxx')
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
                    detail_string[-1] += 'xc-functional unknown!'
                if 'cut_off_energy' in doc:
                    detail_string[-1] += ', ' + "{:4.2f}".format(doc['cut_off_energy']) + ' eV'
                else:
                    detail_string[-1] += 'cutoff unknown'
                if 'external_pressure' in doc:
                    detail_string[-1] += (', ' +
                                          "{:4.2f}".format(doc['external_pressure'][0][0]) +
                                          ' GPa')
                if 'kpoints_mp_spacing' in doc:
                    detail_string[-1] += ', ~' + str(doc['kpoints_mp_spacing']) + ' 1/A'
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
                    source_string[-1] += ' ' + file.split('structure_repository')[-1]
                    if num != len(doc['source'])-1:
                        source_string[-1] += '\n'
        print(len(header_string)*'─')
        print(header_string)
        print(len(header_string)*'─')
        if self.args.get('summary'):
            current_formula = ''
            formula_list = []
            count = 0
            for ind, substring in enumerate(formula_string):
                if substring != current_formula and substring not in formula_list:
                    count += 1
                    print(struct_string[ind])
                    if details:
                        print(detail_string[ind])
                        print(detail_substring[ind])
                    if self.args.get('source'):
                        print(source_string[ind])
                    current_formula = substring
                    formula_list.append(substring)
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
        all the elements taken as input. Passing this
        function a number is a deprecated feature, replaced
        by query_num_species.
        """
        if custom_elem is None:
            elements = self.args.get('composition')
        else:
            elements = custom_elem
        # if there's only one string, try split it by caps
        try:
            if len(elements) == 1:
                elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements[0]) if elem]
                if elements[0].isdigit():
                    raise RuntimeError('Composition string must be a ' +
                                       'list of elements, use --num_species.')
                while '[' in elements:
                    tmp_elements = list(elements)
                    for ind, tmp in enumerate(tmp_elements):
                        if tmp == '[':
                            while tmp_elements[ind+1] != ']':
                                tmp_elements[ind] += tmp_elements[ind+1]
                                del tmp_elements[ind+1]
                            tmp_elements[ind] += ']'
                    try:
                        tmp_elements.remove(']')
                    except:
                        pass
                    try:
                        tmp_elements.remove('')
                    except:
                        pass
                    elements = tmp_elements
                    for i in range(len(elements)):
                        elements = elements[i].split('][')
                for elem in elements:
                    if bool(re.search(r'\d', elem)):
                        raise RuntimeError('Composition string must be a ' +
                                           'list of elements or a single number.')
        except Exception as oops:
            print(type(oops), oops)
            return EmptyCursor()
        query_dict = dict()
        query_dict['$and'] = []
        for ind, elem in enumerate(elements):
            # prototype for chemically motivated searches, e.g. transition metals
            if '[' in elem or ']' in elem:
                types_dict = dict()
                types_dict['$or'] = list()
                elem = elem.strip('[').strip(']')
                for group_elem in self.periodic_table[elem]:
                    types_dict['$or'].append(dict())
                    types_dict['$or'][-1]['atom_types'] = group_elem
            else:
                types_dict = dict()
                types_dict['atom_types'] = dict()
                types_dict['atom_types']['$in'] = [elem]
            query_dict['$and'].append(types_dict)
        if not self.args.get('partial_formula'):
            size_dict = dict()
            size_dict['stoichiometry'] = dict()
            num = len(elements)
            size_dict['stoichiometry']['$size'] = num
            query_dict['$and'].append(size_dict)
        return query_dict

    def query_num_species(self):
        """ Query database for all structures with a
        given number of elements, e.g. binaries, ternaries etc.
        """
        num = self.args.get('num_species')
        if len(num) != 1:
            exit('--num_species takes a single integer')
        else:
            num = num[0]
        query_dict = dict()
        query_dict['stoichiometry'] = dict()
        query_dict['stoichiometry']['$size'] = num
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
        # if missing xc, return dict that will have no matches
        if 'xc_functional' not in doc:
            query_dict.append(dict())
            query_dict[-1]['xc_functional'] = 'missing'
        temp_dict = dict()
        temp_dict['xc_functional'] = doc['xc_functional']
        query_dict.append(temp_dict)
        temp_dict = dict()
        if 'spin_polarized' in doc and doc['spin_polarized'] != []:
            temp_dict['spin_polarized'] = doc['spin_polarized']
            query_dict.append(temp_dict)
        else:
            temp_dict['spin_polarized'] = dict()
            temp_dict['spin_polarized']['$ne'] = True
            query_dict.append(temp_dict)
        if self.args.get('loose'):
            temp_dict = dict()
            query_dict.append(dict())
            query_dict[-1]['cut_off_energy'] = doc['cut_off_energy']
            # temp_dict = dict()
            # query_dict.append(dict())
            # temp_dict['$lte'] = float(doc['kpoints_mp_spacing'])
            # query_dict[-1]['kpoints_mp_spacing'] = temp_dict
        else:
            temp_dict = dict()
            temp_dict['kpoints_mp_spacing'] = dict()
            temp_dict['kpoints_mp_spacing']['$gt'] = doc['kpoints_mp_spacing'] - 0.02
            temp_dict['kpoints_mp_spacing']['$lt'] = doc['kpoints_mp_spacing'] + 0.02
            query_dict.append(temp_dict)
            query_dict.append(dict())
            query_dict[-1]['cut_off_energy'] = doc['cut_off_energy']
        if 'species_pot' in doc:
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
