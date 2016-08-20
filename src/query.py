# coding: utf-8
""" This file implements all queries to the database,
including parsing user inputs, displaying results
and calling other functionality. """
from __future__ import print_function
# import related matador functionality
from export import query2files
from print_utils import print_failure, print_warning
from chem_utils import get_periodic_table
# import external libraries
import pymongo as pm
import numpy as np
from bson.son import SON
from bson.json_util import dumps
# import standard library
import re
from os import uname
from itertools import combinations
from traceback import print_exc
from fractions import gcd


class DBQuery:
    """ Class that implements queries to MongoDB
    structure database.
    """
    def __init__(self, client=False, collections=False, **kwargs):
        """ Parse arguments from matador or API call
        before calling query.
        """
        # read args
        self.args = kwargs
        if client is not False:
            self.client = client
            self.db = client.crystals
        if collections is not False:
            self.collections = collections

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
            self.collections = dict()
            if self.args.get('db') is not None:
                if type(self.args['db']) is not list:
                    self.args['db'] = [self.args['db']]
                for database in self.args['db']:
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

        if self.args.get('summary'):
            self.top = -1
        else:
            self.top = self.args.get('top') if self.args.get('top') is not None else 10

        # define some periodic table macros
        self.periodic_table = get_periodic_table()

        # create the dictionary to pass to MongoDB
        self.construct_query()

        # execute the query
        self.perform_query()

    def construct_query(self):
        """ Set up query dict and perform query depending on
        command-line / API arguments.
        """
        self.cursor = EmptyCursor()

        # initalize query_dict to '$and' all queries
        self.query_dict = dict()
        self.query_dict['$and'] = []
        self.empty_query = True

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
            elif len(self.cursor) >= 1:
                if (self.args.get('cell') or self.args.get('res')) and \
                        not self.args.get('calc_match'):
                    query2files(self.cursor, self.args)
                self.display_results(self.cursor)
                if len(self.cursor) > 1:
                    print_warning('WARNING: matched multiple structures with same text_id. ' +
                                  'The first one will be used.')

            if self.args.get('calc_match') or self.args['subcmd'] == 'hull':
                # save special copy of calc_dict for hulls
                self.calc_dict = dict()
                self.calc_dict['$and'] = []
                # to avoid deep recursion, and since this is always called first
                # don't append, just set
                self.query_dict['$and'] = self.query_calc(self.cursor[0])
                self.calc_dict['$and'] = list(self.query_dict['$and'])
                if self.args['subcmd'] == 'hull':
                    self.args['composition'] = ''
                    for elem in self.cursor[0]['stoichiometry']:
                        self.args['composition'] += elem[0]
                    self.args['composition'] = [self.args['composition']]
                self.empty_query = False

        # create alias for formula for backwards-compatibility
        self.args['stoichiometry'] = self.args.get('formula')
        if self.args.get('stoichiometry') is not None:
            self.query_dict['$and'].append(self.query_stoichiometry())
            self.empty_query = False

        if self.args.get('composition') is not None:
            self.query_dict['$and'].append(self.query_composition())
            self.empty_query = False

        if self.args.get('num_species') is not None:
            self.query_dict['$and'].append(self.query_num_species())
            self.empty_query = False

        if self.args.get('pressure') is not None:
            self.query_dict['$and'].append(self.query_pressure())
            self.empty_query = False

        if self.args.get('space_group') is not None:
            self.query_dict['$and'].append(self.query_space_group())
            self.empty_query = False

        if self.args.get('num_fu') is not None:
            self.query_dict['$and'].append(self.query_num_fu())
            self.empty_query = False

        if self.args.get('encapsulated') is True:
            self.query_dict['$and'].append(self.query_encap())
            self.empty_query = False

        if self.args.get('cnt_radius') is not None:
            self.query_dict['$and'].append(self.query_cnt_radius())
            self.empty_query = False

        if self.args.get('tags') is not None:
            self.query_dict['$and'].append(self.query_tags())
            self.empty_query = False

        # only query quality when making a hull
        if not self.args.get('ignore_warnings'):
            self.query_dict['$and'].append(self.query_quality())

    def perform_query(self):
        """ Find results that match the query_dict
        inside the MongoDB database.
        """
        # if no query submitted, find all
        if self.empty_query:
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
        if not self.empty_query:
            self.cursors = []
            for collection in self.collections:
                self.repo = self.collections[collection]
                if self.args.get('details'):
                    print(dumps(self.query_dict, indent=4))

                # execute query
                self.cursor = self.repo.find(SON(self.query_dict)).sort('enthalpy_per_atom',
                                                                        pm.ASCENDING)
                self.cursors.append(self.cursor.clone())
                cursor_count = self.cursor.count()

                # write query to res or cell with param files
                if self.args.get('cell') or self.args.get('res'):
                    if cursor_count >= 1:
                        cursor = list(self.cursor)
                        if self.args.get('top') is not None:
                            query2files(cursor[:self.args.get('top')], self.args)
                        else:
                            query2files(cursor, self.args)

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
            if self.args.get('id') is None and (self.args.get('subcmd') == 'hull' or
                                                self.args.get('subcmd') == 'voltage' or
                                                self.args.get('hull_cutoff') is not None):
                if 'oqmd' in self.collections:
                    exit('Use --include_oqmd instead of --db, exiting...')
                if len(self.collections.keys()) == 1:
                    self.repo = self.collections[list(self.collections.keys())[0]]
                else:
                    exit('Hulls and voltage curves require just one source or --include_oqmd, \
                          exiting...')
                print('Creating hull from AJM db structures.')
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
                sample = 2
                rand_sample = 5 if self.args.get('biggest') else 5
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
                    if id_cursor.count() > 1:
                        print_warning('WARNING: matched multiple structures with text_id ' +
                                      id_cursor[0]['text_id'][0] + ' ' +
                                      id_cursor[0]['text_id'][1] + '.' +
                                      'Skipping this set...')
                        rand_sample += 1
                    else:
                        self.query_dict = dict()
                        self.query_dict['$and'] = self.query_calc(id_cursor[0])
                        cutoff.append(id_cursor[0]['cut_off_energy'])
                        calc_dicts.append(dict())
                        calc_dicts[-1]['$and'] = list(self.query_dict['$and'])
                        self.query_dict['$and'].append(self.query_composition())
                        if not self.args.get('ignore_warnings'):
                            self.query_dict['$and'].append(self.query_quality())
                        test_query_dict.append(self.query_dict)
                        test_cursors.append(
                            self.repo.find(SON(test_query_dict[-1])).sort('enthalpy_per_atom',
                                                                          pm.ASCENDING))
                        test_cursor_count.append(test_cursors[-1].count())
                        print("{:^24}".format(self.cursor[ind]['text_id'][0] + ' ' +
                                              self.cursor[ind]['text_id'][1]) +
                              ': matched ' + str(test_cursor_count[-1]), 'structures.', end=' -> ')
                        try:
                            print(self.cursor[ind]['xc_functional'] + ',',
                                  self.cursor[ind]['cut_off_energy'], 'eV,',
                                  self.cursor[ind]['kpoints_mp_spacing'], '1/A')
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

    def display_results(self, cursor, hull=False):
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
        if hull:
            header_string += "{:^18}".format('Hull dist./atom')
        else:
            header_string += "{:^18}".format('Enthalpy/fu')
        header_string += "{:^12}".format('Space group')
        header_string += "{:^10}".format('Formula')
        header_string += "{:^8}".format('# fu')
        header_string += "{:^8}".format('Prov.')

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
            if hull and doc['hull_distance'] == 0.0:
                struct_string.append(
                    '* ' + "{:^22}".format(doc['text_id'][0]+' '+doc['text_id'][1]))
            else:
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
                struct_string[-1] += "{:^12.3f}".format(doc['pressure'])
            except:
                struct_string[-1] += "{:^12}".format('xxx')
            try:
                struct_string[-1] += "{:^12.3f}".format(doc['cell_volume']/doc['num_fu'])
            except:
                struct_string[-1] += "{:^12}".format('xxx')
            try:
                if hull:
                    struct_string[-1] += "{:^18.5f}".format(doc['hull_distance'])
                else:
                    struct_string[-1] += "{:^18.5f}".format(doc['enthalpy']/doc['num_fu'] -
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
            try:
                source = 'AIRSS'
                if doc['source'] is str:
                    source_list = [doc['source']]
                else:
                    source_list = doc['source']
                for fname in source_list:
                    if (fname.endswith('.castep') or fname.endswith('.res') or
                            fname.endswith('.history') or 'OQMD' in fname):
                        if 'swap' in fname.lower():
                            source = 'SWAPS'
                        elif 'oqmd' in fname.lower():
                            source = 'OQMD'
                        elif 'collcode' in fname.lower():
                            if fname.split('/')[-1].count('-') == 2:
                                source = 'SWAPS'
                            else:
                                source = 'ICSD'
                        elif '-icsd' in fname.lower():
                            source = 'ICSD'
                struct_string[-1] += "{:^8}".format(source)
            except:
                struct_string[-1] += "{:^8}".format('xxx')

            if last_formula != formula_substring:
                self.gs_enthalpy = doc['enthalpy'] / doc['num_fu']
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
            tmp_stoich = stoich
            for ind, strng in enumerate(stoich):
                tmp_stoich[ind] = [elem for elem in re.split(r'([0-9]*)', strng) if elem]
            stoich = [item for sublist in tmp_stoich for item in sublist]
            while '[' in stoich or '][' in stoich:
                tmp_stoich = list(stoich)
                for ind, tmp in enumerate(tmp_stoich):
                    if tmp == '][':
                        del tmp_stoich[ind]
                        tmp_stoich.insert(ind, '[')
                        tmp_stoich.insert(ind, ']')
                        break
                for ind, tmp in enumerate(tmp_stoich):
                    if tmp == '[':
                        end_bracket = False
                        while not end_bracket:
                            if tmp_stoich[ind+1] == ']':
                                end_bracket = True
                            tmp_stoich[ind] += tmp_stoich[ind+1]
                            del tmp_stoich[ind+1]
                try:
                    tmp_stoich.remove(']')
                except:
                    pass
                try:
                    tmp_stoich.remove('')
                except:
                    pass
                stoich = tmp_stoich

        elements = []
        fraction = []
        for i in range(0, len(stoich), 1):
            if not bool(re.search(r'\d', stoich[i])):
                elements.append(stoich[i])
                try:
                    fraction.append(float(stoich[i+1]))
                except:
                    fraction.append(1.0)
        gcd_val = 0
        for frac in fraction:
            if gcd_val == 0:
                gcd_val = frac
            else:
                gcd_val = gcd(frac, gcd_val)
        fraction = np.asarray(fraction)
        fraction /= gcd_val

        query_dict = dict()
        query_dict['$and'] = []

        for ind, elem in enumerate(elements):
            if '[' in elem or ']' in elem:
                types_dict = dict()
                types_dict['$or'] = list()
                elem = elem.strip('[').strip(']')
                for group_elem in self.periodic_table[elem]:
                    types_dict['$or'].append(dict())
                    types_dict['$or'][-1]['stoichiometry'] = dict()
                    types_dict['$or'][-1]['stoichiometry']['$in'] = [[group_elem, fraction[ind]]]
                query_dict['$and'].append(types_dict)
            else:
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
        for char in elements[0]:
            if char.isdigit():
                print_failure('Composition cannot contain a number.')
                exit()
        try:
            if len(elements) == 1:
                valid = False
                for char in elements[0]:
                    if char.isupper():
                        valid = True
                if not valid:
                    print_failure('Composition must contain at least one upper case character.')
                    exit()
                elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements[0]) if elem]
                while '[' in elements or '][' in elements:
                    tmp_stoich = list(elements)
                    for ind, tmp in enumerate(tmp_stoich):
                        if tmp == '][':
                            del tmp_stoich[ind]
                            tmp_stoich.insert(ind, '[')
                            tmp_stoich.insert(ind, ']')
                            break
                    for ind, tmp in enumerate(tmp_stoich):
                        if tmp == '[':
                            end_bracket = False
                            while not end_bracket:
                                if tmp_stoich[ind+1] == ']':
                                    end_bracket = True
                                tmp_stoich[ind] += tmp_stoich[ind+1]
                                del tmp_stoich[ind+1]
                    try:
                        tmp_stoich.remove(']')
                    except:
                        pass
                    try:
                        tmp_stoich.remove('')
                    except:
                        pass
                    elements = tmp_stoich
        except Exception:
            print_exc()
            exit()

        if self.args.get('intersection'):
            query_dict = dict()
            query_dict['$or'] = []
            # iterate over all combinations
            for rlen in range(1, len(elements)+1):
                for combi in combinations(elements, r=rlen):
                    list_combi = list(combi)
                    types_dict = dict()
                    types_dict['$and'] = list()
                    types_dict['$and'].append(dict())
                    types_dict['$and'][-1]['stoichiometry'] = dict()
                    types_dict['$and'][-1]['stoichiometry']['$size'] = len(list_combi)
                    for elem in list_combi:
                        types_dict['$and'].append(dict())
                        types_dict['$and'][-1]['atom_types'] = dict()
                        types_dict['$and'][-1]['atom_types']['$in'] = [elem]
                    query_dict['$or'].append(types_dict)
        else:
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
                        types_dict['$or'][-1]['atom_types'] = dict()
                        types_dict['$or'][-1]['atom_types']['$in'] = [group_elem]
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

    def query_cnt_radius(self):
        """ Query structures within a nanotube of given radius
        to within a tolerance of 0.01 A.
        """
        query_dict = dict()
        query_dict['$and'] = []
        query_dict['$and'].append(dict())
        query_dict['$and'][-1]['cnt_radius'] = dict()
        query_dict['$and'][-1]['cnt_radius']['$gt'] = self.args.get('cnt_radius') - 0.01
        query_dict['$and'][-1]['cnt_radius']['$lt'] = self.args.get('cnt_radius') + 0.01

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
        if 'spin_polarized' in doc and doc['spin_polarized']:
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
        else:
            temp_dict = dict()
            temp_dict['kpoints_mp_spacing'] = dict()
            temp_dict['kpoints_mp_spacing']['$gte'] = doc['kpoints_mp_spacing'] - 0.02
            temp_dict['kpoints_mp_spacing']['$lte'] = doc['kpoints_mp_spacing'] + 0.01
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
