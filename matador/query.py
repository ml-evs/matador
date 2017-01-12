# coding: utf-8
""" This file implements all queries to the database,
including parsing user inputs, displaying results
and calling other functionality. """

from __future__ import print_function
# matador modules
from .utils.print_utils import print_failure, print_warning, print_success
from .utils.chem_utils import get_periodic_table
from .utils.cursor_utils import display_results
# external libraries
import pymongo as pm
import numpy as np
from bson.son import SON
from bson.json_util import dumps
# standard library
import re
from os import uname
from itertools import combinations
from traceback import print_exc
try:
    from math import gcd
except ImportError:
    from fractions import gcd
    print_warning('Use of fractions.gcd is deprecated - \
                   Python3 is recommended but will try to proceed.')
from sys import exit


class DBQuery:
    """ Class that implements queries to MongoDB
    structure database.
    """
    def __init__(self, client=False, collections=False, subcmd='query', debug=False, **kwargs):
        """ Parse arguments from matador or API call
        before calling query.
        """
        # read args
        self.args = kwargs
        self.debug = debug
        if debug:
            print(self.args)
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd
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
        # improve this clause at some point
        if self.args.get('summary') or self.args.get('subcmd') in ['swaps', 'polish']:
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

        # operate on one structure and related others
        if self.args.get('id') is not None:
            if type(self.args.get('id')) == str:
                self.args['id'] = self.args['id'].strip().split(' ')
            self.cursor = []

            for collection in self.collections:
                query_dict = dict()
                query_dict['$and'] = []
                query_dict['$and'].append(self.query_id())
                if not self.args.get('ignore_warnings'):
                    query_dict['$and'].append(self.query_quality())
                self.repo = self.collections[collection]
                temp_cursor = self.repo.find(query_dict)
                for doc in temp_cursor:
                    self.cursor.append(doc)

            if len(self.cursor) < 1:
                exit('Could not find a match with ' + str(self.args.get('id')) + ' try widening your search.')
            elif len(self.cursor) >= 1:
                display_results(list(self.cursor)[:self.top], args=self.args)
                if len(self.cursor) > 1:
                    print_warning('WARNING: matched multiple structures with same text_id. ' +
                                  'The first one will be used.')
                if self.debug:
                    print(dumps(self.cursor[0], indent=1))

            if self.args.get('calc_match') or self.args['subcmd'] == 'hull':
                # save special copy of calc_dict for hulls
                self.calc_dict = dict()
                self.calc_dict['$and'] = []
                # to avoid deep recursion, and since this is always called first
                # don't append, just set
                self.query_dict['$and'] = self.query_calc(self.cursor[0])
                self.calc_dict['$and'] = list(self.query_dict['$and'])
                if self.args['subcmd'] == 'hull' and self.args.get('composition') is None:
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

        if self.args.get('cutoff') is not None:
            self.query_dict['$and'].append(self.query_cutoff())
            self.empty_query = False

        if self.args.get('sedc') is not None:
            self.query_dict['$and'].append(self.query_sedc())
            self.empty_query = False

        if self.args.get('mp_spacing') is not None:
            self.query_dict['$and'].append(self.query_kpoints())
            self.empty_query = False

        if self.args.get('spin') is not None:
            self.query_dict['$and'].append(self.query_spin())
            self.empty_query = False

        if self.args.get('tags') is not None:
            self.query_dict['$and'].append(self.query_tags())
            self.empty_query = False

        if self.args.get('doi') is not None:
            self.query_dict['$and'].append(self.query_doi())
            self.empty_query = False

        if self.args.get('src_str') is not None:
            self.query_dict['$and'].append(self.query_source())
            self.empty_query = False

        if self.args.get('icsd') is not None:
            self.query_dict['$and'].append(self.query_icsd())
            self.empty_query = False

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
                    if self.debug:
                        print('Empty query, showing all...')
                    self.cursor = self.repo.find().sort('enthalpy_per_atom', pm.ASCENDING)
                    if self.top == -1:
                        self.top = len(self.cursor)
                    display_results(list(self.cursor[:self.top]), args=self.args)

        # if no special query has been made already, begin executing the query
        if not self.empty_query:
            # self.cursors = []
            for collection in self.collections:
                self.repo = self.collections[collection]
                if self.debug:
                    print(dumps(self.query_dict, indent=1))
                # execute query
                self.cursor = list(self.repo.find(SON(self.query_dict)).sort('enthalpy_per_atom',
                                                                             pm.ASCENDING))
                # self.cursors.append(self.cursor)
                cursor_count = len(self.cursor)

                # if called as script, always print results
                if self.args.get('id') is None:
                    print(cursor_count, 'results found for query in', collection+'.')
                if self.args.get('subcmd') != 'hull' and self.args.get('subcmd') != 'voltage':
                    if cursor_count >= 1:
                        if self.top == -1:
                            self.top = cursor_count
                        if cursor_count > self.top:
                            display_results(list(self.cursor)[:self.top], args=self.args)
                        else:
                            display_results(list(self.cursor), args=self.args)

            # building hull from just comp, find best structure to calc_match
            if self.args.get('id') is None and (self.args.get('subcmd') == 'hull' or
                                                self.args.get('subcmd') == 'voltage' or
                                                self.args.get('hull_cutoff') is not None):
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
                rand_sample = 5 if self.args.get('biggest') else 3
                i = 0
                count = len(self.cursor)
                if count <= 0:
                    exit('No structures found for hull.')
                while i < sample+rand_sample:
                    # start with sample/2 lowest enthalpy structures
                    if i < int(sample):
                        ind = i
                    # then do some random samples
                    else:
                        ind = np.random.randint(rand_sample if rand_sample < count-1 else 0, count-1)
                    id_cursor = list(self.repo.find({'text_id': self.cursor[ind]['text_id']}))
                    if len(id_cursor) > 1:
                        print_warning('WARNING: matched multiple structures with text_id ' +
                                      id_cursor[0]['text_id'][0] + ' ' +
                                      id_cursor[0]['text_id'][1] + '.' +
                                      'Skipping this set...')
                        rand_sample += 1
                    else:
                        self.query_dict = dict()
                        try:
                            self.query_dict['$and'] = self.query_calc(id_cursor[0])
                            cutoff.append(id_cursor[0]['cut_off_energy'])
                            calc_dicts.append(dict())
                            calc_dicts[-1]['$and'] = list(self.query_dict['$and'])
                            self.query_dict['$and'].append(self.query_composition())
                            if not self.args.get('ignore_warnings'):
                                self.query_dict['$and'].append(self.query_quality())
                            test_query_dict.append(self.query_dict)
                            test_cursors.append(
                                list(self.repo.find(SON(test_query_dict[-1])).sort('enthalpy_per_atom',
                                                                                   pm.ASCENDING)))
                            test_cursor_count.append(len(test_cursors[-1]))
                            print("{:^24}".format(self.cursor[ind]['text_id'][0] + ' ' +
                                                  self.cursor[ind]['text_id'][1]) +
                                  ': matched ' + str(test_cursor_count[-1]), 'structures.', end='\t-> ')
                            print('S-' if self.cursor[ind].get('spin_polarized') else '',
                                  self.cursor[ind]['sedc_scheme'] + '-' if self.cursor[ind].get('sedc_scheme') is not None else '',
                                  self.cursor[ind]['xc_functional'] + ', ',
                                  self.cursor[ind]['cut_off_energy'], ' eV, ',
                                  self.cursor[ind]['kpoints_mp_spacing'] if self.cursor[ind].get('kpoints_mp_spacing') is not None else 'xxx', ' 1/A', sep='')
                            if test_cursor_count[-1] == count:
                                print('Matched all structures...')
                                break
                            if self.args.get('biggest'):
                                if test_cursor_count[-1] > 2*int(count/3):
                                    print('Matched at least 2/3 of total number, composing hull...')
                                    break
                        except(KeyboardInterrupt, SystemExit):
                            print('Received exit signal, exiting...')
                            exit()
                        except:
                            print_exc()
                            print_warning('Error with ' + id_cursor[0]['text_id'][0] + ' ' + id_cursor[0]['text_id'][1])
                            rand_sample += 1
                    i += 1

                if self.args.get('biggest'):
                    choice = np.argmax(np.asarray(test_cursor_count))
                else:
                    # by default, find highest cutoff hull as first proxy for quality
                    choice = np.argmax(np.asarray(cutoff))
                self.cursor = test_cursors[choice]
                print_success('Composing hull from set containing ' +
                              self.cursor[0]['text_id'][0] + ' ' + self.cursor[0]['text_id'][1])
                self.calc_dict = calc_dicts[choice]

    def __del__(self):
        """ Clean up any temporary databases on garbage
        collection of DBQuery object.
        """
        try:
            self.temp.drop()
        except:
            pass

    def query_stoichiometry(self, custom_stoich=None, partial_formula=None):
        """ Query DB for particular stoichiometry. """
        # alias stoichiometry
        if custom_stoich is None:
            stoich = self.args.get('stoichiometry')
        else:
            stoich = custom_stoich
        if partial_formula is None:
            partial_formula = self.args.get('partial_formula')

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
                gcd_val = gcd(int(frac), int(gcd_val))
        fraction = np.asarray(fraction)
        fraction /= gcd_val

        query_dict = dict()
        query_dict['$and'] = []

        for ind, elem in enumerate(elements):
            if '[' in elem or ']' in elem:
                types_dict = dict()
                types_dict['$or'] = list()
                elem = elem.strip('[').strip(']')
                if elem in self.periodic_table:
                    for group_elem in self.periodic_table[elem]:
                        types_dict['$or'].append(dict())
                        types_dict['$or'][-1]['stoichiometry'] = dict()
                        types_dict['$or'][-1]['stoichiometry']['$in'] = [[group_elem, fraction[ind]]]
                    query_dict['$and'].append(types_dict)
                elif ',' in elem:
                    for group_elem in elem.split(','):
                        types_dict['$or'].append(dict())
                        types_dict['$or'][-1]['stoichiometry'] = dict()
                        types_dict['$or'][-1]['stoichiometry']['$in'] = [[group_elem, fraction[ind]]]
                    query_dict['$and'].append(types_dict)
            else:
                stoich_dict = dict()
                stoich_dict['stoichiometry'] = dict()
                stoich_dict['stoichiometry']['$in'] = [[elem, fraction[ind]]]
                query_dict['$and'].append(stoich_dict)
        if not partial_formula:
            size_dict = dict()
            size_dict['stoichiometry'] = dict()
            size_dict['stoichiometry']['$size'] = len(elements)
            query_dict['$and'].append(size_dict)

        return query_dict

    def query_ratio(self, ratios):
        """ Query DB for ratio of two elements.

        Input, e.g.:

            ratios = [['MoS', 2],
                      ['LiS', 1]]

        """
        query_dict = dict()
        for pair in ratios:
            query_dict['ratios.' + pair[0]] = pair[1]
        return query_dict

    def query_composition(self, custom_elem=None, partial_formula=None):
        """ Query DB for all structures containing
        all the elements taken as input. Passing this
        function a number is a deprecated feature, replaced
        by query_num_species.
        """
        if custom_elem is None:
            elements = list(self.args.get('composition'))
        else:
            elements = custom_elem
        if partial_formula is None:
            partial_formula = self.args.get('partial_formula')
        non_binary = False
        if ':' in elements[0]:
            non_binary = True
        # if there's only one string, try split it by caps
        if not non_binary:
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
            size = len(elements)
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
        elif non_binary:
            query_dict = dict()
            query_dict['$and'] = []
            size = 0
            for ind, elem in enumerate(elements):
                if elem != ':' and not elem.isdigit():
                    query_dict['$and'].append(self.query_composition(custom_elem=[elem], partial_formula=True))
                    size += 1
                if elem == ':':
                    # convert e.g. MoS2 to [['MoS', 2]]
                    # or LiMoS2 to [['LiMo', 1], ['MoS', '2], ['LiS', 2]]
                    ratio_elements = elements[ind+1:]
                    for ind in range(len(ratio_elements)):
                        if ind < len(ratio_elements)-1:
                            if not ratio_elements[ind].isdigit() and not ratio_elements[ind+1].isdigit():
                                ratio_elements.insert(ind+1, '1')
                    if not ratio_elements[-1].isdigit():
                        ratio_elements.append('1')
                    ratios = []
                    for ind in range(0, len(ratio_elements), 2):
                        for jind in range(ind, len(ratio_elements), 2):
                            if ratio_elements[ind] != ratio_elements[jind]:
                                ratios.append([ratio_elements[ind]+ratio_elements[jind],
                                               round(float(ratio_elements[ind+1])/float(ratio_elements[jind+1]), 3)])
                    query_dict['$and'].append(self.query_ratio(ratios))
        else:
            query_dict = dict()
            query_dict['$and'] = []
            size = len(elements)
            for ind, elem in enumerate(elements):
                # prototype for chemically motivated searches, e.g. transition metals
                if '[' in elem or ']' in elem:
                    types_dict = dict()
                    types_dict['$or'] = list()
                    elem = elem.strip('[').strip(']')
                    if elem in self.periodic_table:
                        for group_elem in self.periodic_table[elem]:
                            types_dict['$or'].append(dict())
                            types_dict['$or'][-1]['atom_types'] = dict()
                            types_dict['$or'][-1]['atom_types']['$in'] = [group_elem]
                    elif ',' in elem:
                        for group_elem in elem.split(','):
                            types_dict['$or'].append(dict())
                            types_dict['$or'][-1]['atom_types'] = dict()
                            types_dict['$or'][-1]['atom_types']['$in'] = [group_elem]
                else:
                    types_dict = dict()
                    types_dict['atom_types'] = dict()
                    types_dict['atom_types']['$in'] = [elem]
                query_dict['$and'].append(types_dict)
        if not partial_formula and not self.args.get('intersection'):
            size_dict = dict()
            size_dict['stoichiometry'] = dict()
            size_dict['stoichiometry']['$size'] = size
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

    def query_doi(self):
        """ Find all structures matching given DOI,
        in format xxxx/xxxx.
        """
        query_dict = dict()
        query_dict['$and'] = []
        doi = self.args.get('doi')
        temp_dict = dict()
        temp_dict['doi'] = dict()
        temp_dict['doi']['$eq'] = doi
        query_dict['$and'].append(temp_dict)

        return query_dict

    def query_id(self):
        """ Find all structures matching given tags. """
        query_dict = dict()
        query_dict['text_id'] = self.args.get('id')
        return query_dict

    def query_icsd(self):
        """ Find all structures matching given ICSD CollCode. """
        query_dict = dict()
        query_dict['$and'] = []
        sub_dict = dict()
        sub_dict['icsd'] = dict()
        if self.args.get('icsd') == 0:
            sub_dict['icsd']['$exists'] = True
        else:
            sub_dict['icsd']['$eq'] = str(self.args.get('icsd'))
        query_dict['$and'].append(sub_dict)
        return query_dict

    def query_source(self):
        """ Find all structures with partial source string from args. """
        query_dict = dict()
        query_dict['$and'] = []
        sub_dict = dict()
        sub_dict['source'] = dict()
        sub_dict['source']['$in'] = [self.args.get('src_str')]
        query_dict['$and'].append(sub_dict)
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
        temp_dict['pressure'] = dict()
        temp_dict['pressure']['$lt'] = approx_pressure[1]
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

    def query_cutoff(self):
        """ Query all calculations above given plane-wave cutoff. """
        query_dict = dict()
        query_dict['cut_off_energy'] = dict()
        if len(self.args.get('cutoff')) == 2:
            query_dict['cut_off_energy']['$gte'] = self.args.get('cutoff')[0]
            query_dict['cut_off_energy']['$lte'] = self.args.get('cutoff')[1]
        else:
            query_dict['cut_off_energy']['$gte'] = self.args.get('cutoff')[0]
        return query_dict

    def query_sedc(self):
        """ Query all calculations using given SEDC scheme.

        Use --sedc null to query for no dispersion correction.
        """
        query_dict = dict()
        if self.args.get('sedc') != 'null':
            query_dict['sedc_scheme'] = self.args.get('sedc')
        else:
            query_dict['sedc_scheme'] = dict()
            query_dict['sedc_scheme']['$exists'] = False

        return query_dict

    def query_kpoints(self):
        """ Query all calculations with finer than the given
        kpoint sampling.
        """
        query_dict = dict()
        query_dict['kpoints_mp_spacing'] = dict()
        query_dict['kpoints_mp_spacing']['$lte'] = self.args.get('mp_spacing')
        return query_dict

    def query_spin(self):
        """ Query all calculations with spin polarisation,
        i.e. --spin n!=0, or non-spin-polarization, i.e. --spin 0.
        """
        query_dict = dict()
        if self.args.get('spin') == '0':
            query_dict['spin_polarized'] = dict()
            query_dict['spin_polarized']['$ne'] = True
        else:
            query_dict['spin_polarized'] = True
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
        if 'sedc_scheme' in doc:
            temp_dict['sedc_scheme'] = doc['sedc_scheme']
            query_dict.append(temp_dict)
        else:
            temp_dict['sedc_scheme'] = dict()
            temp_dict['sedc_scheme']['$exists'] = False
            query_dict.append(temp_dict)

        db = self.args.get('db')
        if db is not None:
            db = db[0]
        else:
            db = ''
        if self.args.get('loose') or 'oqmd' in db:
            return query_dict
            # temp_dict = dict()
            # query_dict.append(dict())
            # query_dict[-1]['cut_off_energy'] = doc['cut_off_energy']
        else:
            temp_dict = dict()
            temp_dict['kpoints_mp_spacing'] = dict()
            if self.args.get('kpoint_tolerance') is not None:
                try:
                    tol = float(self.args.get('kpoint_tolerance'))
                except:
                    print_warning('Failed to read custom kpoint tolerance.')
            else:
                tol = 0.01
            temp_dict['kpoints_mp_spacing']['$gte'] = doc['kpoints_mp_spacing'] - tol
            temp_dict['kpoints_mp_spacing']['$lte'] = doc['kpoints_mp_spacing'] + tol
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
        if len(cursor) != 0:
            self.temp.insert(cursor)
        else:
            self.temp.drop()
            exit('No structures found.')

        return self.temp


class EmptyCursor:
    """ Empty cursor class for failures. """
    def count(self):
        return 0
