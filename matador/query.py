# coding: utf-8
""" This file implements all queries to the database,
including parsing user inputs, displaying results
and calling other functionality.
"""

from os import devnull
import sys
from itertools import combinations
from traceback import print_exc

import pymongo as pm
import numpy as np
from bson.son import SON
from bson.json_util import dumps
from bson.objectid import ObjectId

from matador.utils.print_utils import print_failure, print_warning, print_success, print_notify
from matador.utils.chem_utils import get_periodic_table, get_formula_from_stoich
from matador.utils.chem_utils import parse_element_string, get_stoich_from_formula
from matador.utils.cursor_utils import display_results
from matador.utils.db_utils import make_connection_to_collection, load_custom_settings


class DBQuery:
    """ Class that implements queries to MongoDB
    structure database.

    Attributes:
        cursor (list of dict): list of structures matching query.
        args (dict): contains all keyword arguments used to construct
            the query (see matador query --help) for list.
        query_dict (dict): dictionary passed to database for query
        calc_dict (dict): if performing a matching query (e.g. self.args['subcmd'] = 'hull'), this
            dictionary contains the parameters used to match to other structures
        repo (pymongo.collection.Collection): the pymongo collection that is being queried.
        top (int): number of structures to print/export set by self.args.get('top') (DEFAULT: 10).

    """
    def __init__(self, client=False, collections=False, subcmd='query', debug=False, quiet=False, **kwargs):
        """ Parse arguments from matador or API call
        before calling query.
        """
        # read args and set housekeeping
        self.args = kwargs
        self.debug = debug
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd
        if self.args.get('testing') is None:
            self.args['testing'] = False

        # public attributes
        self.cursor = EmptyCursor()
        self.query_dict = None
        self.calc_dict = None
        self.repo = None

        # private attributes to be set later
        self._empty_query = None
        self._gs_enthalpy = None

        if debug:
            print(self.args)

        if quiet:
            f = open(devnull, 'w')
            sys.stdout = f

        # connect to db or use passed client
        if client:
            self._client = client
            self._db = client.crystals
        if collections is not False:
            self._collections = collections

        if (not collections or not client) and not self.args.get('testing'):
            self.mongo_settings = load_custom_settings(config_fname=self.args.get('config_fname'))
            result = make_connection_to_collection(self.args.get('db'),
                                                   mongo_settings=self.mongo_settings)
            self._client, self._db, self._collections = result

        # define some periodic table macros
        self._periodic_table = get_periodic_table()

        # set default top value to 10
        if self.args.get('summary') or self.args.get('subcmd') in ['swaps', 'polish']:
            self.top = None
        else:
            self.top = self.args.get('top') if self.args.get('top') is not None else 10

        # create the dictionary to pass to MongoDB
        self._construct_query()

        if not self.args.get('testing'):

            # execute the query
            self.perform_query()

            # only filter for uniqueness if not eventually making a hull
            if self.args.get('uniq') and self.args.get('subcmd') not in ['hull', 'hulldiff', 'voltage']:
                from matador.similarity.similarity import get_uniq_cursor
                print_notify('Filtering for unique structures...')
                if self.top is not None:
                    # filter for uniqueness
                    unique_set, _, _, _ = get_uniq_cursor(self.cursor[:self.args.get('top')],
                                                          debug=self.args.get('debug'),
                                                          sim_tol=self.args.get('uniq'))
                    print('Filtered {} down to {}'.format(len(self.cursor[:self.args.get('top')]),
                                                          len(unique_set)))
                    self.cursor = [self.cursor[:self.args.get('top')][ind] for ind in unique_set]
                else:
                    unique_set, _, _, _ = get_uniq_cursor(self.cursor,
                                                          debug=self.args.get('debug'),
                                                          sim_tol=self.args.get('uniq'))
                    print('Filtered {} down to {}'.format(len(self.cursor),
                                                          len(unique_set)))
                    self.cursor = [self.cursor[ind] for ind in unique_set]

                display_results(self.cursor, hull=None, args=self.args)

            if not client and not self.args.get('testing'):
                self._client.close()

        if quiet:
            f.close()
            sys.stdout = sys.__stdout__

    def _construct_query(self):
        """ Set up query dict and perform query depending on
        command-line / API arguments.
        """

        self.cursor = EmptyCursor()
        # initalize query_dict to '$and' all queries
        self.query_dict = dict()
        self.query_dict['$and'] = []
        self._empty_query = True

        # benchmark enthalpy to display (set by calc_match)
        self._gs_enthalpy = 0.0

        # operate on one structure and related others
        if self.args.get('id') is not None:
            if isinstance(self.args.get('id'), str):
                self.args['id'] = self.args['id'].strip().split(' ')

            for collection in self._collections:
                query_dict = dict()
                query_dict['$and'] = []
                query_dict['$and'].append(self._query_id())
                if not self.args.get('ignore_warnings'):
                    query_dict['$and'].append(self._query_quality())
                self.repo = self._collections[collection]
                temp_cursor = self.repo.find(query_dict)
                for doc in temp_cursor:
                    self.cursor.append(doc)

            if len(self.cursor) < 1:
                sys.exit('Could not find a match with {} try widening your search.'
                         .format(self.args.get('id')))

            elif len(self.cursor) >= 1:
                display_results(list(self.cursor)[:self.top], args=self.args)
                if len(self.cursor) > 1:
                    print_warning('WARNING: matched multiple structures with same text_id. ' +
                                  'The first one will be used.')
                if self.debug:
                    print(dumps(self.cursor[0], indent=1))

            if self.args.get('calc_match') or \
                    self.args['subcmd'] in ['hull', 'hulldiff', 'voltage']:
                # save special copy of calc_dict for hulls
                self.calc_dict = dict()
                self.calc_dict['$and'] = []
                # to avoid deep recursion, and since this is always called first
                # don't append, just set
                self.query_dict = self._query_calc(self.cursor[0])
                self.calc_dict['$and'] = list(self.query_dict['$and'])
                if self.args['subcmd'] in ['hull', 'hulldiff'] and \
                        self.args.get('composition') is None:
                    self.args['composition'] = ''
                    for elem in self.cursor[0]['stoichiometry']:
                        self.args['composition'] += elem[0]
                    self.args['composition'] = [self.args['composition']]
                self._empty_query = False

        # create alias for formula for backwards-compatibility
        self.args['stoichiometry'] = self.args.get('formula')
        if self.args.get('stoichiometry') is not None:
            self.query_dict['$and'].append(self._query_stoichiometry())
            self._empty_query = False

        if self.args.get('composition') is not None:
            self.query_dict['$and'].append(self._query_composition())
            self._empty_query = False

        if self.args.get('num_species') is not None:
            self.query_dict['$and'].append(self._query_num_species())
            self._empty_query = False

        if self.args.get('space_group') is not None:
            self.query_dict['$and'].append(self._query_space_group())
            self._empty_query = False

        if self.args.get('num_fu') is not None:
            self.query_dict['$and'].append(self._query_num_fu())
            self._empty_query = False

        if self.args.get('tags') is not None:
            self.query_dict['$and'].append(self._query_tags())
            self._empty_query = False

        if self.args.get('doi') is not None:
            self.query_dict['$and'].append(self._query_doi())
            self._empty_query = False

        if self.args.get('icsd') is not None:
            self.query_dict['$and'].append(self._query_icsd())
            self._empty_query = False

        if self.args.get('cutoff') is not None:
            self.query_dict['$and'].append(self._query_cutoff())
            self._empty_query = False

        if self.args.get('geom_force_tol') is not None:
            self.query_dict['$and'].append(self._query_geom_force_tol())
            self._empty_query = False

        if self.args.get('src_str') is not None:
            self.query_dict['$and'].append(self._query_source())
            self._empty_query = False

        if self.args.get('root_src') is not None:
            self.query_dict['$and'].append(self._query_root_source())
            self._empty_query = False

        if self.args.get('pressure') is not None:
            self.query_dict['$and'].append(self._query_pressure())
            self._empty_query = False

        if self.args.get('encapsulated') is True:
            self.query_dict['$and'].append(self._query_encap())
            self._empty_query = False

        if self.args.get('cnt_radius') is not None:
            self.query_dict['$and'].append(self._query_cnt_radius())
            self._empty_query = False

        if self.args.get('cnt_vector') is not None:
            self.query_dict['$and'].append(self._query_cnt_vector())
            self._empty_query = False

        if self.args.get('sedc') is not None:
            self.query_dict['$and'].append(self._query_sedc())
            self._empty_query = False

        if self.args.get('xc_functional') is not None:
            self.query_dict['$and'].append(self._query_xc_functional())
            self._empty_query = False

        if self.args.get('mp_spacing') is not None:
            self.query_dict['$and'].append(self._query_kpoints())
            self._empty_query = False

        if self.args.get('spin') is not None:
            tmp_dict = self._query_spin()
            if tmp_dict:
                self.query_dict['$and'].append(tmp_dict)
            self._empty_query = False

        if not self.args.get('ignore_warnings'):
            self.query_dict['$and'].append(self._query_quality())

        if self.args.get('time') is not None:
            self.query_dict['$and'].append(self._query_time())
            self._empty_query = False

    def perform_query(self):
        """ Find results that match the query_dict
        inside the MongoDB database.
        """
        # if no query submitted, find all
        if self._empty_query:
            if self.args.get('id') is None:
                for collection in self._collections:
                    self.repo = self._collections[collection]
                    if self.debug:
                        print('Empty query, showing all...')
                    self.cursor = self.repo.find().sort('enthalpy_per_atom', pm.ASCENDING)
                    if self.top == -1 or self.top is None:
                        self.top = self.cursor.count()
                    self.cursor = list(self.cursor)
                    display_results(self.cursor[:self.top], args=self.args)

        # if no special query has been made already, begin executing the query
        if not self._empty_query:
            for collection in self._collections:
                self.repo = self._collections[collection]
                if self.debug:
                    print(dumps(self.query_dict, indent=1))
                # execute query
                self.cursor = list(self.repo.find(SON(self.query_dict))
                                   .sort('enthalpy_per_atom', pm.ASCENDING))

                # self.cursors.append(self.cursor)
                cursor_count = len(self.cursor)

                # if called as script, always print results
                if self.args.get('id') is None:
                    print(cursor_count, 'results found for query in', collection+'.')
                if self.args.get('subcmd') not in ['hull', 'hulldiff', 'voltage', 'swaps']:
                    if cursor_count >= 1:
                        self._num_to_display = cursor_count
                        if self.args.get('delta_E') is not None:
                            self.cursor = list(self.cursor)
                            if len(set([get_formula_from_stoich(doc['stoichiometry']) for doc in self.cursor])) != 1:
                                print('Multiple stoichiometries in cursor, unable to filter by energy.')
                            else:
                                gs_enthalpy = self.cursor[0]['enthalpy_per_atom']
                                if self.debug:
                                    print('Filtering by {} eV/atom'.format(self.args.get('delta_E')))
                                    print('gs_enthalpy = {}'.format(gs_enthalpy))
                                num_to_display = 1
                                for doc in self.cursor[1:]:
                                    if doc['enthalpy_per_atom'] - gs_enthalpy > self.args.get('delta_E'):
                                        break
                                    else:
                                        num_to_display += 1
                                self._num_to_display = num_to_display
                                cursor_count = self._num_to_display
                                if self.debug:
                                    print('Displaying top {}'.format(self._num_to_display))
                        elif self.top == -1 or self.top is None:
                            self._num_to_display = cursor_count
                            self.top = cursor_count
                        elif cursor_count > self.top:
                            self._num_to_display = self.top

                        display_results(list(self.cursor)[:self._num_to_display], args=self.args)

                if self.args.get('delta_E') is not None:
                    self.cursor = self.cursor[:self._num_to_display]

            # building hull from just comp, find best structure to calc_match
            if self.args.get('id') is None and (self.args.get('subcmd') == 'hull' or
                                                self.args.get('subcmd') == 'hulldiff' or
                                                self.args.get('subcmd') == 'voltage' or
                                                self.args.get('hull_cutoff') is not None):
                if len(self._collections) == 1:
                    self.repo = self._collections[list(self._collections.keys())[0]]
                else:
                    sys.exit('Hulls and voltage curves require just one source or --include_oqmd, exiting...')
                print('Creating hull from AJM db structures.')
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
                    sys.exit('No structures found for hull.')
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
                                      ' Skipping this set...')
                        rand_sample += 1
                    else:
                        self.query_dict = dict()
                        try:
                            self.query_dict = self._query_calc(id_cursor[0])
                            cutoff.append(id_cursor[0]['cut_off_energy'])
                            calc_dicts.append(dict())
                            calc_dicts[-1]['$and'] = list(self.query_dict['$and'])
                            self.query_dict['$and'].append(self._query_composition())
                            if not self.args.get('ignore_warnings'):
                                self.query_dict['$and'].append(self._query_quality())
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
                                  self.cursor[ind]['geom_force_tol'] if self.cursor[ind].get('geom_force_tol') is not None else 'xxx', ' eV/A, ',
                                  self.cursor[ind]['kpoints_mp_spacing'] if self.cursor[ind].get('kpoints_mp_spacing') is not None else 'xxx', ' 1/A', sep='')
                            if test_cursor_count[-1] == count:
                                print('Matched all structures...')
                                break
                            if test_cursor_count[-1] > 2*int(count/3):
                                print('Matched at least 2/3 of total number, composing hull...')
                                break
                        except(KeyboardInterrupt, SystemExit):
                            print('Received exit signal, exiting...')
                            sys.exit()
                        except Exception:
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

    def _query_stoichiometry(self, custom_stoich=None, partial_formula=None):
        """ Query DB for particular stoichiometry. """
        # alias stoichiometry
        if custom_stoich is None:
            stoich = self.args.get('stoichiometry')
            if isinstance(stoich, str):
                stoich = [stoich]
        else:
            stoich = custom_stoich
        if partial_formula is None:
            partial_formula = self.args.get('partial_formula')
        if ':' in stoich[0]:
            sys.exit('Formula cannot contain ":", you probably meant to query composition.')

        stoich = get_stoich_from_formula(stoich[0])

        query_dict = dict()
        query_dict['$and'] = []

        for ind, _ in enumerate(stoich):
            elem = stoich[ind][0]
            fraction = stoich[ind][1]
            if '[' in elem or ']' in elem:
                types_dict = dict()
                types_dict['$or'] = list()
                elem = elem.strip('[').strip(']')
                if elem in self._periodic_table:
                    for group_elem in self._periodic_table[elem]:
                        types_dict['$or'].append(dict())
                        types_dict['$or'][-1]['stoichiometry'] = dict()
                        types_dict['$or'][-1]['stoichiometry']['$in'] = [[group_elem, fraction]]
                    query_dict['$and'].append(types_dict)
                elif ',' in elem:
                    for group_elem in elem.split(','):
                        types_dict['$or'].append(dict())
                        types_dict['$or'][-1]['stoichiometry'] = dict()
                        types_dict['$or'][-1]['stoichiometry']['$in'] = [[group_elem, fraction]]
                    query_dict['$and'].append(types_dict)
            else:
                stoich_dict = dict()
                stoich_dict['stoichiometry'] = dict()
                stoich_dict['stoichiometry']['$in'] = [[elem, fraction]]
                query_dict['$and'].append(stoich_dict)
        if not partial_formula:
            size_dict = dict()
            size_dict['stoichiometry'] = dict()
            size_dict['stoichiometry']['$size'] = len(stoich)
            query_dict['$and'].append(size_dict)

        return query_dict

    def _query_ratio(self, ratios):
        """ Query DB for ratio of two elements.

        Input, e.g.:

            ratios = [['MoS', 2],
                      ['LiS', 1]]

        """
        query_dict = dict()
        for pair in ratios:
            query_dict['ratios.' + pair[0]] = pair[1]
        return query_dict

    def _query_composition(self, custom_elem=None, partial_formula=None, elem_field='elems'):
        """ Query DB for all structures containing all the elements
        taken as input. Passing this function a number is a deprecated
        feature, replaced by query_num_species.

        Args:

            | custom_elem     : str, use to query custom string, rather than CLI args
            | partial_formula : bool, remove stoich size from query if True
            | elem_field      : str, which field to query for elems, either `atom_types` or `elems`

        """
        if custom_elem is None:
            if isinstance(self.args.get('composition'), str):
                elements = [self.args.get('composition')]
            else:
                elements = self.args.get('composition')
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
                    sys.exit()

        elements = parse_element_string(elements[0])

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
                        types_dict['$and'][-1][elem_field] = dict()
                        types_dict['$and'][-1][elem_field]['$in'] = [elem]
                    query_dict['$or'].append(types_dict)
        elif non_binary:
            query_dict = dict()
            query_dict['$and'] = []
            size = 0
            for ind, elem in enumerate(elements):
                if elem != ':' and not elem.isdigit():
                    query_dict['$and'].append(self._query_composition(custom_elem=[elem], partial_formula=True))
                    size += 1
                if elem == ':':
                    # convert e.g. MoS2 to [['MoS', 2]]
                    # or LiMoS2 to [['LiMo', 1], ['MoS', '2], ['LiS', 2]]
                    ratio_elements = elements[ind+1:]
                    for _ind, _ in enumerate(ratio_elements):
                        if _ind < len(ratio_elements)-1:
                            if not ratio_elements[_ind].isdigit() and \
                                    not ratio_elements[_ind+1].isdigit():
                                ratio_elements.insert(_ind+1, '1')
                    if not ratio_elements[-1].isdigit():
                        ratio_elements.append('1')
                    ratios = []
                    for _ind in range(0, len(ratio_elements), 2):
                        for _jind in range(_ind, len(ratio_elements), 2):
                            if ratio_elements[_ind] != ratio_elements[_jind]:
                                ratios.append([ratio_elements[_ind] + ratio_elements[_jind],
                                               round(float(ratio_elements[_ind+1]) /
                                                     float(ratio_elements[_jind+1]),
                                                     3)])
                    query_dict['$and'].append(self._query_ratio(ratios))
        else:
            # expand group macros
            query_dict = dict()
            query_dict['$and'] = []
            size = len(elements)
            for ind, elem in enumerate(elements):
                # prototype for chemically motivated searches, e.g. transition metals
                if '[' in elem or ']' in elem:
                    types_dict = dict()
                    types_dict['$or'] = list()
                    elem = elem.strip('[').strip(']')
                    if elem in self._periodic_table:
                        for group_elem in self._periodic_table[elem]:
                            types_dict['$or'].append(dict())
                            types_dict['$or'][-1][elem_field] = dict()
                            types_dict['$or'][-1][elem_field]['$in'] = [group_elem]
                    elif ',' in elem:
                        for group_elem in elem.split(','):
                            types_dict['$or'].append(dict())
                            types_dict['$or'][-1][elem_field] = dict()
                            types_dict['$or'][-1][elem_field]['$in'] = [group_elem]
                else:
                    types_dict = dict()
                    types_dict[elem_field] = dict()
                    types_dict[elem_field]['$in'] = [elem]
                query_dict['$and'].append(types_dict)
        if not partial_formula and not self.args.get('intersection'):
            size_dict = dict()
            size_dict['stoichiometry'] = dict()
            size_dict['stoichiometry']['$size'] = size
            query_dict['$and'].append(size_dict)

        return query_dict

    def _query_num_species(self):
        """ Query database for all structures with a
        given number of elements, e.g. binaries, ternaries etc.
        """
        num = self.args.get('num_species')
        if not isinstance(num, list):
            num = num
        elif isinstance(num, list):
            num = num[0]
        else:
            sys.exit('--num_species takes a single integer or list containing a single integer')
        query_dict = dict()
        query_dict['stoichiometry'] = dict()
        query_dict['stoichiometry']['$size'] = num

        return query_dict

    def _query_space_group(self):
        """ Query DB for all structures with given
        space group.
        """
        query_dict = dict()
        if not isinstance(self.args.get('space_group'), list):
            spg = [self.args.get('space_group')]
        else:
            spg = self.args.get('space_group')
        query_dict['space_group'] = str(spg[0])

        return query_dict

    def _query_num_fu(self):
        """ Query DB for all structures with more than a
        given number of formula units in the simulation.
        """
        query_dict = dict()
        num = self.args.get('num_fu')
        if isinstance(num, list):
            num = num[0]
        query_dict['num_fu'] = dict()
        query_dict['num_fu']['$gte'] = num

        return query_dict

    def _query_tags(self):
        """ Find all structures matching given tags. """
        query_dict = dict()
        query_dict['$and'] = []
        for tag in self.args.get('tags'):
            temp_dict = dict()
            temp_dict['tags'] = dict()
            temp_dict['tags']['$in'] = [tag]
            query_dict['$and'].append(temp_dict)

        return query_dict

    def _query_doi(self):
        """ Find all structures matching given DOI,
        in format xxxx/xxxx.
        """
        doi = self.args.get('doi')
        if not isinstance(doi, list):
            doi = [doi]
        query_dict = dict()
        query_dict['doi'] = dict()
        query_dict['doi']['$in'] = doi

        return query_dict

    def _query_id(self):
        """ Find all structures matching given tags. """
        query_dict = dict()
        query_dict['text_id'] = self.args.get('id')
        return query_dict

    def _query_icsd(self):
        """ Find all structures matching given ICSD CollCode. """
        if not isinstance(self.args.get('icsd'), list):
            icsd = [self.args.get('icsd')]
        else:
            icsd = self.args.get('icsd')
        query_dict = dict()
        query_dict['icsd'] = dict()
        if self.args.get('icsd') == 0:
            query_dict['icsd']['$exists'] = True
        else:
            query_dict['icsd']['$eq'] = str(icsd[0])
        return query_dict

    def _query_source(self):
        """ Find all structures with source string from args. """
        import re
        src_str = self.args.get('src_str')
        if not isinstance(src_str, list):
            src_str = [src_str]
        query_dict = dict()
        query_dict['source'] = dict()
        query_dict['source']['$in'] = [re.compile(src) for src in src_str]
        return query_dict

    def _query_root_source(self):
        """ Find all structures with root source string from args. """
        root_src = self.args.get('root_src')
        if not isinstance(root_src, list):
            root_src = [root_src]
        query_dict = dict()
        for src in root_src:
            query_dict['$or'] = []
            query_dict['$or'].append(dict())
            query_dict['$or'][-1]['root_source'] = src
        return query_dict

    def _query_quality(self):
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

    def _query_pressure(self):
        """ Query pressure, either by an exact match on external_pressure
        or an approximate match on the pressure on the cell, with tolerance
        of either 0.05 GPa or 10%.
        """
        input_pressure = self.args.get('pressure')

        if input_pressure < 0:
            approx_pressure = [1.1*input_pressure-0.05, 0.9*input_pressure+0.05]
        else:
            approx_pressure = [0.9*input_pressure-0.05, 1.1*input_pressure+0.05]

        query_dict = dict()
        query_dict['pressure'] = dict()
        query_dict['pressure']['$lt'] = approx_pressure[1]
        query_dict['pressure']['$gt'] = approx_pressure[0]

        return query_dict

    def _query_encap(self):
        """ Query only CNT encapsulated structures. """
        query_dict = dict()
        query_dict['encapsulated'] = dict()
        query_dict['encapsulated']['$exists'] = True

        return query_dict

    def _query_cnt_radius(self):
        """ Query structures within a nanotube of given radius
        to within a tolerance of 0.01 A.
        """
        query_dict = dict()
        if not isinstance(self.args.get('cnt_radius'), list):
            cnt_rad = [self.args.get('cnt_radius')]
        query_dict['cnt_radius'] = dict()
        query_dict['cnt_radius']['$gt'] = cnt_rad[0] - 0.01
        query_dict['cnt_radius']['$lt'] = cnt_rad[0] + 0.01

        return query_dict

    def _query_cnt_vector(self):
        """ Query structures within a nanotube of given chiral vector. """
        query_dict = dict()
        if not isinstance(self.args.get('cnt_vector'), list) or len(self.args.get('cnt_vector')) != 2:
            sys.exit('CNT vector query needs to be of form [n, m]')
        else:
            chiral_vec = self.args.get('cnt_vector')
        query_dict['cnt_chiral'] = dict()
        query_dict['cnt_chiral']['$eq'] = chiral_vec

        return query_dict

    def _query_cutoff(self):
        """ Query all calculations above given plane-wave cutoff. """
        query_dict = dict()
        query_dict['cut_off_energy'] = dict()
        if not isinstance(self.args.get('cutoff'), list):
            cutoffs = [self.args.get('cutoff')]
        else:
            cutoffs = self.args.get('cutoff')
        if len(cutoffs) == 2:
            if cutoffs[0] > cutoffs[1]:
                sys.exit('Cutoff query needs to be of form [min, max]')
            query_dict['cut_off_energy']['$gte'] = cutoffs[0]
            query_dict['cut_off_energy']['$lte'] = cutoffs[1]
        else:
            query_dict['cut_off_energy']['$gte'] = cutoffs[0]
        return query_dict

    def _query_geom_force_tol(self):
        """ Query all calculations with the correct relaxed force tolerance. """
        query_dict = dict()
        query_dict['geom_force_tol'] = dict()
        if not isinstance(self.args.get('geom_force_tol'), list):
            tols = [self.args.get('geom_force_tol')]
        else:
            tols = self.args.get('geom_force_tol')
        if len(tols) == 2:
            if tols[0] > tols[1]:
                sys.exit('Force tol needs to be of form [min, max]')
            query_dict['geom_force_tol']['$gte'] = tols[0]
            query_dict['geom_force_tol']['$lte'] = tols[1]
        else:
            query_dict['geom_force_tol']['$lte'] = tols[0]
        return query_dict

    def _query_sedc(self):
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

    def _query_xc_functional(self, xc_functional=None):
        """ Query all calculations with specified xc-functional.

        Args:

            | xc_functional: str, CASTEP string for xc-functional to override CLI.

        """
        query_dict = dict()
        if xc_functional is None:
            if isinstance(self.args.get('xc_functional'), list):
                xc_functional = self.args.get('xc_functional')[0]
            else:
                xc_functional = self.args.get('xc_functional')
        if xc_functional is not None:
            query_dict['xc_functional'] = xc_functional.upper()
        return query_dict

    def _query_kpoints(self):
        """ Query all calculations with finer than the given
        kpoint sampling.
        """
        query_dict = dict()
        if not isinstance(self.args.get('mp_spacing'), list):
            mp_spacing = [self.args.get('mp_spacing')]
        else:
            mp_spacing = self.args.get('mp_spacing')
        tol = 0.01
        if self.args.get('kpoint_tolerance') is not None:
            try:
                tol = float(self.args.get('kpoint_tolerance'))
            except Exception:
                print_warning('Failed to read custom kpoint tolerance.')
        query_dict['kpoints_mp_spacing'] = dict()
        query_dict['kpoints_mp_spacing']['$lte'] = mp_spacing[0] + tol
        query_dict['kpoints_mp_spacing']['$gte'] = mp_spacing[0] - tol
        return query_dict

    def _query_spin(self):
        """ Query all calculations with spin polarisation,
        i.e. --spin n!=0, or non-spin-polarization, i.e. --spin 0.
        """
        query_dict = dict()
        if isinstance(self.args.get('spin'), list):
            spin = self.args.get('spin')[0]
        else:
            spin = self.args.get('spin')
        if spin == 'any':
            query_dict = dict()
        elif int(spin) == 0:
            query_dict['spin_polarized'] = dict()
            query_dict['spin_polarized']['$ne'] = True
        elif int(spin) > 0:
            query_dict['spin_polarized'] = True
        return query_dict

    def _query_calc(self, doc):
        """ Find all structures with matching
        accuracy to specified structure.
        """
        self._gs_enthalpy = doc['enthalpy_per_atom']

        query_dict = {}
        query_dict['$and'] = []
        query_dict['$and'].append(self._query_xc_functional(xc_functional=doc.get('xc_functional')))
        if self.args.get('time') is not None:
            query_dict['$and'].append(self._query_time())
        if 'spin_polarized' in doc and doc['spin_polarized']:
            if self.args.get('spin') != 'any':
                temp_dict = dict()
                temp_dict['spin_polarized'] = doc['spin_polarized']
                query_dict['$and'].append(temp_dict)
        else:
            if self.args.get('spin') != 'any':
                temp_dict = dict()
                temp_dict['spin_polarized'] = dict()
                temp_dict['spin_polarized']['$ne'] = True
                query_dict['$and'].append(temp_dict)
        if 'geom_force_tol' in doc and doc['geom_force_tol'] != 0.05:
            temp_dict = dict()
            temp_dict['geom_force_tol'] = doc['geom_force_tol']
            query_dict['$and'].append(temp_dict)
        else:
            temp_dict = dict()
            temp_dict['$or'] = dict()
            temp_dict['$or'] = []
            temp_dict['$or'].append({'geom_force_tol': {'$exists': False}})
            temp_dict['$or'].append({'geom_force_tol': {'$eq': 0.05}})
            query_dict['$and'].append(temp_dict)
        if 'sedc_scheme' in doc:
            temp_dict = dict()
            temp_dict['sedc_scheme'] = doc['sedc_scheme']
            query_dict['$and'].append(temp_dict)
        else:
            temp_dict = dict()
            temp_dict['sedc_scheme'] = dict()
            temp_dict['sedc_scheme']['$exists'] = False
            query_dict['$and'].append(temp_dict)

        db = self.args.get('db')
        if isinstance(db, list):
            db = db[0]

        if self.args.get('loose') or (db is not None and 'oqmd' in db):
            return query_dict
            # temp_dict = dict()
            # query_dict.append(dict())
            # query_dict[-1]['cut_off_energy'] = doc['cut_off_energy']
        else:
            temp_dict = dict()
            temp_dict['kpoints_mp_spacing'] = dict()
            if self.args.get('kpoint_tolerance') is not None:
                tol = float(self.args.get('kpoint_tolerance'))
            else:
                tol = 0.01
            temp_dict['kpoints_mp_spacing']['$gte'] = doc['kpoints_mp_spacing'] - tol
            temp_dict['kpoints_mp_spacing']['$lte'] = doc['kpoints_mp_spacing'] + tol
            query_dict['$and'].append(temp_dict)
            query_dict['$and'].append(dict())
            query_dict['$and'][-1]['cut_off_energy'] = doc['cut_off_energy']
        if 'species_pot' in doc:
            for species in doc['species_pot']:
                temp_dict = dict()
                temp_dict['$or'] = []
                temp_dict['$or'].append(dict())
                temp_dict['$or'][-1]['species_pot.'+species] = dict()
                temp_dict['$or'][-1]['species_pot.'+species]['$exists'] = False
                temp_dict['$or'].append(dict())
                temp_dict['$or'][-1]['species_pot.'+species] = doc['species_pot'][species]
                query_dict['$and'].append(temp_dict)

        if self.debug:
            print('Calc match dict:')
            print(dumps(query_dict, indent=2))

        return query_dict

    def _query_time(self, preceding=True):
        """ Only include structures added before
        the date given in args['time'].
        """
        from datetime import datetime, timedelta
        from time import mktime
        query_dict = dict()
        time_period = timedelta(days=int(self.args.get('time')))
        time = (datetime.today() - time_period).timetuple()
        elapsed = str(hex(int(mktime(time))))[2:]
        cutoff_id = ObjectId(elapsed + '0000000000000000')
        query_dict['_id'] = dict()
        if preceding:
            query_dict['_id']['$lte'] = cutoff_id

        return query_dict


class EmptyCursor:
    """ Empty cursor class for failures. """
    def count(self):
        return 0
