# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements all queries to the database, including parsing
user inputs, displaying results and calling other functionality.

"""


import sys
import random
from os import devnull
from itertools import combinations
from traceback import print_exc

import pymongo as pm
import numpy as np
from bson.json_util import dumps
from bson.objectid import ObjectId

from matador.utils.print_utils import print_warning, print_success, print_notify
from matador.utils.chem_utils import get_periodic_table
from matador.utils.chem_utils import parse_element_string, get_stoich_from_formula
from matador.utils.cursor_utils import display_results, filter_cursor_by_chempots
from matador.db import make_connection_to_collection
from matador.config import load_custom_settings


class DBQuery:
    """ Class that implements queries to MongoDB
    structure database.

    Attributes:
        cursor (list of dict or :obj:`pymongo.Cursor`): list or cursor of structures
            matching query.
        args (dict): contains all keyword arguments used to construct
            the query (see matador query --help) for list.
        query_dict (dict): dictionary passed to database for query
        calc_dict (dict): if performing a matching query (e.g. self.args['subcmd'] = 'hull'), this
            dictionary contains the parameters used to match to other structures
        repo (pymongo.collection.Collection): the pymongo collection that is being queried.
        top (int): number of structures to print/export set by self.args.get('top') (DEFAULT: 10).
        cursor_min_limit (int): if a query returns more structures than this, do not implicitly convert to a list.

    """

    # below this number of documents,
    # all queries will return a list rather than a pymongo Cursor
    cursor_min_limit = 1000

    def __init__(
        self,
        client=False,
        collections=False,
        subcmd='query',
        debug=False,
        quiet=False,
        mongo_settings=None,
        **kwargs
    ):
        """ Parse arguments from matador or API call before calling
        query.

        Keyword arguments:
            client (pm.MongoClient): the MongoClient to connect to.
            collections (dict of pm.collections.Collection): dictionary of pymongo Collections.
            subcmd (str): either 'query' or 'hull', 'voltage', 'hulldiff'.
                These will decide whether calcuation accuracies are matched
                in the final results.

        """
        # read args and set housekeeping
        self.args = kwargs
        self.debug = debug
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd
        if self.args.get('testing') is None:
            self.args['testing'] = False
        if self.args.get('as_crystal') is None:
            self.args['as_crystal'] = False

        if subcmd in ['hull', 'hulldiff', 'voltage'] and self.args.get('composition') is None:
            raise RuntimeError('{} requires composition query'.format(subcmd))

        self._create_hull = (self.args.get('subcmd') in ['hull', 'hulldiff', 'voltage'] or
                             self.args.get('hull_cutoff') is not None)

        # public attributes
        self.cursor = EmptyCursor()
        self.query_dict = None
        self.calc_dict = None
        self.repo = None

        # private attributes to be set later
        self._empty_query = None
        self._gs_enthalpy = None
        self._non_elemental = None
        self._chempots = None
        self._num_to_display = None

        if debug:
            print(self.args)

        if quiet:
            f = open(devnull, 'w')
            sys.stdout = f

        # if testing keyword is used, all database operations are ignored
        if not self.args.get('testing'):

            # connect to db or use passed client
            if client:
                self._client = client
                self._db = client.crystals
            if collections is not False:
                _collections = collections

            if (not collections or not client):
                # use passed settings or load from config file
                if mongo_settings:
                    self.mongo_settings = mongo_settings
                else:
                    self.mongo_settings = load_custom_settings(
                        config_fname=self.args.get('config'), debug=self.args.get('debug')
                    )

                result = make_connection_to_collection(
                    self.args.get('db'), mongo_settings=self.mongo_settings
                )
                # ideally this would be rewritten to use a context manager to ensure
                # that connections are _always_ cleaned up
                self._client, self._db, _collections = result

            if len(_collections) > 1:
                raise NotImplementedError("Querying multiple collections is no longer supported.")
            else:
                for collection in _collections:
                    self._collection = _collections[collection]
                    break

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

            if self.args.get('id') is not None and (self._create_hull or self.args.get('calc_match')):
                # if we've requested and ID and hull/calc_match, do the ID query
                self.perform_id_query()

            self.perform_query()

            if self._create_hull and self.args.get('id') is None:
                # if we're making a normal hull, find the sets of calculations to use
                self.perform_hull_query()

            if not self._create_hull:
                # only filter for uniqueness if not eventually making a hull
                if self.args.get('uniq'):
                    from matador.utils.cursor_utils import filter_unique_structures
                    print_notify('Filtering for unique structures...')

                    if isinstance(self.cursor, pm.cursor.Cursor):
                        raise RuntimeError("Unable to filter pymongo cursor for uniqueness directly.")

                    if self.args.get('top') is not None:
                        top = self.args['top']
                    else:
                        top = len(self.cursor)

                    self.cursor = filter_unique_structures(
                        self.cursor[:top],
                        debug=self.args.get('debug'),
                        sim_tol=self.args.get('uniq'),
                        energy_tol=1e20
                    )

            if self.args.get('available_values') is not None:
                print('Querying available values...')
                self._query_available_values(self.args.get('available_values'), self.cursor)

            # if no client was passed, then we need to close the one we made
            if not client and not self.args.get('testing'):
                self._client.close()

        if quiet:
            f.close()
            sys.stdout = sys.__stdout__

    def _construct_query(self):
        """ Set up query dict and perform query depending on
        command-line / API arguments. Sets self.query_dict.

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
            if not self._create_hull and not self.args.get('calc_match'):
                self.query_dict['$and'].append(self._query_id())
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

        if self.args.get('field') is not None:
            try:
                for ind, field in enumerate(self.args.get('field')):
                    _filter = self.args.get('filter')[ind]
                    try:
                        for i, value in enumerate(_filter):
                            _filter[i] = float(value)
                        filter_type = 'float'
                    except ValueError:
                        filter_type = 'string'

                    if filter_type == 'float':
                        self.query_dict['$and'].append(self._query_float_range(
                            field, _filter))
                    else:
                        self.query_dict['$and'].append(self._query_string(
                            field, _filter))
            except Exception:
                raise RuntimeError(
                    "Unexpected field/filter format. Both must be "
                    "provided as lists, even if only one field is being filtered."
                )
            self._empty_query = False

        if self.args.get('cutoff') is not None:
            self.query_dict['$and'].append(self._query_float_range(
                'cut_off_energy', self.args.get('cutoff')))
            self._empty_query = False

        if self.args.get('geom_force_tol') is not None:
            self.query_dict['$and'].append(self._query_float_range(
                'geom_force_tol', self.args.get('geom_force_tol')))
            self._empty_query = False

        if self.args.get('grid_scale') is not None:
            self.query_dict['$and'].append(self._query_float_range(
                'grid_scale', self.args.get('grid_scale')))
            self._empty_query = False

        if self.args.get('fine_grid_scale') is not None:
            self.query_dict['$and'].append(self._query_float_range(
                'fine_grid_scale', self.args.get('fine_grid_scale')))
            self._empty_query = False

        if self.args.get('src_str') is not None:
            self.query_dict['$and'].append(self._query_source())
            self._empty_query = False

        if self.args.get('root_src') is not None:
            self.query_dict['$and'].append(self._query_root_source())
            self._empty_query = False

        if self.args.get('pressure') is not None:
            self.query_dict['$and'].append(self._query_float_range(
                'pressure', self.args.get('pressure') or 0.0, tolerance=self.args.get('pressure_tolerance') or 0.5))
            self._empty_query = False

        elif self.args['subcmd'] in ['hull', 'hulldiff', 'voltage']:
            self.query_dict['$and'].append(self._query_float_range(
                'pressure', 0.0, tolerance=self.args.get('pressure_tolerance') or 0.5))

        if self.args.get('encapsulated') is True:
            self.query_dict['$and'].append(self._query_encap())
            self._empty_query = False

        if self.args.get('cnt_radius') is not None:
            self.query_dict['$and'].append(self._query_float_range(
                'cnt_radius', self.args.get('cnt_radius'), tolerance=0.01))
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
            self.query_dict['$and'].append(self._query_float_range(
                'kpoints_mp_spacing', self.args.get('mp_spacing'), tolerance=self.args.get('kpoint_tolerance') or 0.01))
            self._empty_query = False

        if self.args.get('spin') is not None:
            tmp_dict = self._query_spin()
            if tmp_dict:
                self.query_dict['$and'].append(tmp_dict)
            self._empty_query = False

        if not self.args.get('ignore_warnings'):
            self.query_dict['$and'].append(self._query_quality())

        if self.args.get('time') is not None:
            self.query_dict['$and'].append(self._query_time(self.args.get('since') or False))
            self._empty_query = False

    def _perform_empty_query(self, as_list=False):
        """ No parameters were asked for, so just return a cursor
        that contains the entire collection currently stored as
        :attr:`repo`.

        Returns:
            list or :obj:`pymongo.cursor.Cursor`: the results of the query.


        """
        num_documents = self.repo.count_documents({})
        cursor = self.repo.find().sort('enthalpy_per_atom', pm.ASCENDING)
        if self.debug:
            print('Empty query, showing all...')

        if num_documents < self.cursor_min_limit or as_list:
            return list(cursor), num_documents

        return cursor, num_documents

    def perform_query(self):
        """ Find results that match the query_dict
        inside the MongoDB database.
        """
        # if no query submitted, find all
        if self._empty_query and self.args.get('id') is None:
            self.repo = self._collection
            self.cursor, cursor_count = self._perform_empty_query()

        # if no special query has been made already, begin executing the query
        if not self._empty_query:
            self.repo = self._collection
            if self.debug:
                print('Query dict:')
                print(dumps(self.query_dict, indent=1))

            # execute query
            self.cursor, cursor_count = self._find_and_sort(self.query_dict)
            if self._non_elemental:
                self.cursor = filter_cursor_by_chempots(self._chempots, self.cursor)

        print('{} results found for query in {}.'.format(cursor_count, self.repo.name))
        self._num_to_display = cursor_count
        if self.args.get('subcmd') != 'swaps' and not self._create_hull:
            self._set_filter_display_results(cursor_count)

        # if a summary has been requested, cursor must be converted to list
        if self.args.get('summary'):
            self.cursor = list(self.cursor)

        if self.args.get('subcmd') != 'swaps' and not self._create_hull:
            if self._num_to_display >= 1 or self._num_to_display is None:
                if self._num_to_display == cursor_count:
                    display_results(self.cursor, **self.args)
                else:
                    display_results(self.cursor[:self._num_to_display], **self.args)

        if isinstance(self.cursor, pm.cursor.Cursor):
            self.cursor.rewind()

    def _set_filter_display_results(self, cursor_count):
        """ Filter and display the results based on the command line parameters. """
        # by default, show the top structures only
        # if delta_E requested, count how many exist below that energy
        if self.args.get('delta_E') is not None:
            if isinstance(self.cursor, pm.cursor.Cursor) and len(self.cursor.distinct('stoichiometry')) > 1:
                print('Multiple stoichiometries in cursor, unable to filter by energy with --delta_E.')
            else:
                self.cursor = list(self.cursor)
                gs_enthalpy = self.cursor[0]['enthalpy_per_atom']
                for ind, doc in enumerate(self.cursor[1:]):
                    if abs(doc['enthalpy_per_atom'] - gs_enthalpy) > self.args.get('delta_E'):
                        self._num_to_display = ind + 1
                        break

        elif self.top == -1 or self.top is None or cursor_count <= self.top:
            self._num_to_display = cursor_count
            self.top = cursor_count
        elif cursor_count > self.top:
            self._num_to_display = self.top

    def _find_and_sort(self, query_filter=None, as_list=False, **kwargs):
        """ Query `self.repo` using Pymongo arguments/kwargs. Sorts based
        on enthalpy_per_atom and optionally returns list of Crystals.

        Keyword arguments:
            query_filter (dict): the query to use. If None, perform a blank query.
            as_list (bool): whether to return a list of a pm.cursor.Cursor object.

        Returns:
            list/pm.cursor.Cursor: the results of the query.
            int: the number of results in the query.

        """
        from matador.crystal import Crystal
        if query_filter is None:
            query_filter = {}
        count = self.repo.count_documents(query_filter, **kwargs)
        cursor = self.repo.find(query_filter, **kwargs).sort('enthalpy_per_atom', pm.ASCENDING)

        if self.args.get('as_crystal'):
            return [Crystal(doc) for doc in cursor], count
        if count < self.cursor_min_limit or as_list:
            return list(cursor), count

        return cursor, count

    def perform_hull_query(self):
        """ Perform the multiple queries necessary to find possible
        calculation sets to create a convex hull from.

        Raises:
            SystemExit: if no structures are found for hull.

        """
        if self._collection is not None:
            self.repo = self._collection
            print('Creating hull from structures in query results.')
            if self.args.get('biggest'):
                print('\nFinding biggest calculation set for hull...\n')
            else:
                print('\nFinding the best calculation set for hull...')

            test_cursors = []
            test_cursor_count = []
            text_ids = []
            calc_dicts = []
            cutoff = []

            num_sample = 2
            num_rand_sample = 5 if self.args.get('biggest') else 3

            if isinstance(self.cursor, pm.cursor.Cursor):
                count = self.cursor.count()
            else:
                count = len(self.cursor)

            if count <= 0:
                raise SystemExit('No structures found for hull.')

            # generate some random indices to match to, make sure they are in order
            # so can be accessed without cursor rewinds
            sampling_indices = list(range(num_sample)) + sorted(random.sample(range(2, count), num_rand_sample))

            for ind in sampling_indices:
                doc = self.cursor[ind]
                text_ids.append(doc['text_id'])
                try:
                    self.query_dict = self._query_calc(doc)
                    cutoff.append(doc['cut_off_energy'])
                    calc_dicts.append(dict())
                    calc_dicts[-1]['$and'] = list(self.query_dict['$and'])
                    self.query_dict['$and'].append(self._query_composition())
                    if not self.args.get('ignore_warnings'):
                        self.query_dict['$and'].append(self._query_quality())

                    probe_cursor, probe_count = self._find_and_sort(self.query_dict)

                    if self._non_elemental:
                        probe_cursor = filter_cursor_by_chempots(self._chempots, probe_cursor)
                        probe_count = len(probe_cursor)

                    test_cursors.append(probe_cursor)
                    test_cursor_count.append(probe_count)

                    print("{:^24}: matched {} structures."
                          .format(' '.join(doc['text_id']), probe_count),
                          end='\t-> ')
                    print('{spin}{sedc}{functional} {cutoff} eV, {geom_force_tol} eV/A, {kpoints} 1/A.'
                          .format(spin="S-" if doc.get('spin_polarized') else '',
                                  sedc="+" + doc.get('sedc') + "+" if doc.get('sedc') else "",
                                  functional=doc["xc_functional"],
                                  cutoff=doc["cut_off_energy"],
                                  geom_force_tol=doc.get('geom_force_tol', 'xxx'),
                                  kpoints=doc.get('kpoints_mp_spacing', 'xxx')))

                    if test_cursor_count[-1] == count:
                        print('Matched all structures...')
                        break
                    if test_cursor_count[-1] > 2 * int(count / 3):
                        print('Matched at least 2/3 of total number, composing hull...')
                        break

                except Exception:
                    print_exc()
                    print_warning('Error with {}'.format(' '.join(doc['text_id'])))

            if self.args.get('biggest'):
                choice = np.argmax(np.asarray(test_cursor_count))
            else:
                # by default, find highest cutoff hull as first proxy for quality
                choice = np.argmax(np.asarray(cutoff))

            text_id = text_ids[choice]
            self.cursor = test_cursors[choice]
            self.calc_dict = calc_dicts[choice]
            if not test_cursor_count[choice]:
                raise RuntimeError('No structures found that match chemical potentials.')

            print_success('Composing hull from set containing {}'.format(' '.join(text_id)))

    def perform_id_query(self):
        """ Query the `text_id` field for the ID provided in the args for a calc_match
        or hull/voltage query. Use the results of the text_id query to match to other
        entries that have the same calculation parameters. Sets self.query_dict and
        self.calc_dict.

        Raises:
            RuntimeError: if no structures are found.

        """

        self.cursor = []
        query_dict = dict()
        query_dict['$and'] = []
        query_dict['$and'].append(self._query_id())
        if not self.args.get('ignore_warnings'):
            query_dict['$and'].append(self._query_quality())
        self.repo = self._collection
        self.cursor = list(self._find_and_sort(query_dict))

        if not self.cursor:
            raise RuntimeError('Could not find a match with {} try widening your search.'.format(self.args.get('id')))

        if len(self.cursor) >= 1:
            display_results(list(self.cursor)[:self.top], **self.args)

            if len(self.cursor) > 1:
                print_warning('Matched multiple structures with same text_id. The first one will be used.')

            # save special copy of calc_dict for hulls
            self.calc_dict = dict()
            self.calc_dict['$and'] = []
            # to avoid deep recursion, and since this is always called first
            # don't append, just set
            self.query_dict = self._query_calc(self.cursor[0])
            if self.args.get('composition'):
                self.args['intersection'] = True
                self.query_dict['$and'].append(self._query_composition())
            self.calc_dict['$and'] = list(self.query_dict['$and'])

    def query_stoichiometry(self, **kwargs):
        """ Alias for private function of the same name. """
        return self._query_stoichiometry(**kwargs)

    def query_composition(self, **kwargs):
        """ Alias for private function of the same name. """
        return self._query_composition(**kwargs)

    def query_tags(self, **kwargs):
        """ Alias for private function of the same name. """
        return self._query_tags(**kwargs)

    def query_quality(self, **kwargs):
        """ Alias for private function of the same name. """
        return self._query_quality(**kwargs)

    @staticmethod
    def _query_float_range(field, values, tolerance=None):
        """ Query all entries with field between float value range,
        or with float value.

        Parameters:
            field (str): the field to query.
            values (float/list of float): either single value, or list
                of 2 floats.

        Keyword arguments:
            tolerance (float): tolerance to add and subtract if single value is provided.

        Returns:
            dict: the constructed query.

        """
        query_dict = dict()
        query_dict[field] = dict()
        if not isinstance(values, list):
            values = [values]
        if len(values) == 2:
            if values[0] > values[1]:
                tmp = values[0]
                values[0] = values[1]
                values[1] = tmp

            query_dict[field]['$gte'] = values[0]
            query_dict[field]['$lte'] = values[1]
        else:
            if tolerance is None:
                query_dict[field]['$eq'] = values[0]
            else:
                query_dict[field]['$gte'] = round(values[0] - tolerance, 8)
                query_dict[field]['$lte'] = round(values[0] + tolerance, 8)
        return query_dict

    @staticmethod
    def _query_string(field, values):
        """ Query all entries for an exact string match on field.

        Parameters:
            field (str): the field to query.
            values (list or str): strings to query with $or joins.

        Returns:
            dict: the constructed query.

        """
        query_dict = dict()
        if not isinstance(values, list):
            values = [values]
        if len(values) > 1:
            query_dict['$or'] = []
            for value in values:
                query_dict['$or'].append({field: value})
        else:
            query_dict[field] = values[0]

        return query_dict

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
            raise RuntimeError('Formula cannot contain ":", you probably meant to query composition.')

        stoich = get_stoich_from_formula(stoich[0], sort=False)

        query_dict = dict()
        query_dict['$and'] = []

        for ind, _ in enumerate(stoich):
            elem = stoich[ind][0]
            fraction = int(stoich[ind][1])

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

    @staticmethod
    def _query_ratio(ratios):
        """ Query DB for ratio of two elements. Ratios must be integers.

        Parameters:
            ratios (list): e.g.  ratios = [['MoS', 2], ['LiS', 1]]

        """
        query_dict = dict()
        for pair in ratios:
            query_dict['ratios.' + pair[0]] = int(pair[1])
        return query_dict

    def _query_composition(self, custom_elem=None, partial_formula=None, elem_field='elems'):
        """ Query DB for all structures containing all the elements
        taken as input. Passing this function a number is a deprecated
        feature, replaced by query_num_species.

        Keyword arguments:
            custom_elem (str): use to query custom string, rather than CLI args
            partial_formula (bool): remove stoich size from query if True
            elem_field (str): which field to query for elems, either `atom_types` or `elems`

        Returns:
            dict: dictionary containing database query.

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

        self._non_elemental = False
        if ':' in elements[0]:
            self._non_elemental = True
            self.args['intersection'] = True
            self._chempots = elements[0].split(':')
            elements = [parse_element_string(elem) for elem in self._chempots]
            elements = list(dict.fromkeys([char for elem in elements for char in elem if char.isalpha()]))
        # if there's only one string, try split it by caps
        if not self._non_elemental:
            for char in elements[0]:
                if char.isdigit():
                    raise SystemExit('Composition cannot contain a number.')
            elements = parse_element_string(elements[0])

        or_preference = False
        for _, elem in enumerate(elements):
            if '{' in elem or '}' in elem:
                or_preference = True

        elements_tmp = [element for ind, element in enumerate(elements)
                        if element not in elements[:ind]]
        if len(elements_tmp) < len(elements):
            print('Ignoring duplicate element...')
        elements = elements_tmp

        if self.args.get('intersection'):
            if or_preference:
                raise RuntimeError('Intersection not implemented for overlapping sets, e.g. {}')

            query_dict = dict()
            query_dict['$or'] = []
            size = len(elements)
            # iterate over all combinations, limited by num species
            if self.args.get('num_species'):
                max_num = self.args.get('num_species')
                min_num = max_num
            else:
                max_num = 8
                min_num = 1
                if len(elements) > max_num:
                    print('Limiting query to up to {} elements per structure...'.format(max_num))
            for rlen in range(min_num, max_num+1):
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
        else:
            # expand group macros
            query_dict = dict()
            query_dict['$and'] = []
            size = len(elements)

            if or_preference:
                element_slots = []
                for elem in elements:
                    if '[' in elem or '{' in elem:
                        elem = elem.strip('{').strip('}').strip('[').strip(']')
                        if elem in self._periodic_table:
                            element_slots.append(self._periodic_table[elem])
                        elif ',' in elem:
                            element_slots.append(elem.split(','))
                        else:
                            element_slots.append([elem])
                    else:
                        element_slots.append([elem])

                from itertools import product
                slots = [list(config) for config in product(*element_slots)]
                types_dict = dict()
                types_dict['$or'] = list()
                for slot in slots:
                    if len({elem for elem in slot}) == len(slot):
                        types_dict['$or'].append(dict())
                        types_dict['$or'][-1]['$and'] = []
                        for elem in slot:
                            types_dict['$or'][-1]['$and'].append(dict())
                            types_dict['$or'][-1]['$and'][-1][elem_field] = dict()
                            types_dict['$or'][-1]['$and'][-1][elem_field]['$in'] = [elem]
                query_dict['$and'].append(types_dict)

            else:
                for elem in elements:
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
        if isinstance(self.args.get('id'), str):
            self.args['id'] = self.args['id'].strip().split(' ')
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
        if isinstance(icsd[0], bool):
            query_dict['icsd'] = dict()
            query_dict['icsd']['$exists'] = icsd[0]
        elif icsd[0] == 0:
            query_dict['icsd'] = dict()
            query_dict['icsd']['$exists'] = True
        else:
            query_dict['$or'] = [{'icsd': {'$eq': str(icsd[0])}}, {'icsd': {'$eq': icsd[0]}}]
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

    @staticmethod
    def _query_available_values(field, cursor):
        """ Query the values stored under a particular field and
        print the information.

        Parameters:
            field (str): the field to query.
            cursor (list): the cursor to query.

        """
        supported_fields = [
            'doi',
            'tags',
            'root_source',
            'cnt_vector',
            'castep_version',
            'cut_off_energy'
        ]
        number_containing_field = sum([1 for doc in cursor if field in doc])
        if field in supported_fields and number_containing_field != 0:
            value_degeneracy = dict()
            for doc in cursor:
                if doc.get(field) is not None:
                    values = doc.get(field)
                    if isinstance(values, list):
                        for value in values:
                            if value in value_degeneracy:
                                value_degeneracy[value] += 1
                            else:
                                value_degeneracy[value] = 1

                    else:
                        if values in value_degeneracy:
                            value_degeneracy[values] += 1
                        else:
                            value_degeneracy[values] = 1

            print('Set of values under key {}:'.format(field))
            for value in sorted(value_degeneracy, key=value_degeneracy.get):
                print('{:<10} -> {:<10}'.format(value_degeneracy[value], value))
        else:
            print('Field {} unsupported for finding all possible values, must be one of {}'
                  .format(field, supported_fields))
            print('{}/{} contain field {}'.format(number_containing_field, len(cursor), field))

    @staticmethod
    def _query_quality():
        """ Find all structures with non-zero or non-existent (e.g.
        OQMD) quality.

        """
        query_dict = dict()
        query_dict['$or'] = []
        query_dict['$or'].append(dict())
        query_dict['$or'][-1]['quality'] = dict()
        query_dict['$or'][-1]['quality']['$gt'] = 0
        query_dict['$or'].append(dict())
        query_dict['$or'][-1]['quality'] = dict()
        query_dict['$or'][-1]['quality']['$exists'] = False

        return query_dict

    @staticmethod
    def _query_encap():
        """ Query only CNT encapsulated structures. """
        query_dict = dict()
        query_dict['encapsulated'] = dict()
        query_dict['encapsulated']['$exists'] = True

        return query_dict

    def _query_cnt_vector(self):
        """ Query structures within a nanotube of given chiral vector. """
        query_dict = dict()
        if not isinstance(self.args.get('cnt_vector'), list) or len(self.args.get('cnt_vector')) != 2:
            raise SystemExit('CNT vector query needs to be of form [n, m]')

        chiral_vec = self.args.get('cnt_vector')
        query_dict['cnt_chiral'] = dict()
        query_dict['cnt_chiral']['$eq'] = chiral_vec

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

        Keyword arguments:

            xc_functional (str): CASTEP string for xc-functional to
                override CLI.

        """
        query_dict = dict()
        if xc_functional is None:
            if isinstance(self.args.get('xc_functional'), list):
                xc_functional = self.args.get('xc_functional')[0]
            else:
                xc_functional = self.args.get('xc_functional')
        if xc_functional is not None:
            query_dict['xc_functional'] = xc_functional
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
        query_dict['$and'].append(self._query_float_range(
            'pressure', doc.get('pressure', 0.0), tolerance=self.args.get('pressure_tolerance') or 0.5))

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
        if doc.get('grid_scale', 1.75) == 1.75:
            temp_dict = dict()
            temp_dict['$or'] = []
            temp_dict['$or'].append({'grid_scale': {'$exists': False}})
            temp_dict['$or'].append({'grid_scale': {'$eq': 1.75}})
            query_dict['$and'].append(temp_dict)
        else:
            temp_dict = dict()
            temp_dict['$or'] = []
            temp_dict['$or'].append({'grid_scale': {'$eq': doc.get('grid_scale')}})
            query_dict['$and'].append(temp_dict)
        if doc.get('fine_grid_scale', 1.75) == 1.75:
            temp_dict = dict()
            temp_dict['$or'] = []
            temp_dict['$or'].append({'fine_grid_scale': {'$exists': False}})
            temp_dict['$or'].append({'fine_grid_scale': {'$eq': 1.75}})
            query_dict['$and'].append(temp_dict)
        else:
            temp_dict = dict()
            temp_dict['$or'] = []
            temp_dict['$or'].append({'fine_grid_scale': {'$eq': doc.get('fine_grid_scale')}})
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

        query_dict['$and'].append(self._query_float_range(
            'kpoints_mp_spacing', doc.get('kpoints_mp_spacing'), tolerance=self.args.get('kpoint_tolerance') or 0.01))
        query_dict['$and'].append(dict())
        query_dict['$and'][-1]['cut_off_energy'] = doc['cut_off_energy']

        if 'species_pot' in doc:
            for species in doc['species_pot']:
                temp_dict = dict()
                temp_dict['$or'] = []
                temp_dict['$or'].append(dict())
                temp_dict['$or'][-1]['species_pot.' + species] = dict()
                temp_dict['$or'][-1]['species_pot.' + species]['$exists'] = False
                temp_dict['$or'].append(dict())
                temp_dict['$or'][-1]['species_pot.' + species] = doc['species_pot'][species]
                query_dict['$and'].append(temp_dict)

        if self.debug:
            print('Calc match dict:')
            print(dumps(query_dict, indent=2))

        return query_dict

    def _query_time(self, since=False):
        """ Only include structures added before or after (depending on
        since) the date given in args['time'].

        Keyword arguments:
            since (bool): query before or after this time.

        """
        from datetime import datetime, timedelta
        from time import mktime
        query_dict = dict()
        time_period = timedelta(days=int(self.args.get('time')))
        time = (datetime.today() - time_period).timetuple()
        elapsed = str(hex(int(mktime(time))))[2:]
        cutoff_id = ObjectId(elapsed + '0000000000000000')
        query_dict['_id'] = dict()
        if since:
            query_dict['_id']['$gte'] = cutoff_id
        else:
            query_dict['_id']['$lte'] = cutoff_id

        return query_dict


class EmptyCursor:
    """ Empty cursor class for failures. """

    @staticmethod
    def count_documents(*args, **kwargs):
        """ Dummy function always returns 0. """
        return 0
