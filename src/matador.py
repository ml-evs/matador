#!/usr/bin/python
# coding: utf-8
""" This file implements all queries to the database,
including parsing user inputs, displaying results
and calling other functionality. """
from __future__ import print_function
# matador modules
from query import DBQuery
from hull import QueryConvexHull
from scrapers.spatula import Spatula
# import external libraries
import pymongo as pm
import numpy as np
import argparse
import string
from sys import argv
from os import uname
from copy import deepcopy


class Matador:
    """ Class that implements the interface
    to MongoDB structure repository.
    """
    def __init__(self, *args, **kwargs):
        """ Initialise the query with command line
        arguments and return results.
        """
        # read args
        self.kwargs = kwargs
        self.args = vars(args[0])
        # connect to MonoDB structure repository
        local = uname()[1]
        if local == 'cluster2':
            remote = 'node1'
        else:
            remote = None
        self.client = pm.MongoClient(remote)
        self.db = self.client.crystals
        # choose desired collections
        self.collections = dict()
        if self.args['db'] is not None:
            for database in self.args['db']:
                if database == 'all':
                    self.collections = dict()
                    self.collections['ajm'] = self.db['repo']
                    self.collections['oqmd'] = self.db['oqmd']
                elif database == 'ajm':
                    database = 'repo'
                    self.collections['ajm'] = self.db['repo']
                else:
                    self.collections[database] = self.db[database]
        else:
            self.collections['repo'] = self.db['repo']
        # print last spatula report
        self.report = self.client.crystals.spatula
        self.print_report()
        if self.args['subcmd'] == 'stats':
            self.stats()
        if self.args['subcmd'] == 'import':
            self.importer = Spatula(dryrun=self.args['dryrun'],
                                    debug=self.args['debug'],
                                    scratch=self.args['scratch'],
                                    verbosity=self.args['verbosity'],
                                    tags=self.args['tags'])
        if self.args['subcmd'] == 'query':
            self.query = DBQuery(self.client, self.collections, self.args)
        if self.args['subcmd'] == 'hull':
            self.query = DBQuery(self.client, self.collections, self.args)
            self.hull = QueryConvexHull(self.query)

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
                swap = np.random.randint(0, len(swapAtoms) - 1, size=2)
                if swap[0] != swap[1] and swapAtoms[swap[0]] != swapAtoms[swap[1]]:
                        valid = True
            swapAtoms[swap[1]], swapAtoms[swap[0]] = swapAtoms[swap[0]], swapAtoms[swap[1]]
        swapPos = np.asarray(swapDoc['positions_frac'])
        for i in range(len(swapAtoms)):
            swapPos[i] += np.random.rand(3) * (0.1 / 7.9)
        swapDoc['positions_frac'] = swapPos
        hash = self.generate_hash(8)
        self.doc2cell(swapDoc, name + '-' + hash)
        self.doc2param(swapDoc, name + '-' + hash, template_param)

    def generate_hash(self, hashLen=6):
        """ Quick hash generator, based on implementation in PyAIRSS by J. Wynn. """
        hashChars = [str(x) for x in range(0, 10)] + [x for x in string.ascii_lowercase]
        hash = ''
        for i in range(hashLen):
            hash += np.random.choice(hashChars)
        return hash

    def print_report(self):
        """ Print spatula report on current database. """
        try:
            report = self.report.find_one()
            print('Database last modified on', report['last_modified'], 'with spatula',
                  report['version'], 'changeset (' + report['git_hash'] + ').')
        except:
            print('Failed to print database report: spatula is probably running!')

    def stats(self):
        """ Print some useful stats about the database. """
        comp_list = dict()
        overall_stats_dict = dict()
        overall_stats_dict['count'] = 0
        overall_stats_dict['avgObjSize'] = 0
        overall_stats_dict['storageSize'] = 0
        for collection in self.collections:
            db_stats_dict = self.db.command('collstats', collection)
            overall_stats_dict['count'] += db_stats_dict['count']
            overall_stats_dict['avgObjSize'] += db_stats_dict['avgObjSize']
            overall_stats_dict['storageSize'] += db_stats_dict['storageSize']
        print('The collection(s) queried in', self.db.name, 'contain',
              overall_stats_dict['count'], 'structures at',
              "{:.1f}".format(overall_stats_dict['avgObjSize'] / (1024 * len(self.collections))),
              'kB each, totalling', "{:.1f}".format(overall_stats_dict['storageSize'] / (1024**2)),
              'MB when padding is included.')
        for collname in self.collections:
            cursor = self.collections[collname].find()
            for doc in cursor:
                temp = ''
                for ind, elem in enumerate(doc['stoichiometry']):
                    temp += str(elem[0])
                    if ind != len(doc['stoichiometry']) - 1:
                        temp += '+'
                if temp not in comp_list:
                    comp_list[temp] = 0
                comp_list[temp] += 1
        keys = list(comp_list.keys())
        vals = list(comp_list.values())
        comp_list = zip(keys, vals)
        comp_list.sort(key=lambda t: t[1], reverse=True)
        small_count = 0
        first_ind = 1000
        cutoff = 100
        for ind, comp in enumerate(comp_list):
            if comp[1] < cutoff:
                if ind < first_ind:
                    first_ind = ind
                small_count += comp[1]
        comp_list = comp_list[:first_ind]
        comp_list.append(['others < ' + str(cutoff), small_count])
        comp_list.sort(key=lambda t: t[1], reverse=True)
        try:
            from ascii_graph import Pyasciigraph
            from ascii_graph.colors import Gre, Blu, Red
            from ascii_graph.colordata import hcolor
        except:
            exit('Pyascii graph missing; not printing detailed stats.')
        graph = Pyasciigraph(line_length=80, multivalue=False)
        thresholds = {int(overall_stats_dict['count'] / 40): Gre,
                      int(overall_stats_dict['count'] / 10): Blu,
                      int(overall_stats_dict['count'] / 4): Red}
        data = hcolor(comp_list, thresholds)
        for line in graph.graph(label=None, data=data):
            print(line)
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='MATADor',
        description='MATerial and Atomic Database Of Refined structures.',
        epilog='Written by Matthew Evans (2016).')
    # define subparsers for subcommands
    subparsers = parser.add_subparsers(title='subcommands',
                                       description='valid sub-commands',
                                       dest='subcmd')
    # define parent parser for global arguments
    global_parser = argparse.ArgumentParser(add_help=False)
    # common arguments to all subcommands
    global_parser.add_argument('--db', choices=['ajm', 'oqmd', 'all', 'scratch'],
                               help='choose which databases to query. \
                                     NB: "all" does not include scratch',
                               nargs='+')
    # define structure parser for structure query strings
    structure_parser = argparse.ArgumentParser(add_help=False)
    structure_parser.add_argument('-c', '--composition', type=str, nargs='+',
                                  help='find all structures containing exclusively the given elements, \
                                        e.g. LiSi.')
    structure_parser.add_argument('-f', '--formula', type=str, nargs='+',
                                  help='query a particular chemical formula, e.g. GeTeSi3')
    structure_parser.add_argument('-i', '--id', type=str, nargs='+',
                                  help='specify a particular structure by its text_id')
    # define subcommand parsers and their arguments
    stat_parser = subparsers.add_parser('stats',
                                        help='print some stats about the database.',
                                        parents=[global_parser])
    query_parser = subparsers.add_parser('query',
                                         help='query and extract structures from the database',
                                         parents=[global_parser, structure_parser])
    query_parser.add_argument('-s', '--summary', action='store_true',
                              help='show only the ground state for each stoichiometry.')
    query_parser.add_argument('-t', '--top', type=int,
                              help='number of structures to show (DEFAULT: 10)')
    query_parser.add_argument('-d', '--details', action='store_true',
                              help='show as much detail about calculation as possible')
    query_parser.add_argument('-sg', '--space_group',
                              help='query a particular space group')
    query_parser.add_argument('-p', '--pressure', type=float,
                              help='specify an isotropic external pressure to search for \
                                   , e.g. 10 (GPa)')
    query_parser.add_argument('--source', action='store_true',
                              help='print filenames from which structures were wrangled')
    query_parser.add_argument('-ac', '--calc-match', action='store_true',
                              help='display calculations of the same accuracy as specified id')
    query_parser.add_argument('-pf', '--partial-formula', action='store_true',
                              help=('stoichiometry/composition queries will include other \
                                    unspecified species, e.g. -pf search for Li will query \
                                    any structure containing Li, not just pure Li.'))
    query_parser.add_argument('--encap', action='store_true',
                              help='query only structures encapsulated in a carbon nanotube.')
    query_parser.add_argument('--tags', nargs='+', type=str,
                              help=('search for up to 3 manual tags at once'))
    query_parser.add_argument('--cell', action='store_true',
                              help='export query to .cell files in folder name from query string')
    query_parser.add_argument('--res', action='store_true',
                              help='export query to .res files in folder name from query string')
    import_parser = subparsers.add_parser('import',
                                          help='import structures into database; does not \
                                                care about non-unique structures')
    import_parser.add_argument('-d', '--dryrun', action='store_true',
                               help='run the importer without connecting to the database')
    import_parser.add_argument('-v', '--verbosity', action='count',
                               help='enable verbose output')
    import_parser.add_argument('-t', '--tags', nargs='+', type=str,
                               help='set user tags, e.g. nanotube, project name')
    import_parser.add_argument('--debug', action='store_true',
                               help='enable debug output to print every dict')
    import_parser.add_argument('-s', '--scratch', action='store_true',
                               help='import to scratch collection')
    hull_parser = subparsers.add_parser('hull',
                                        help='create a convex hull from query results \
                                        (currently limited to binaries)',
                                        parents=[global_parser, structure_parser])
    voltage_parser = subparsers.add_parser('voltage',
                                           help='plot a voltage curve from query results \
                                           (currently limited to binaries)',
                                           parents=[global_parser, structure_parser])
    swaps_parser = subparsers.add_parser('swaps',
                                         help='perform atomic swaps on query results',
                                         parents=[global_parser, structure_parser])
    # parser.add_argument('--strict', action='store_true',
                        # help=('strictly matches with calc_match,'
                              # 'useful for hulls where convergence is rough'))
    # parser.add_argument('--loose', action='store_true',
                        # help=('loosely matches with calc_match, i.e. only matches' +
                              # 'pspot and xc_functional'))
    # parser.add_argument('--ignore_warnings', action='store_true',
                        # help='includes possibly bad structures')
    # parser.add_argument('--dis', action='store_true',
                        # help='smear hull with local stoichiometry')
    # parser.add_argument('--scratch', action='store_true',
                        # help='query local scratch collection')
    # parser.add_argument('--cell', action='store_true',
                        # help='export query to .cell files in folder name from query string')
    # parser.add_argument('--res', action='store_true',
                        # help='export query to .res files in folder name from query string')
    # parser.add_argument('--debug', action='store_true',
                        # help='print some useful (to me) debug info')
    # parser.add_argument('--write_pressure', nargs='+', type=str,
                        # help=('pressure to add to new cell file, either one float' +
                              # 'for isotropic or 6 floats for anisotropic.'))
    args = parser.parse_args()
    matador = Matador(args, argstr=argv[1:])
