#!/usr/bin/python
# coding: utf-8
""" This file implements all queries to the database,
including parsing user inputs, displaying results
and calling other functionality. """
from __future__ import print_function
# matador modules
from query import DBQuery
from hull import QueryConvexHull
from polish import Polisher
from spatula import Spatula
# import external libraries
import pymongo as pm
import argparse
from sys import argv
from os import uname


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
        if self.args['subcmd'] not in ['import', 'rebuild']:
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
        if self.args['subcmd'] in ['import', 'rebuild']:
            self.importer = Spatula(self.args)
        if self.args['subcmd'] == 'query':
            self.query = DBQuery(self.client, self.collections, self.args)
        if self.args['subcmd'] == 'swaps':
            self.query = DBQuery(self.client, self.collections, self.args)
            if self.args['hull_cutoff'] is not None:
                self.hull = QueryConvexHull(self.query, self.args)
                self.swaps = Polisher(self.hull.hull_cursor, self.args)
            else:
                self.swaps = Polisher(self.query.cursor, self.args)
        if self.args['subcmd'] == 'polish':
            self.query = DBQuery(self.client, self.collections, self.args)
            self.polish = Polisher(self.query.cursor, self.args)
        if self.args['subcmd'] == 'hull':
            self.query = DBQuery(self.client, self.collections, self.args)
            self.hull = QueryConvexHull(self.query, self.args)
            if len(self.hull.hull_cursor) == 0:
                print('No structures on hull with chosen chemical potentials.')
            else:
                print(len(self.hull.hull_cursor), 'structures within', self.hull.hull_cutoff,
                      'eV of the hull with chosen chemical potentials.')
                self.query.display_results(self.hull.hull_cursor)
        if self.args['subcmd'] == 'voltage':
            self.query = DBQuery(self.client, self.collections, self.args)
            self.hull = QueryConvexHull(self.query, self.args)

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
    global_parser.add_argument('--db',
                               help='choose which databases to query',
                               nargs='+')
    # define structure parser for structure query strings
    structure_parser = argparse.ArgumentParser(add_help=False)
    structure_parser.add_argument('-c', '--composition', type=str, nargs='+',
                                  help='find all structures containing exclusively the given elements, \
                                        e.g. LiSi. Macros defined for groups [I]-[VII], [Tran] \
                                        [Lan] and [Act], used with square brackets.')
    structure_parser.add_argument('-n', '--num_species', type=int, nargs='+',
                                  help='find all structures containing a certain \
                                        number of species.')
    structure_parser.add_argument('-f', '--formula', type=str, nargs='+',
                                  help='query a particular chemical formula, e.g. GeTeSi3')
    structure_parser.add_argument('-i', '--id', type=str, nargs='+',
                                  help='specify a particular structure by its text_id')
    structure_parser.add_argument('-ac', '--calc-match', action='store_true',
                                  help='display calculations of the same accuracy as specified id')
    structure_parser.add_argument('-z', '--num_fu', type=int,
                                  help='query a calculations with more than n formula units')
    structure_parser.add_argument('-sg', '--space_group',
                                  help='query a particular space group')
    structure_parser.add_argument('-p', '--pressure', type=float,
                                  help='specify an isotropic external pressure to search for \
                                        , e.g. 10 (GPa)')
    structure_parser.add_argument('-pf', '--partial-formula', action='store_true',
                                  help='stoichiometry/composition queries will include other \
                                        unspecified species, e.g. -pf search for Li will query \
                                        any structure containing Li, not just pure Li.')
    structure_parser.add_argument('--tags', nargs='+', type=str,
                                  help=('search for up to 3 manual tags at once'))
    structure_parser.add_argument('--encap', action='store_true',
                                  help='query only structures encapsulated in a carbon nanotube.')
    structure_parser.add_argument('--loose', action='store_true',
                                  help='loosely matches with calc_match, i.e. only matches \
                                        pspot and xc_functional')
    # define material parser for hull/voltage arguments
    material_parser = argparse.ArgumentParser(add_help=False)
    material_parser.add_argument('--include_oqmd', action='store_true',
                                 help='include OQMD structures on hull and voltage curve.')
    material_parser.add_argument('-hc', '--hull_cutoff', type=float,
                                 help='return only structures within a certain distance from hull')
    material_parser.add_argument('--biggest', action='store_true',
                                 help='try to use the largest subset of structures to create a hull')
    plot_parser = argparse.ArgumentParser(add_help=False)
    plot_parser.add_argument('--png', action='store_true',
                             help='save png rather than showing plot in X')
    plot_parser.add_argument('--subplot', action='store_true',
                             help='plot combined hull and voltage graph')
    spatula_parser = argparse.ArgumentParser(add_help=False)
    spatula_parser.add_argument('-d', '--dryrun', action='store_true',
                                help='run the importer without connecting to the database')
    spatula_parser.add_argument('-v', '--verbosity', action='count',
                                help='enable verbose output')
    spatula_parser.add_argument('-t', '--tags', nargs='+', type=str,
                                help='set user tags, e.g. nanotube, project name')
    spatula_parser.add_argument('--debug', action='store_true',
                                help='enable debug output to print every dict')
    spatula_parser.add_argument('-s', '--scan', action='store_true',
                                help='only scan the database for new structures, do not dictify')
    # define parser for output of new files for swaps/polishes
    collection_parser = argparse.ArgumentParser(add_help=False)
    collection_parser.add_argument('--to', type=str,
                                   help='the text_id of a structure with the desired parameters')
    collection_parser.add_argument('--with', type=str,
                                   help=('the seedname (must be within pwd) of cell and param ' +
                                         'files to use for polishing/swaps'))
    collection_parser.add_argument('--prefix', type=str,
                              help='add a prefix to all file names to write out (auto-appended \
                              with an underscore')
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
    query_parser.add_argument('--source', action='store_true',
                              help='print filenames from which structures were wrangled')
    query_parser.add_argument('--cell', action='store_true',
                              help='export query to .cell files in folder name from query string')
    query_parser.add_argument('--res', action='store_true',
                              help='export query to .res files in folder name from query string')
    query_parser.add_argument('--ignore_warnings', action='store_true',
                              help='includes possibly bad structures')
    import_parser = subparsers.add_parser('import',
                                          help='import new structures in folder into database',
                                          parents=[global_parser, spatula_parser])
    rebuild_parser = subparsers.add_parser('rebuild',
                                           help='rebuild whole database.',
                                           parents=[spatula_parser])
    hull_parser = subparsers.add_parser('hull',
                                        help='create a convex hull from query results \
                                        (currently limited to binaries)',
                                        parents=[global_parser, structure_parser,
                                                 material_parser, plot_parser])
    voltage_parser = subparsers.add_parser('voltage',
                                           help='plot a voltage curve from query results \
                                           (currently limited to binaries)',
                                           parents=[global_parser, structure_parser,
                                                    material_parser, plot_parser])
    swaps_parser = subparsers.add_parser('swaps',
                                         help='perform atomic swaps on query results',
                                         parents=[global_parser, collection_parser,
                                                  structure_parser, material_parser])
    swaps_parser.add_argument('-s', '--swap', type=str, nargs='+',
                              help='swap all atoms in structures from a query from the first n-1 \
                                    species to the nth, e.g. --swaps N P As will swap all N, P \
                                    atoms for As. Uses the same macros  as --composition.')
    polish_parser = subparsers.add_parser('polish',
                                          help='re-relax a series of structures with \
                                          new parameters.',
                                          parents=[global_parser, collection_parser,
                                                   structure_parser])
    # parser.add_argument('--dis', action='store_true',
                        # help='smear hull with local stoichiometry')
    # parser.add_argument('--write_pressure', nargs='+', type=str,
                        # help=('pressure to add to new cell file, either one float' +
                              # 'for isotropic or 6 floats for anisotropic.'))
    args = parser.parse_args()
    if vars(args)['include_oqmd']:
        exit('--include_oqmd is currently disabled, please try again soon...')
    matador = Matador(args, argstr=argv[1:])
