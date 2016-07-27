#!/usr/bin/python
# coding: utf-8
""" This file implements all queries to the database,
including parsing user inputs, displaying results
and calling other functionality. """
from __future__ import print_function
# matador modules
from query import DBQuery
from print_utils import print_failure, print_warning
from polish import Polisher
# import external libraries
import pymongo as pm
# import standard library
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
            from spatula import Spatula
            self.importer = Spatula(self.args)
        if self.args['subcmd'] == 'query':
            self.query = DBQuery(self.client, self.collections, self.args)
        if self.args['subcmd'] == 'swaps':
            self.query = DBQuery(self.client, self.collections, self.args)
            if self.args['hull_cutoff'] is not None:
                from hull import QueryConvexHull
                self.hull = QueryConvexHull(self.query, self.args)
                self.swaps = Polisher(self.hull.hull_cursor, self.args)
            else:
                self.swaps = Polisher(self.query.cursor, self.args)
        if self.args['subcmd'] == 'polish':
            self.query = DBQuery(self.client, self.collections, self.args)
            if self.args['hull_cutoff'] is not None:
                from hull import QueryConvexHull
                self.hull = QueryConvexHull(self.query, self.args)
                self.polish = Polisher(self.hull.hull_cursor, self.args)
            else:
                self.polish = Polisher(self.query.cursor, self.args)
        if self.args['subcmd'] == 'hull' or self.args['subcmd'] == 'voltage':
            from hull import QueryConvexHull
            self.query = DBQuery(self.client, self.collections, self.args)
            self.hull = QueryConvexHull(self.query, self.args)

    def print_report(self):
        """ Print spatula report on current database. """
        try:
            report = self.report.find_one()
            print('Database last modified on', report['last_modified'], 'with spatula',
                  report['version'], 'changeset (' + report['git_hash'] + ').')
        except:
            print_warning('Failed to print database report: spatula is probably running!')

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
            print_failure('Pyascii graph missing; not printing detailed stats.')
            exit()
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
            print_failure('No structures found.')
            exit()
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
    global_flags = argparse.ArgumentParser(add_help=False)
    # common arguments to all subcommands
    global_flags.add_argument('--db',
                              help='choose which databases to query',
                              nargs='+')
    # define all other flags by group
    # define flags for query contents
    structure_flags = argparse.ArgumentParser(add_help=False)
    structure_flags.add_argument('-c', '--composition', type=str, nargs='+',
                                 help='find all structures containing exclusively the given elements, \
                                       e.g. LiSi. Macros defined for groups [I]-[VII], [Tran] \
                                       [Lan] and [Act], used with square brackets.')
    structure_flags.add_argument('-int', '--intersection', action='store_true',
                                 help='query the intersection of compositions instead of the union \
                                       e.g. -c LiSnS -int queries Li, Sn, S, LiSn, LiS & LiSnS.')
    structure_flags.add_argument('-n', '--num_species', type=int, nargs='+',
                                 help='find all structures containing a certain \
                                       number of species.')
    structure_flags.add_argument('-f', '--formula', type=str, nargs='+',
                                 help='query a particular chemical formula, e.g. GeTeSi3')
    structure_flags.add_argument('-i', '--id', type=str, nargs='+',
                                 help='specify a particular structure by its text_id')
    structure_flags.add_argument('-ac', '--calc-match', action='store_true',
                                 help='display calculations of the same accuracy as specified id')
    structure_flags.add_argument('-z', '--num_fu', type=int,
                                 help='query a calculations with more than n formula units')
    structure_flags.add_argument('-sg', '--space_group',
                                 help='query a particular space group')
    structure_flags.add_argument('-p', '--pressure', type=float,
                                 help='specify an isotropic external pressure to search for \
                                       , e.g. 10 (GPa)')
    structure_flags.add_argument('-pf', '--partial-formula', action='store_true',
                                 help='stoichiometry/composition queries will include other \
                                       unspecified species, e.g. -pf search for Li will query \
                                       any structure containing Li, not just pure Li.')
    structure_flags.add_argument('--tags', nargs='+', type=str,
                                 help=('search for up to 3 manual tags at once'))
    structure_flags.add_argument('-encap', '--encapsulated', action='store_true',
                                 help='query only structures encapsulated in a carbon nanotube.')
    structure_flags.add_argument('-cntr', '--cnt_radius', type=float,
                                 help='specify the radius of the encapsulating nanotube \
                                       to within 0.01 A')
    structure_flags.add_argument('--loose', action='store_true',
                                 help='loosely matches with calc_match, i.e. only matches \
                                       pspot and xc_functional')
    structure_flags.add_argument('--ignore_warnings', action='store_true',
                                 help='includes possibly bad structures')
    # define material flags for hull/voltage arguments
    material_flags = argparse.ArgumentParser(add_help=False)
    material_flags.add_argument('-hc', '--hull_cutoff', type=float,
                                help='return only structures within a certain distance from hull \
                                      in eV/atom')
    material_flags.add_argument('-hT', '--hull_temp', type=float,
                                help='return only structures within a certain distance from hull \
                                      in K')
    material_flags.add_argument('--biggest', action='store_true',
                                help='use the largest subset of structures to create a hull')
    material_flags.add_argument('--volume', action='store_true',
                                help='plot a volume curve from convex hull\
                                      (currently limited to binaries)')
    material_flags.add_argument('--chempots', type=float, nargs='+',
                                help='manually specify chem pots as enthalpy per atom for \
                                      a rough hull.')
    # define flags for plotting output
    plot_flags = argparse.ArgumentParser(add_help=False)
    plot_flags.add_argument('--pdf', action='store_true',
                            help='save pdf rather than showing plot in X')
    plot_flags.add_argument('--subplot', action='store_true',
                            help='plot combined hull and voltage graph')
    plot_flags.add_argument('--bokeh', action='store_true',
                            help='plot using bokeh')
    plot_flags.add_argument('--no_plot', action='store_true',
                            help='suppress plotting')
    # define flags for spatula scraping
    spatula_flags = argparse.ArgumentParser(add_help=False)
    spatula_flags.add_argument('-d', '--dryrun', action='store_true',
                               help='run the importer without connecting to the database')
    spatula_flags.add_argument('-v', '--verbosity', action='count',
                               help='enable verbose output')
    spatula_flags.add_argument('-t', '--tags', nargs='+', type=str,
                               help='set user tags, e.g. nanotube, project name')
    spatula_flags.add_argument('--debug', action='store_true',
                               help='enable debug output to print every dict')
    spatula_flags.add_argument('-s', '--scan', action='store_true',
                               help='only scan the database for new structures, do not dictify')
    # define parser for output of new files for swaps/polishes
    collection_flags = argparse.ArgumentParser(add_help=False)
    collection_flags.add_argument('--to', type=str,
                                  help='the text_id of a structure with the desired parameters')
    collection_flags.add_argument('--with', type=str,
                                  help=('the seedname (must be within pwd) of cell and param ' +
                                        'files to use for polishing/swaps'))
    collection_flags.add_argument('--prefix', type=str,
                                  help='add a prefix to all file names to write out (auto-appended \
                                        with an underscore')
    query_flags = argparse.ArgumentParser(add_help=False)
    query_flags.add_argument('-s', '--summary', action='store_true',
                             help='show only the ground state for each stoichiometry.')
    query_flags.add_argument('-t', '--top', type=int,
                             help='number of structures to show (DEFAULT: 10)')
    query_flags.add_argument('-d', '--details', action='store_true',
                             help='show as much detail about calculation as possible')
    query_flags.add_argument('--source', action='store_true',
                             help='print filenames from which structures were wrangled')
    query_flags.add_argument('--cell', action='store_true',
                             help='export query to .cell files in folder name from query string')
    query_flags.add_argument('--res', action='store_true',
                             help='export query to .res files in folder name from query string')
    swap_flags = argparse.ArgumentParser(add_help=False)
    swap_flags.add_argument('-s', '--swap', type=str, nargs='+',
                            help='swap all atoms in structures from a query from the first n-1 \
                                  species to the nth, e.g. --swaps N P As will swap all N, P \
                                  atoms for As. Uses the same macros  as --composition.')
    # define subcommand parsers and their arguments
    stat_parser = subparsers.add_parser('stats',
                                        help='print some stats about the database.',
                                        parents=[global_flags])
    query_parser = subparsers.add_parser('query',
                                         help='query and extract structures from the database',
                                         parents=[global_flags, query_flags, structure_flags])
    import_parser = subparsers.add_parser('import',
                                          help='import new structures in folder into database',
                                          parents=[global_flags, spatula_flags])
    rebuild_parser = subparsers.add_parser('rebuild',
                                           help='rebuild whole database.',
                                           parents=[spatula_flags])
    hull_parser = subparsers.add_parser('hull',
                                        help='create a convex hull from query results \
                                        (currently limited to binaries)',
                                        parents=[global_flags, structure_flags,
                                                 material_flags, plot_flags, query_flags])
    voltage_parser = subparsers.add_parser('voltage',
                                           help='plot a voltage curve from query results \
                                           (currently limited to binaries)',
                                           parents=[global_flags, structure_flags,
                                                    material_flags, plot_flags, query_flags])
    swaps_parser = subparsers.add_parser('swaps',
                                         help='perform atomic swaps on query results',
                                         parents=[global_flags, collection_flags,
                                                  structure_flags, material_flags, swap_flags])
    polish_parser = subparsers.add_parser('polish',
                                          help='re-relax a series of structures with \
                                          new parameters.',
                                          parents=[global_flags, collection_flags,
                                                   structure_flags, material_flags,
                                                   query_flags])
    args = parser.parse_args()
    if vars(args).get('include_oqmd'):
        print_failure('--include_oqmd is currently disabled, please try again soon...')
        exit()
    if vars(args).get('intersection') and vars(args).get('composition') is None:
        print_failure('--intersection requires --composition.')
        exit()
    if vars(args).get('subcmd') == 'hull' and vars(args).get('formula') is not None:
        print_failure('hull not compatible with --formula, please use --composition')
        exit()
    if vars(args).get('hull_cutoff') and vars(args).get('hull_temp'):
        print_failure('hull_cutoff and hull_temp both specified, exiting...')
        exit()
    matador = Matador(args, argstr=argv[1:])
