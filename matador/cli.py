# coding: utf-8
""" This file implements the "matador" command line, as initialised by
`../bin/matador`.

"""

from matador.query import DBQuery
from matador.hull import QueryConvexHull
from matador.utils.print_utils import print_failure, print_warning, print_notify
from matador.utils.db_utils import make_connection_to_collection, load_custom_settings


class MatadorCommandLine(object):
    """ Class that implements the command-line interface
    to a MongoDB structure repository.
    """

    def __init__(self, *args, **kwargs):
        """ Initialise the query with command line
        arguments and return results.
        """
        # read args
        self.kwargs = kwargs
        self.args = vars(args[0])
        self.argstr = kwargs.get('argstr')

        file_exts = ['cell', 'res', 'pdb', 'markdown', 'latex', 'param', 'xsf']
        self.export = any([self.args.get(ext) for ext in file_exts])

        if self.args['subcmd'] != 'import':
            self.settings = load_custom_settings(
                config_fname=self.args.get('config'))
            result = make_connection_to_collection(self.args.get('db'),
                                                   check_collection=True,
                                                   mongo_settings=self.settings)
            self.client, self.db, self.collections = result

        if self.args['subcmd'] == 'stats':
            self.stats()

        try:

            if self.args['subcmd'] == 'import':
                from matador.importer import Spatula
                self.importer = Spatula(self.args)

            if self.args['subcmd'] == 'query':
                self.query = DBQuery(
                    self.client, self.collections, **self.args)
                self.cursor = self.query.cursor

            if self.args['subcmd'] in ['swaps']:
                from matador.polish import Polisher
                self.query = DBQuery(
                    self.client, self.collections, **self.args)
                if self.args.get('hull_cutoff') is not None:
                    self.hull = QueryConvexHull(self.query, **self.args)
                    self.polisher = Polisher(self.hull.hull_cursor, self.args)
                else:
                    self.polisher = Polisher(self.query.cursor, self.args)
                self.cursor = self.polisher.cursor

            if self.args['subcmd'] == 'refine':
                from matador.refine import Refiner
                self.query = DBQuery(
                    self.client, self.collections, **self.args)
                if self.args.get('hull_cutoff') is not None:
                    self.hull = QueryConvexHull(self.query, **self.args)
                    self.refiner = Refiner(
                        self.hull.cursor, self.query.repo, **self.args)
                else:
                    self.refiner = Refiner(
                        self.query.cursor, self.query.repo, **self.args)

            if self.args['subcmd'] == 'pdffit':
                self.query = DBQuery(
                    self.client, self.collections, **self.args)
                self.cursor = list(self.query.cursor)
                if self.args.get('hull_cutoff') is not None:
                    self.hull = QueryConvexHull(self.query, **self.args)
                    self.cursor = self.hull.hull_cursor
                    self.top = len(self.cursor)
                if self.args.get('top') is not None:
                    self.top = self.args.get('top')
                if not self.cursor[:self.top]:
                    print_notify('Performing PDF fit for ' +
                                 str(len(self.cursor[:self.top])) +
                                 ' structures.')
                    from matador.pdffit import PDFFitter
                    self.pdffit = PDFFitter(
                        self.cursor[:self.top], **self.args)
                    try:
                        self.pdffit.spawn()
                    except(KeyboardInterrupt, RuntimeError, SystemExit) as oops:
                        raise oops('Exiting top-level...')
                else:
                    exit('No structure match query.')

            if self.args['subcmd'] == 'hull' or self.args['subcmd'] == 'voltage':
                self.query = DBQuery(
                    self.client, self.collections, **self.args)
                self.hull = QueryConvexHull(self.query, **self.args)
                self.cursor = self.hull.hull_cursor

            if self.args['subcmd'] == 'changes':
                from matador.changes import DatabaseChanges
                if len(self.collections) != 1:
                    exit('Cannot view changes of more than one collection at once.')
                if self.args.get('undo'):
                    action = 'undo'
                else:
                    action = 'view'
                changeset = self.args.get('changeset')
                if changeset is None:
                    changeset = 0
                DatabaseChanges([key for key in self.collections][0],
                                changeset_ind=changeset,
                                action=action,
                                mongo_settings=self.settings)

            if self.args['subcmd'] == 'hulldiff':
                from matador.hull_diff import diff_hulls
                if self.args.get('compare') is None:
                    exit('Please specify which hulls to query with --compare.')
                diff_hulls(self.client, self.collections, **self.args)

            # perform any extra filtration
            if self.args.get('filter'):
                from matador.utils.cursor_utils import filter_cursor
                self.cursor = filter_cursor(self.cursor,
                                            self.args.get('filter'),
                                            self.args.get('values'))

            if self.export and self.cursor:
                from matador.export import query2files
                if self.args.get('write_n') is not None:
                    self.cursor = [doc for doc in self.cursor if len(
                        doc['stoichiometry']) == self.args.get('write_n')]
                if len(self.cursor) < 1:
                    print_failure('No structures left to export.')
                query2files(self.cursor, self.args, argstr=self.argstr)

            if self.args.get('view'):
                from matador.viz import viz
                if self.args.get('top') is None:
                    self.top = len(self.cursor)
                else:
                    self.top = self.args.get('top')
                if len(self.cursor[:self.top]) > 10:
                    from time import sleep
                    print_warning(
                        'WARNING: opening {} files with ase-gui...'.format(len(self.cursor)))
                    print_warning(
                        'Please kill script within 3 seconds if undesired...')
                    sleep(3)
                if len(self.cursor[:self.top]) > 20:
                    print_failure('You will literally be opening that many windows, ' +
                                  'I\'ll give you another 5 seconds to reconsider...')

                    sleep(5)
                    print_notify('It\'s your funeral...')
                    sleep(1)
                for doc in self.cursor[:self.top]:
                    viz(doc)

            if self.args.get('subcmd') != 'import':
                self.client.close()

        except SystemExit as oops:
            print(oops)
            print('Trying to nicely close connection...')
            try:
                self.client.close()
            except AttributeError:
                pass

    def print_report(self):
        """ Print spatula report on current database. """
        try:
            report = self.report.find_one()
            print('Database last modified on', report['last_modified'], 'with matador',
                  report['version'] + '.')
        except Exception:
            print_warning('Failed to print database report: spatula is probably running!')

    def stats(self):
        """ Print some useful stats about the database. """
        if self.args.get('list'):
            print_notify(str(len(self.db.collection_names())) +
                         ' collections found in database:\n')
            collstats_list = []
            for name in self.db.collection_names():
                collstats_list.append(self.db.command('collstats', name))
                collstats_list[-1]['name'] = name
            collstats_list = sorted(
                collstats_list, key=lambda k: k['count'], reverse=True)
            print("\t{:^20}\t{:^20}".format('Name', 'Number of structures'))
            for collection in collstats_list:
                if not collection['name'].startswith('__'):
                    print("\t{:<20}\t{:>20d}".format(
                        collection['name'], collection['count']))
            print('\n')
        elif self.args.get('delete'):
            target = self.args.get('db')
            if isinstance(target, list) and len(target) == 1:
                target = target[0]
            else:
                exit('I will only delete one collection at a time...')
            if target is None:
                exit('Please specify a collection to delete.')
            elif target not in self.db.collection_names():
                exit('No collection named {} was found'.format(target))
            else:
                from getpass import getuser
                user = getuser()
                if user not in target:
                    exit('I cannot delete a collection that\'s name does not start with \
                         your username, {}'.format(user))
                stats = self.db.command('collstats', target)
                answer = input('Are you sure you want to delete collection {} containing {} \
                               structures? [y/n])'.format(target, stats['count']))
                if answer.lower() == 'y':
                    if target == 'repo':
                        exit('I\'m sorry Dave, I\'m afraid I can\'t do that...')
                    else:
                        print('Deleting collection {}...'.format(target))
                        self.db[target].drop()
                        print('and its changelog...')
                        self.db['__changelog_{}'.format(target)].drop()
                else:
                    exit('Nevermind then!')
        else:
            comp_list = dict()
            stats_dict = dict()
            stats_dict['count'] = 0
            stats_dict['avgObjSize'] = 0
            stats_dict['storageSize'] = 0
            stats_dict['totalIndexSize'] = 0
            for collection in self.collections:
                db_stats_dict = self.db.command('collstats', collection)
                stats_dict['count'] += db_stats_dict['count']
                stats_dict['avgObjSize'] += db_stats_dict['avgObjSize']
                stats_dict['storageSize'] += db_stats_dict['storageSize']
                stats_dict['totalIndexSize'] += db_stats_dict['totalIndexSize']
            print(("The collection(s) queried in {} contain {} structures at {:.1f} kB each "
                   "totalling {} when padding and indices are included.")
                  .format(self.db.name,
                          stats_dict['count'],
                          stats_dict['avgObjSize'] /
                          (1024 * len(self.collections)),
                          (stats_dict['totalIndexSize'] + stats_dict['storageSize']) /
                          (1024**2)))
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
            comp_list = list(zip(keys, vals))
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
            except ImportError:
                exit('Pyascii graph missing; not printing detailed stats.')

            graph = Pyasciigraph(line_length=80, multivalue=False)
            thresholds = {int(stats_dict['count'] / 40): Gre,
                          int(stats_dict['count'] / 10): Blu,
                          int(stats_dict['count'] / 4): Red}
            data = hcolor(comp_list, thresholds)
            for line in graph.graph(label=None, data=data):
                print(line)
            print('\n')
