# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the "matador" command line. """


from matador.query import DBQuery
from matador.hull import QueryConvexHull
from matador.utils.print_utils import print_failure, print_warning, print_notify
from matador.utils.db_utils import make_connection_to_collection, load_custom_settings


class MatadorCommandLine(object):
    """ Class that implements the command-line interface to a MongoDB
    structure repository.

    """
    def __init__(self, *args, **kwargs):
        """ Initialise the query with command line arguments and return
        results.

        """
        # read args
        self.kwargs = kwargs
        self.args = vars(args[0])
        self.args['testing'] = self.kwargs.get('testing')
        self.argstr = kwargs.get('argstr')

        file_exts = ['cell', 'res', 'pdb', 'markdown', 'latex', 'param', 'xsf']
        self.export = any([self.args.get(ext) for ext in file_exts])

        if self.args['subcmd'] != 'import':
            self.settings = load_custom_settings(config_fname=self.args.get('config'))
            result = make_connection_to_collection(self.args.get('db'),
                                                   check_collection=True,
                                                   mongo_settings=self.settings)
            self.client, self.db, self.collections = result

        if self.args['subcmd'] == 'stats':
            self.stats()

        try:
            if self.args['subcmd'] == 'import':
                from matador.db import Spatula
                self.importer = Spatula(self.args)

            if self.args['subcmd'] == 'query':
                self.query = DBQuery(self.client, self.collections, **self.args)
                self.cursor = self.query.cursor

            if self.args['subcmd'] == 'swaps':
                from matador.swaps import AtomicSwapper
                self.query = DBQuery(self.client, self.collections, **self.args)
                if self.args.get('hull_cutoff') is not None:
                    self.hull = QueryConvexHull(self.query, **self.args)
                    self.swapper = AtomicSwapper(self.hull.hull_cursor, self.args)
                else:
                    self.swapper = AtomicSwapper(self.query.cursor, self.args)
                self.cursor = self.swapper.cursor

            if self.args['subcmd'] == 'refine':
                from matador.db import Refiner
                self.query = DBQuery(self.client, self.collections, **self.args)
                if self.args.get('hull_cutoff') is not None:
                    self.hull = QueryConvexHull(self.query, **self.args)
                    self.refiner = Refiner(self.hull.cursor, self.query.repo, **self.args)
                else:
                    self.refiner = Refiner(self.query.cursor, self.query.repo, **self.args)

            if self.args['subcmd'] == 'pdffit':
                self.query = DBQuery(self.client, self.collections, **self.args)
                self.cursor = list(self.query.cursor)
                if self.args.get('hull_cutoff') is not None:
                    self.hull = QueryConvexHull(self.query, **self.args)
                    self.cursor = self.hull.hull_cursor
                    self.top = len(self.cursor)
                if self.args.get('top') is not None:
                    self.top = self.args.get('top')
                if not self.cursor[:self.top]:
                    print_notify('Performing PDF fit for ' + str(len(self.cursor[:self.top])) + ' structures.')
                    from matador.plugins.pdffit.pdffit import PDFFitter
                    self.pdffit = PDFFitter(self.cursor[:self.top], **self.args)
                    try:
                        self.pdffit.spawn()
                    except (KeyboardInterrupt, RuntimeError, SystemExit) as oops:
                        raise oops('Exiting top-level...')
                else:
                    exit('No structure match query.')

            if self.args['subcmd'] == 'hull' or self.args['subcmd'] == 'voltage':
                self.query = DBQuery(self.client, self.collections, **self.args)
                self.hull = QueryConvexHull(self.query, **self.args)
                self.cursor = self.hull.hull_cursor

            if self.args['subcmd'] == 'changes':
                from matador.db import DatabaseChanges
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
                from matador.hull.hull_diff import diff_hulls
                if self.args.get('compare') is None:
                    exit('Please specify which hulls to query with --compare.')
                diff_hulls(self.client, self.collections, **self.args)

            # perform any extra filtration
            if self.args.get('filter'):
                from matador.utils.cursor_utils import filter_cursor
                self.cursor = filter_cursor(self.cursor, self.args.get('filter'), self.args.get('values'))

            if self.export and self.cursor:
                from matador.export import query2files
                if self.args.get('write_n') is not None:
                    self.cursor = [doc for doc in self.cursor if len(doc['stoichiometry']) == self.args.get('write_n')]
                if len(self.cursor) < 1:
                    print_failure('No structures left to export.')
                query2files(self.cursor, self.args, argstr=self.argstr)

            if self.args.get('view'):
                from matador.utils.viz_utils import viz
                if self.args.get('top') is None:
                    self.top = len(self.cursor)
                else:
                    self.top = self.args.get('top')
                if len(self.cursor[:self.top]) > 10:
                    from time import sleep
                    print_warning('WARNING: opening {} files with ase-gui...'.format(len(self.cursor)))
                    print_warning('Please kill script within 3 seconds if undesired...')
                    sleep(3)
                if len(self.cursor[:self.top]) > 20:
                    print_failure(
                        'You will literally be opening that many windows, ' +
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
            print('Database last modified on', report['last_modified'], 'with matador', report['version'] + '.')
        except Exception:
            print_warning('Failed to print database report: spatula is probably running!')

    def stats(self):
        """ Print some useful stats about the database. """
        if self.args.get('list'):
            print_notify(str(len(self.db.collection_names())) + ' collections found in database:\n')
            collstats_list = []
            for name in self.db.collection_names():
                collstats_list.append(self.db.command('collstats', name))
                collstats_list[-1]['name'] = name
            collstats_list = sorted(collstats_list, key=lambda k: k['count'], reverse=True)
            print("\t{:^20}\t{:^20}".format('Name', 'Number of structures'))
            for collection in collstats_list:
                if not collection['name'].startswith('__'):
                    print("\t{:<20}\t{:>20d}".format(collection['name'], collection['count']))
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
                    exit('I cannot delete a collection that\'s name does not start with '
                         'your username, {}'.format(user))
                stats = self.db.command('collstats', target)
                answer = input('Are you sure you want to delete collection {} containing {} '
                               'structures? [y/n]\n'.format(target, stats['count']))
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


def main(testing=False):
    """ Parse all user args and construct a MatadorCommandLine object. """

    import argparse
    import os
    from sys import argv
    from matador import __version__

    parser = argparse.ArgumentParser(
        prog='matador',
        description='MATerial and Atomic Database Of Refined structures.',
        epilog='Written and maintained by Matthew Evans (me388@cam.ac.uk) 2016-2017, version {}.'
        .format(__version__.strip()))
    parser.add_argument('--version', action='version', version='matador version ' + __version__ + '.')

    # define subparsers for subcommands
    subparsers = parser.add_subparsers(title='subcommands', description='valid sub-commands', dest='subcmd')

    # define parent parser for global arguments
    global_flags = argparse.ArgumentParser(add_help=False)

    # common arguments to all subcommands
    global_flags.add_argument('--db', nargs='+', help='choose which collection to query')
    global_flags.add_argument('--debug', action='store_true', help='enable debug printing throughout code.')
    global_flags.add_argument('-conf', '--config', type=str,
                              help='specify custom location of matador config file.'
                                   '(DEFAULT: $MATADOR_ROOT/config/matador_conf.json)')
    global_flags.add_argument('--devel', action='store_true', help='test devel code.')
    global_flags.add_argument('--profile', action='store_true', help='run code profiler.')
    global_flags.add_argument('-q', '--quiet', action='store_true', help='redirect most output to /dev/null.')

    # define all other flags by group
    structure_flags = argparse.ArgumentParser(add_help=False)
    structure_flags.add_argument('-c', '--composition', type=str, nargs='+',
                                 help='find all structures containing exclusively the given '
                                      'elements, e.g. LiSi. Macros defined for groups [I]-[VII] '
                                      '[Tran] [Lan] and [Act], used with square brackets.')
    structure_flags.add_argument('-int', '--intersection', action='store_true',
                                 help='query the intersection of compositions instead of the union '
                                      'e.g. -c LiSnS -int queries Li, Sn, S, LiSn, LiS & LiSnS.')
    structure_flags.add_argument('-n', '--num_species', type=int,
                                 help='find all structures containing a certain number of species.')
    structure_flags.add_argument('-f', '--formula', type=str, nargs='+',
                                 help='query a particular chemical formula, e.g. GeTeSi3')
    structure_flags.add_argument('-i', '--id', type=str, nargs='+',
                                 help='specify a particular structure by its text_id')
    structure_flags.add_argument('-ac', '--calc_match', action='store_true',
                                 help='display calculations of the same accuracy as specified id')
    structure_flags.add_argument('-kpttol', '--kpoint-tolerance', type=float,
                                 help='kpoint tolerance for calculation matches (DEFAULT: +/- 0.01 1/Å)')
    structure_flags.add_argument('-z', '--num_fu', type=int,
                                 help='query a calculations with more than n formula units')
    structure_flags.add_argument('-sg', '--space_group', help='query a particular space group')
    structure_flags.add_argument('-u', '--uniq', type=float, nargs='?', const=0.1,
                                 help='float, return only unique structures (filtered by PDF '
                                      'overlap), to this tolerance (DEFAULT: 0.1)')
    structure_flags.add_argument('-p', '--pressure', type=float,
                                 help='specify an isotropic external pressure to search for, e.g. 10 (GPa)')
    structure_flags.add_argument('-pf', '--partial-formula', action='store_true',
                                 help='stoichiometry/composition queries will include other unspecified species, e.g. '
                                      '-pf search for Li will query any structure containing Li, not just pure Li.')
    structure_flags.add_argument('--tags', nargs='+', type=str, help=('search for manual tags'))
    structure_flags.add_argument('--doi', type=str, help=('search for DOI in format xxxx/xxxx'))
    structure_flags.add_argument('-icsd', '--icsd', type=int, const=0, nargs='?', help=('search for an ICSD CollCode'))
    structure_flags.add_argument('-ss', '--src_str', type=str,
                                 help=('search for a string inside the structure sources'))
    structure_flags.add_argument('-root', '--root_src', type=str,
                                 help=('search for a root_source string of the structure'))
    structure_flags.add_argument('-encap', '--encapsulated', action='store_true',
                                 help='query only structures encapsulated in a carbon nanotube.')
    structure_flags.add_argument('-cntr', '--cnt_radius', type=float,
                                 help='specify the radius of the encapsulating nanotube to within 0.01 Å')
    structure_flags.add_argument('-cntv', '--cnt_vector', type=int, nargs='+',
                                 help='specify the chiral vector of the encapsulating nanotube')
    structure_flags.add_argument('-ecut', '--cutoff', type=float, nargs='+',
                                 help='specify the min. and optionally max. planewave cutoff.')
    structure_flags.add_argument('-geom', '--geom_force_tol', type=float, nargs='+',
                                 help='force tolerance in eV/Å to query for calc matches.')
    structure_flags.add_argument('--sedc', type=str, help='specify the dispersion correction scheme, e.g. TS or null.')
    structure_flags.add_argument('-xc', '--xc_functional', type=str,
                                 help='specify an xc-functional to query (case-insensitive).')
    structure_flags.add_argument('-kpts', '--mp_spacing', type=float,
                                 help='specify an MP grid spacing in 2π/Å units, e.g. 0.05, will return all values '
                                      'structures with value within --kpt_tol')
    structure_flags.add_argument('--spin', type=str,
                                 help='specifiy whether to query non-spin-polarized (0) calcs or spin polarized calcs '
                                      '(!=1), or lump them both together with `any`')
    structure_flags.add_argument('--loose', action='store_true',
                                 help='loosely matches with calc_match, i.e. only matches pspot and xc_functional')
    structure_flags.add_argument('--ignore_warnings', action='store_true', help='includes possibly bad structures')
    structure_flags.add_argument('--filter', type=str,
                                 help='specify a simple float field to filter. Requires --values')
    structure_flags.add_argument('--values', nargs='+', type=float,
                                 help='specify the minimum floats, or [min, max] values of field')

    material_flags = argparse.ArgumentParser(add_help=False)
    material_flags.add_argument('-hc', '--hull_cutoff', type=float,
                                help='return only structures within a certain distance from hull in eV/atom')
    material_flags.add_argument('-lc', '--label_cutoff', type=float,
                                help='label only structures within a certain distance from hull in eV/atom')
    material_flags.add_argument('-hT', '--hull_temp', type=float,
                                help='return only structures within a certain distance from hull in K')
    material_flags.add_argument('--biggest', action='store_true',
                                help='use the largest subset of structures to create a hull')
    material_flags.add_argument('--volume', action='store_true',
                                help='plot a volume curve from convex hull (currently limited to binaries)')
    material_flags.add_argument('--chempots', type=float, nargs='+',
                                help='manually specify chem pots as enthalpy per atom for a rough hull.')

    plot_flags = argparse.ArgumentParser(add_help=False)
    plot_flags.add_argument('--pdf', action='store_true', help='save pdf rather than showing plot in X')
    plot_flags.add_argument('--png', action='store_true', help='save png rather than showing plot in X')
    plot_flags.add_argument('--csv', action='store_true', help='save plotting data to separate csv files')
    plot_flags.add_argument('--labels', action='store_true', help='label hull plots')
    plot_flags.add_argument('--svg', action='store_true', help='save svg rather than showing plot in X')
    plot_flags.add_argument('--subplot', action='store_true', help='plot combined hull and voltage graph')
    plot_flags.add_argument('--no_plot', action='store_true', help='suppress plotting')
    plot_flags.add_argument('--capmap', action='store_true', help='plot heat map of gravimetric capacity')
    plot_flags.add_argument('--sampmap', action='store_true', help='plot heat map of concentration sampling')
    plot_flags.add_argument('--efmap', action='store_true', help='plot heat map of formation energy')
    plot_flags.add_argument('--pathways', action='store_true',
                            help='plot line from stable B_x C_y to pure A in ABC ternary.')
    plot_flags.add_argument('--expt', type=str, help='enter experimental voltage curve .csv file for plotting.')
    plot_flags.add_argument('--expt_label', type=str, help='label for experimental data on voltage curve.')

    spatula_flags = argparse.ArgumentParser(add_help=False)
    spatula_flags.add_argument('-d', '--dryrun', action='store_true',
                               help='run the importer without connecting to the database')
    spatula_flags.add_argument('-v', '--verbosity', action='count', help='enable verbose output')
    spatula_flags.add_argument('-f', '--force', action='store_true', help='override main database protection')
    spatula_flags.add_argument('-t', '--tags', nargs='+', type=str, help='set user tags, e.g. nanotube, project name')
    spatula_flags.add_argument('-s', '--scan', action='store_true',
                               help='only scan the database for new structures, do not dictify')

    changes_flags = argparse.ArgumentParser(add_help=False)
    changes_flags.add_argument('-c', '--changeset', type=int, help='changeset number to query')
    changes_flags.add_argument('-r', '--revert', type=int, help='revert database to specified changeset')
    changes_flags.add_argument('-u', '--undo', action='store_true', help='undo changeset')

    collection_flags = argparse.ArgumentParser(add_help=False)
    collection_flags.add_argument('--to', type=str, help='the text_id of a structure with the desired parameters')
    collection_flags.add_argument('--with', type=str,
                                  help=('the seedname (must be within pwd) of cell and param ' +
                                        'files to use for swaps'))
    collection_flags.add_argument('--prefix', type=str,
                                  help='add a prefix to all file names to write out (auto-appended with an underscore')

    query_flags = argparse.ArgumentParser(add_help=False)

    query_flags.add_argument('-s', '--summary', action='store_true',
                             help='show only the ground state for each stoichiometry.')
    query_flags.add_argument('-t', '--top', type=int, help='number of structures to show/write (DEFAULT: 10)')
    query_flags.add_argument('-dE', '--delta_E', type=float,
                             help='maximum distance from ground state structure to show/write in eV/atom')
    query_flags.add_argument('-d', '--details', action='store_true',
                             help='show as much detail about calculation as possible')
    query_flags.add_argument('-pa', '--per_atom', action='store_true', help='show quantities per atom not per fu.')
    query_flags.add_argument('-dt', '--time', type=int, help='query only structures added before this time in days')
    query_flags.add_argument('--source', action='store_true',
                             help='print filenames from which structures were wrangled')
    query_flags.add_argument('-v', '--view', action='store_true',
                             help='quickly view a structure/structures with ase-gui')
    query_flags.add_argument('--cell', action='store_true',
                             help='export query to .cell files in folder name from query string')
    query_flags.add_argument('--param', action='store_true',
                             help='export query to .param files in folder name from query string')
    query_flags.add_argument('--res', action='store_true',
                             help='export query to .res files in folder name from query string')
    query_flags.add_argument('--pdb', action='store_true',
                             help='export query to .pdb files in folder name from query string')
    query_flags.add_argument('--xsf', action='store_true',
                             help='export query to .xsf files in folder name from query string')
    query_flags.add_argument('--markdown', action='store_true', help='export query summary to a markdown file')
    query_flags.add_argument('--latex', action='store_true', help='export query summary to a LaTeX table')
    query_flags.add_argument('--write_n', type=int, help='export only those structures with n species')

    swap_flags = argparse.ArgumentParser(add_help=False)
    swap_flags.add_argument('-sw', '--swap', type=str, nargs='+',
                            help='swap all atoms in structures from a query from the first n-1 species to the nth, '
                                 'e.g. -sw NAs will swap all N to As, -sw NAs:LiNa will swap all N to As, and all Li '
                                 'to Na, and -sw [V]As:[Li,K,Rb]Na will swap all group V elements to As and all of Li,'
                                 'K and Rb to Na.')
    diff_flags = argparse.ArgumentParser(add_help=False)
    diff_flags.add_argument('-cmp', '--compare', type=str, nargs='+',
                            help='diff phase diagrams between two different times, in standard time format, '
                                 'e.g. `--compare 1y2m5d3h` will compare the present hull with that of 1 year, 2 '
                                 'months, 5 days and 3 hours ago, and `--compare 3d 2d` will compare three days ago '
                                 'to two days ago.')
    pdffit_flags = argparse.ArgumentParser(add_help=False)
    pdffit_flags.add_argument('-file', '--file', type=str, help='experimental input file to fit structures to.')
    pdffit_flags.add_argument('-min', '--xmin', type=float,
                              help='minimum value to compute the PDF (DEFAULT: 1 Angstrom)')
    pdffit_flags.add_argument('-max', '--xmax', type=float,
                              help='maximum value to compute the PDF (DEFAULT: 50 Angstrom')
    pdffit_flags.add_argument('-dx', '--dx', type=float, help='spacing to compute PDF at')
    pdffit_flags.add_argument('-2', '--two_phase', type=float, help='fit two phases to experimental PDF')
    pdffit_flags.add_argument('-np', '--num_processes', type=int, help='number of concurrent fits to perform.')

    refine_flags = argparse.ArgumentParser(add_help=False)
    refine_flags.add_argument('-task', '--task', type=str, help='refine subtask to perform: options are spg or sub')
    refine_flags.add_argument('-mode', '--mode', type=str,
                              help='mode of refinement: options are display, set and overwrite')
    refine_flags.add_argument('-symprec', '--symprec', type=float, help='spglib symmetry precision for refinement')
    refine_flags.add_argument('--new_tag', type=str, help='new tag to add to structures in query')
    refine_flags.add_argument('--new_doi', type=str, help='new doi to add to structures in query')

    stats_flags = argparse.ArgumentParser(add_help=False)
    stats_flags.add_argument('-l', '--list', action='store_true', help='list all collections, their sizes, and owners')
    stats_flags.add_argument('--delete', action='store_true', help='try to delete collection specified by --db')

    # define subcommand parsers and their arguments

    # matador stats
    subparsers.add_parser('stats', help='print some stats about the database.', parents=[global_flags, stats_flags])

    # matador query
    subparsers.add_parser('query',
                          help='query and extract structures from the database',
                          parents=[global_flags, query_flags, structure_flags])

    # matador import
    subparsers.add_parser('import',
                          help='import new structures in folder into database',
                          parents=[global_flags, spatula_flags])

    # matador pdffit
    subparsers.add_parser('pdffit',
                          help='provide experimental .gr file and fit to calculated PDF of structures in query',
                          parents=[global_flags, query_flags, material_flags,
                                   structure_flags, pdffit_flags])

    # matador hull
    subparsers.add_parser('hull',
                          help='create a convex hull from query results (currently limited to binaries and ternaries)',
                          parents=[global_flags, structure_flags,
                                   material_flags, plot_flags, query_flags])

    # matador voltage
    subparsers.add_parser('voltage',
                          help='plot a voltage curve from query results (currently limited to binaries and ternaries)',
                          parents=[global_flags, structure_flags,
                                   material_flags, plot_flags, query_flags])

    # matador changes
    subparsers.add_parser('changes',
                          help='view database changelog or undo additions to database (NB: not deletions!)',
                          parents=[global_flags, changes_flags])

    # matador hulldiff
    subparsers.add_parser('hulldiff',
                          help='diff two convex hulls with the --compare flag.',
                          parents=[global_flags, structure_flags,
                                   material_flags, plot_flags, query_flags, diff_flags])

    # matador swaps
    subparsers.add_parser('swaps',
                          help='perform atomic swaps on query results',
                          parents=[global_flags, collection_flags, query_flags,
                                   structure_flags, material_flags, swap_flags])

    # matador refine
    subparsers.add_parser('refine',
                          help='update structures in the database according to specified --task',
                          parents=[global_flags, query_flags, structure_flags,
                                   refine_flags, material_flags])

    parsed_args = parser.parse_args()

    # check for inconsistent argument combinations
    if vars(parsed_args).get('intersection') and vars(parsed_args).get('composition') is None:
        raise SystemExit('--intersection requires --composition.')
    if vars(parsed_args).get('subcmd') == 'stats' and vars(parsed_args).get('list') and vars(parsed_args).get(
            'delete'):
        raise SystemExit('Cannot use -l/--list and --delete')
    # if vars(parsed_args).get('formula') and vars(parsed_args).get('composition'):
    # raise SystemExit('Cannot use -f/--formula and -c/--composition together.')
    if vars(parsed_args).get('filter') and vars(parsed_args).get('values') is None:
        raise SystemExit('--filter requires --values.')
    if vars(parsed_args).get('values') and vars(parsed_args).get('filter') is None:
        print('Ignoring redundant supplied values...')
    if vars(parsed_args).get('subcmd') == 'hull' and vars(parsed_args).get('composition') is None:
        raise SystemExit('hull requires --composition')
    if vars(parsed_args).get('subcmd') == 'pdffit':
        if vars(parsed_args).get('file') is None:
            raise SystemExit('pdffit requires specified --file, exiting...')
        if not os.path.isfile(vars(parsed_args).get('file')):
            raise SystemExit('specified --file does not exist, exiting...')
    if vars(parsed_args).get('hull_cutoff') and vars(parsed_args).get('hull_temp'):
        raise SystemExit('hull_cutoff and hull_temp both specified, exiting...')
    if vars(parsed_args).get('calc_match') and vars(parsed_args).get('id') is None:
        raise SystemExit('calc_match requires specification of a text_id with -i, exiting...')
    if vars(parsed_args).get('profile'):
        import cProfile
        import pstats
        from sys import version_info
        profiler = cProfile.Profile()
        profiler.enable()

    MatadorCommandLine(parsed_args, argstr=argv[1:], testing=testing)

    if vars(parsed_args).get('profile'):
        profiler.disable()
        fname = 'matador-{}-{}.{}.{}'.format(__version__, version_info.major, version_info.minor, version_info.micro)
        profiler.dump_stats(fname + '.prof')
        with open(fname + '.pstats', 'w') as fp:
            stats = pstats.Stats(profiler, stream=fp).sort_stats('cumulative')
            stats.print_stats()


if __name__ == '__main__':
    main()
