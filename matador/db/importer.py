# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the base class Spatula that calls the scrapers
and interfaces with the MongoDB client.

"""


import os
import collections
import tempfile
import random
import datetime
import copy
import traceback as tb
import logging
import sys

import pymongo as pm

from matador.scrapers.castep_scrapers import castep2dict, param2dict, cell2dict
from matador.scrapers.castep_scrapers import res2dict
from matador.utils.cell_utils import calc_mp_spacing
from matador.utils.chem_utils import get_root_source
from matador.utils.cursor_utils import recursive_get
from matador.db import make_connection_to_collection
from matador.utils.db_utils import WORDS, NOUNS
from matador.config import load_custom_settings
from matador.utils.cursor_utils import loading_bar
from matador import __version__


class Spatula:
    """ The Spatula class implements methods to scrape folders and
    individual files for crystal structures and create a MongoDB
    document for each.

    Files types that can be read are:

        * CASTEP `.castep` output
        * SHELX (from airss.pl / pyAIRSS) `.res` output
        * CASTEP `.param`, `.cell` input

    This class will recursively scan directories from the cwd to
    find the files types above. Base filenames will be matched
    to prevent duplication of data from e.g. `.castep` and `.res`
    files. The following directory structures are recommended:

    - One `.res` file per structure and template `.cell` and
    `.param` files that provide all CASTEP parameters that structures
    in this folder were ran at.

    - One `.castep` file per structure, containing all information. If
    pseudopotential information is not present in the CASTEP file, this
    class will check for the corresponding `.usp` files and try to scrape
    those.

    """
    def __init__(self, *args, settings=None):
        """ Set up arguments and initialise DB client.

        Notes:
            Several arguments can be passed to this class from the command-line,
            and here are interpreted through *args:

        Parameters:
            db (str): the name of the collection to import to.
            scan (bool): whether or not to just scan the directory, rather
                than importing (automatically sets dryrun to true).
            dryrun (bool): perform whole process, up to actually importing to
                the database.
            tags (str): apply this tag to each structure added to database.
            force (bool): override rules about which folders can be imported
                into main database.
            recent_only (bool): if true, sort file lists by modification
                date and stop scanning when a file that already exists in
                database is found.

        """
        self.args = args[0]
        self.dryrun = self.args.get('dryrun')
        self.scan = self.args.get('scan')
        self.recent_only = self.args.get('recent_only')
        if self.scan:
            self.dryrun = True
        self.debug = self.args.get('debug')
        self.verbosity = self.args.get('verbosity') or 0
        self.config_fname = self.args.get('config')
        self.tags = self.args.get('tags')
        self.prototype = self.args.get('prototype')
        self.tag_dict = dict()
        self.tag_dict['tags'] = self.tags
        self.import_count = 0
        self.skipped = 0
        self.exclude_patterns = ['bad_castep', 'input']
        self.errors = 0
        self.struct_list = []
        self.path_list = []

        self.log = logging.getLogger('spatula')
        loglevel = {3: "DEBUG", 2: "INFO", 1: "WARN", 0: "ERROR"}
        self.log.setLevel(loglevel.get(self.verbosity, "WARN"))

        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.log.addHandler(handler)

        # I/O files
        if not self.dryrun:
            self.num_words = len(WORDS)
            self.num_nouns = len(NOUNS)

        if not self.scan:
            self.logfile = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
            self.manifest = tempfile.NamedTemporaryFile(mode='w+t', delete=False)

        if settings is None:
            self.settings = load_custom_settings(config_fname=self.config_fname,
                                                 debug=self.debug,
                                                 no_quickstart=self.args.get('no_quickstart'))
        else:
            self.settings = settings

        result = make_connection_to_collection(self.args.get('db'),
                                               check_collection=False,
                                               import_mode=True,
                                               override=self.args.get('no_quickstart'),
                                               mongo_settings=self.settings)

        self.client, self.db, self.collections = result
        # perform some relevant collection-dependent checks
        assert len(self.collections) == 1, 'Can only import to one collection.'
        self.repo = list(self.collections.values())[0]

        # if trying to import to the default repo, without doing a dryrun or forcing it, then
        # check if we're in the protected directory, i.e. the only one that is allowed to import
        # to the default collection
        default_collection = recursive_get(self.settings, ['mongo', 'default_collection'])
        try:
            default_collection_file_path = recursive_get(self.settings, ['mongo', 'default_collection_file_path'])
        except KeyError:
            default_collection_file_path = None

        if self.args.get('db') is None or self.args.get('db') == default_collection:
            if not self.dryrun and not self.args.get('force'):
                # if using default collection, check we are in the correct path
                if default_collection_file_path is not None:
                    if not os.getcwd().startswith(os.path.expanduser(default_collection_file_path)):
                        print(80 * '!')
                        print('You shouldn\'t be importing to the default database from this folder! '
                              'Please use --db <YourDBName> to create a new collection, '
                              'or copy these files to the correct place!')
                        print(80 * '!')
                        raise RuntimeError('Failed to import to default collection from '
                                           'current directory, import must be called from {}'
                                           .format(default_collection_file_path))
        else:
            if self.args.get('db') is not None:
                if any(['oqmd' in db for db in self.args.get('db')]):
                    exit('Cannot import directly to oqmd repo')
                elif len(self.args.get('db')) > 1:
                    exit('Can only import to one collection.')

        num_prototypes_in_db = self.repo.count_documents({'prototype': True})
        num_objects_in_db = self.repo.count_documents({})
        self.log.info("Found {} existing objects in database.".format(num_objects_in_db))
        if self.args.get('prototype'):
            if num_prototypes_in_db != num_objects_in_db:
                raise SystemExit('I will not import prototypes to a non-prototype database!')
        else:
            if num_prototypes_in_db != 0:
                raise SystemExit('I will not import DFT calculations into a prototype database!')

        if not self.dryrun:
            # either drop and recreate or create spatula report collection
            self.db.spatula.drop()
            self.report = self.db.spatula

        # scan directory on init
        self.file_lists = self._scan_dir()
        # if import, as opposed to rebuild, scan for duplicates and remove from list
        if not self.args.get('subcmd') == 'rebuild':
            self.file_lists, skipped = self._scan_dupes(self.file_lists)
            self.skipped += skipped

        # print number of files found
        self._display_import()

        # only create dicts if not just scanning
        if not self.scan:
            # convert to dict and db if required
            self._files2db(self.file_lists)
        if self.import_count == 0:
            print('No new structures imported!')
        if not self.dryrun and self.import_count > 0:
            print('Successfully imported', self.import_count, 'structures!')
            # index by enthalpy for faster/larger queries
            count = 0
            for _ in self.repo.list_indexes():
                count += 1
            # ignore default id index
            if count > 1:
                self.log.info('Index found, rebuilding...')
                self.repo.reindex()
            else:
                self.log.info('Building index...')
                self.repo.create_index([('enthalpy_per_atom', pm.ASCENDING)])
                self.repo.create_index([('stoichiometry', pm.ASCENDING)])
                self.repo.create_index([('cut_off_energy', pm.ASCENDING)])
                self.repo.create_index([('species_pot', pm.ASCENDING)])
                self.repo.create_index([('kpoints_mp_spacing', pm.ASCENDING)])
                self.repo.create_index([('xc_functional', pm.ASCENDING)])
                self.repo.create_index([('elems', pm.ASCENDING)])
                # index by source for rebuilds
                self.repo.create_index([('source', pm.ASCENDING)])
                self.log.info('Done!')
        elif self.dryrun:
            self.log.info('Dryrun complete!')

        if not self.scan:
            self.logfile.seek(0)
            errors = sum(1 for line in self.logfile)
            self.errors += errors
            if errors == 1:
                self.log.warning('There is 1 error to view in {}'.format(self.logfile.name))
            elif errors == 0:
                self.log.warning('There are no errors to view in {}'.format(self.logfile.name))
            elif errors > 1:
                self.log.warning('There are {} errors to view in {}'.format(errors, self.logfile.name))

        try:
            self.logfile.close()
        except Exception:
            pass
        try:
            self.manifest.close()
        except Exception:
            pass

        if not self.dryrun:
            # construct dictionary in spatula_report collection to hold info
            report_dict = dict()
            report_dict['last_modified'] = datetime.datetime.utcnow().replace(microsecond=0)
            report_dict['num_success'] = self.import_count
            report_dict['num_errors'] = errors
            report_dict['version'] = __version__
            self.report.insert_one(report_dict)

    def _struct2db(self, struct):
        """ Insert completed Python dictionary into chosen
        database, with generated text_id. Add quality factor
        for any missing data.

        Parameters:
            struct (dict): dictionary containing structure.

        Returns:
            int: 1 if successfully inputted into database, 0 otherwise.

        """
        try:
            plain_text_id = [WORDS[random.randint(0, self.num_words-1)].strip(),
                             NOUNS[random.randint(0, self.num_nouns-1)].strip()]
            struct['text_id'] = plain_text_id
            if 'tags' in self.tag_dict:
                struct['tags'] = self.tag_dict['tags']
            struct['quality'] = 5
            # if any missing info at all, score = 0
            # include elem set for faster querying
            if 'elems' not in struct:
                struct['elems'] = sorted(list(set(struct['atom_types'])))

            # check basic DFT params if we're not in a prototype DB
            if not self.args.get('prototype'):
                failed_checks = []
                if 'species_pot' not in struct:
                    struct['quality'] = 0
                    failed_checks.append('missing all pspots')
                else:
                    specified = []
                    for elem in struct['stoichiometry']:
                        # remove all points for a missing pseudo
                        if elem[0] not in struct['species_pot']:
                            struct['quality'] = 0
                            failed_checks.append('missing pspot for {}'.format(elem[0]))
                        else:
                            specified.append(elem[0])
                            # remove a point for a generic OTF pspot
                            if 'OTF' in struct['species_pot'][elem[0]].upper():
                                struct['quality'] -= 1
                                failed_checks.append('pspot not fully specified for {}'.format(elem[0]))
                    struct['species_pot'] = {species: struct['species_pot'][species] for species in specified}

                if 'xc_functional' not in struct:
                    struct['quality'] = 0
                    failed_checks.append('missing xc functional')

                if 'cut_off_energy' not in struct:
                    struct['quality'] -= 1
                if 'kpoints_mp_spacing' not in struct:
                    struct['quality'] -= 1

            else:
                struct['prototype'] = True
                struct['xc_functional'] = 'xxx'
                struct['enthalpy_per_atom'] = 0
                struct['enthalpy'] = 0
                struct['pressure'] = 0
                struct['species_pot'] = {}

            root_src = get_root_source(struct)
            struct['root_source'] = root_src
            exts = ['.castep', '.res', '.history', '.history.gz']
            for ext in exts:
                for src in struct['source']:
                    if src.endswith(ext):
                        expanded_root_src = src

            if struct['quality'] == 5:
                struct_id = self.repo.insert_one(struct).inserted_id
                self.struct_list.append((struct_id, root_src))
                self.manifest.write('+ {}\n'.format(expanded_root_src))
                if self.debug:
                    print('Inserted', struct_id)
            else:
                self.logfile.write('? {} failed quality checks: {}\n'.format(expanded_root_src, failed_checks))
                if self.debug:
                    print('Error with', root_src)
                return 0

        except Exception as exc:
            # this shouldn't fail, but if it does, fail loudly but cleanly
            self.logfile.write('! Importer produced an unexpected error {}. Final state of struct: {}\n'.format(exc, struct))
            tb.print_exc()
            return 0

        return 1

    def _files2db(self, file_lists):
        """ Take all files found by scan and appropriately create dicts
        holding all available data; optionally push to database. There are
        three possible routes this function will take:

            - AIRSS-style results, e.g. folder of res files, cell and param,
                will be scraped into a single dictionary per structure.
            - Straight CASTEP results will be scraped into one structure per
                .castep file.
            - Prototype-style results will take res files and turn them
                into simple structural data entries, without DFT parameters.

        Parameters:
            file_lists (dict): filenames and filetype counts stored by directory name key.

        """
        print('\n{:^52}'.format('###### RUNNING IMPORTER ######') + '\n')
        for _, root in enumerate(file_lists):
            root_str = root
            if root_str == '.':
                root_str = os.getcwd().split('/')[-1]
            self.log.info("Dictifying directory: {}".format(root_str))

            if self.args.get('prototype'):
                self.import_count += self._scrape_prototypes(file_lists, root)

            else:
                # default to only scraping castep files
                style = 'castep'
                # if there are multiple res files, per cell, assume we are working in "airss" mode
                if file_lists[root]['res_count'] > 0:
                    if (file_lists[root]['castep_count'] < file_lists[root]['res_count'] and
                            file_lists[root]['cell_count'] <= file_lists[root]['res_count']):
                        style = 'airss'
                if style == 'airss':
                    imported = self._scrape_multi_file_results(file_lists, root)
                # otherwise, we are just in a folder of CASTEP files
                elif style == 'castep':
                    imported = self._scrape_single_file_structures(file_lists, root)

            self.import_count += imported
            if imported > 0:
                print('Imported {} structures from {}'.format(imported, root))
                self.path_list.append(root)

        if self.struct_list and not self.dryrun:
            self._update_changelog()

    def _scrape_multi_file_results(self, file_lists, root):
        """ Add structures to database by parsing .res or .castep files., with DFT data
        scraped from .castep/.cell/.param files in the same folder, i.e. data from multiple
        files.

        Parameters:
            file_lists: dictionary containing counts of file types in each
                sub-directory.
            root: name of sub-directory.

        Returns:
            int: number of structures successfully imported.

        """
        multi = False  # are there multiple param/cell files?
        cell = False  # was the cell file successfully scraped?
        param = False  # was the param file successfully scraped?
        import_count = 0  # how many files have been successfully imported?
        success = False

        if file_lists[root]['param_count'] == 1:
            param_dict, success = param2dict(file_lists[root]['param'][0],
                                             debug=self.debug, noglob=True,
                                             verbosity=self.verbosity)
            param = success
            if not success:
                self.logfile.write(param_dict)
        elif file_lists[root]['param_count'] > 1:
            self.log.warning('Multiple param files found: {}'.format(file_lists[root]['param']))
            multi = True
        if file_lists[root]['cell_count'] == 1:
            cell_dict, success = cell2dict(file_lists[root]['cell'][0],
                                           db=True,
                                           debug=self.debug,
                                           noglob=True,
                                           verbosity=self.verbosity)
            cell = success
            if not success:
                self.logfile.write(str(cell_dict))
        elif file_lists[root]['cell_count'] > 1:
            multi = True
            self.log.warning('Multiple param files found: {}'.format(file_lists[root]['cell']))
        if multi:
            found_multi = False
            for param_name in file_lists[root]['param']:
                for cell_name in file_lists[root]['cell']:
                    if param_name.split('.')[0] in cell_name:
                        cell_dict, success = cell2dict(cell_name,
                                                       debug=self.debug,
                                                       db=True,
                                                       verbosity=self.verbosity)
                        cell = success
                        if not success:
                            self.logfile.write(str(cell_dict))
                            continue

                        param_dict, success = param2dict(param_name,
                                                         debug=self.debug,
                                                         verbosity=self.verbosity)
                        param = success
                        if not success:
                            self.logfile.write(param_dict)

                        if success:
                            found_multi = True
                            self.log.info('Found matching cell and param files: {}'.format(param_name))
                            break

            if not found_multi:
                self.log.warning("Unable to find matching cell and param files for {}".format(root))

        # combine cell and param dicts for folder
        input_dict = dict()
        if cell and param:
            input_dict.update(cell_dict)
            input_dict.update(param_dict)
            input_dict['source'] = cell_dict['source'] + param_dict['source']
        else:
            self.logfile.write('! {} failed to scrape any cell and param \n'.format(root))

        # create res dicts and combine them with input_dict
        for _, _file in enumerate(loading_bar(file_lists[root]['res'], verbosity=self.verbosity)):
            exts_with_precedence = ['.castep', '.history', 'history.gz']
            # check if a castep-like file exists instead of scraping res
            if any([_file.replace('.res', ext) in file_lists[root]['castep']
                    for ext in exts_with_precedence]):
                for ext in exts_with_precedence:
                    if _file.replace('.res', ext) in file_lists[root]['castep']:
                        struct_dict, success = castep2dict(_file.replace('.res', ext),
                                                           debug=False, noglob=True,
                                                           dryrun=self.args.get('dryrun'),
                                                           verbosity=self.verbosity)
                        break
            # otherwise, scrape res file
            else:
                struct_dict, success = res2dict(_file, verbosity=self.verbosity, noglob=True)

            if not success:
                self.logfile.write('! {}'.format(struct_dict))
            else:
                try:
                    final_struct = copy.deepcopy(input_dict)
                    final_struct.update(struct_dict)
                    # calculate kpoint spacing if not found
                    if 'lattice_cart' not in final_struct and 'lattice_abc' not in final_struct:
                        msg = '! {} missing lattice'.format(_file)
                        self.logfile.write(msg)

                    if 'kpoints_mp_spacing' not in final_struct and 'kpoints_mp_grid' in final_struct:
                        final_struct['kpoints_mp_spacing'] = calc_mp_spacing(final_struct['lattice_cart'],
                                                                             final_struct['mp_grid'])
                    final_struct['source'] = struct_dict['source']
                    if 'source' in input_dict:
                        final_struct['source'] += input_dict['source']

                    if not self.dryrun:
                        final_struct.update(self.tag_dict)
                        import_count += self._struct2db(final_struct)

                except Exception as exc:
                    self.log.error('Unexpected error for {}, {}'.format(_file, final_struct))
                    raise exc

        return import_count

    def _scrape_single_file_structures(self, file_lists, root):
        """ Scrape data from one of CASTEP-like, .synth and .expt files,
        and push to database.

        Parameters:
            file_lists (dict): filenames and filetype counts stored by directory name key.
            root (str): directory name to scrape.

        Returns:
            int: number of successfully structures imported to database.

        """
        import_count = 0
        for _, file in enumerate(loading_bar(file_lists[root]['castep'], verbosity=self.verbosity)):
            castep_dict, success = castep2dict(file, debug=False, verbosity=self.verbosity)
            if not success:
                self.logfile.write('! {}'.format(castep_dict))
            else:
                final_struct = castep_dict
                if not self.dryrun:
                    final_struct.update(self.tag_dict)
                    import_count += self._struct2db(final_struct)

        return import_count

    def _scrape_prototypes(self, file_lists, root):
        """ Scrape prototype data, i.e. structures with no DFT data, and push
        to database.

        Parameters:
            file_lists (dict): filenames and filetype counts stored by directory name key.
            root (str): directory name to scrape.

        Returns:
            int: number of successfully structures imported to database.

        """
        import_count = 0
        for _, file in enumerate(file_lists[root]['res']):
            res_dict, success = res2dict(file, db=False, verbosity=self.verbosity, noglob=True)
            if not success:
                self.logfile.write('! {}'.format(res_dict))
            else:
                final_struct = res_dict
                if not self.dryrun:
                    final_struct.update(self.tag_dict)
                    import_count += self._struct2db(final_struct)

        return import_count

    def _update_changelog(self):
        """ Add a list of ObjectIds to a collection called __changelog_{collection_name},
        storing "commits" that can be undone or reverted to.

        """
        changes = {'date': datetime.datetime.today(),
                   'count': len(self.struct_list),
                   'id_list': [struct[0] for struct in self.struct_list],
                   'src_list': [struct[1] for struct in self.struct_list],
                   'path_list': self.path_list}
        self.db['__changelog_{}'.format(self.repo.name)].insert_one(changes)

    def _scan_dir(self):
        """ Scans folder topdir recursively, returning list of
        CASTEP/AIRSS input/output files.

        Returns:
            dict: dictionary keyed by directory name that stores filename
                and filetype counts.

        """
        self.log.info("Scanning for directories...")
        file_lists = dict()
        topdir = '.'
        topdir_string = os.getcwd().split('/')[-1]
        self.log.info('Scanning {} for CASTEP/AIRSS output files... '.format(topdir_string))
        for root, _, files in os.walk(topdir, followlinks=True, topdown=True):
            for pattern in self.exclude_patterns:
                if pattern in root:
                    self.log.debug('Skipping directory {} as it matched exclude pattern {}'.format(root, pattern))
            # get absolute path for rebuilds
            root = os.path.abspath(root)
            file_lists[root] = collections.defaultdict(list)
            file_lists[root]['res_count'] = 0
            file_lists[root]['cell_count'] = 0
            file_lists[root]['param_count'] = 0
            file_lists[root]['castep_count'] = 0
            file_lists[root]['synth_count'] = 0
            file_lists[root]['expt_count'] = 0
            for _file in files:
                _file = root + '/' + _file
                self.log.debug(_file)
                if _file.endswith('.res'):
                    file_lists[root]['res'].append(_file)
                    file_lists[root]['res_count'] += 1
                elif (_file.endswith('.castep') or _file.endswith('.history') or
                      _file.endswith('.history.gz')):
                    file_lists[root]['castep'].append(_file)
                    file_lists[root]['castep_count'] += 1
                elif _file.endswith('.cell'):
                    if _file.endswith('-out.cell'):
                        continue
                    else:
                        file_lists[root]['cell'].append(_file)
                        file_lists[root]['cell_count'] += 1
                elif _file.endswith('.param'):
                    file_lists[root]['param'].append(_file)
                    file_lists[root]['param_count'] += 1
                elif _file.endswith('.synth'):
                    file_lists[root]['synth'].append(_file)
                    file_lists[root]['synth_count'] += 1
                elif _file.endswith('.expt'):
                    file_lists[root]['expt'].append(_file)
                    file_lists[root]['expt_count'] += 1

        # finally, filter the castep list so that history/castep/history.gz duplicates
        # are not included
        for root in file_lists:
            cached_list = [entry.split('/')[-1].split('.')[0] for entry in file_lists[root]['castep']]
            duplicates = []
            for ind, file in enumerate(file_lists[root]['castep']):
                if not file.endswith('.castep') and cached_list.count(file.split('/')[-1].split('.')[0]) > 1:
                    duplicates.append(ind)
                    file_lists[root]['castep_count'] -= 1
            file_lists[root]['castep'] = [entry for ind, entry in
                                          enumerate(file_lists[root]['castep']) if ind not in duplicates]

        self.log.info("Scanning completed.")
        return file_lists

    def _scan_dupes(self, file_lists):
        """ Scan the file_lists made by scan_dir and remove
        structures already in the database by matching sources.

        Parameters:
            file_lists (dict): dict with directory name keys containing
                file names and filetype counts.

        Returns:
            dict: the input dict minus duplicates.

        """
        self.log.info("Scanning for duplicates...")
        skipped = 0
        new_file_lists = copy.deepcopy(file_lists)
        for _, root in enumerate(file_lists):
            # per folder delete list
            # do not import castep or res if seed name in db already
            delete_list = {}
            delete_list['castep'] = set()
            delete_list['res'] = set()
            types = ['castep', 'res']
            for structure_type in types:
                if self.recent_only:
                    try:
                        new_file_lists[root][structure_type].sort(key=lambda x: os.stat(x).st_ctime)
                    except Exception:
                        self.log.warning('Unable to sort files by creation date, ignoring `recent_only=True`')

                for file_ind, _file in enumerate(new_file_lists[root][structure_type]):
                    # find number of entries with same root filename in database
                    structure_exts = ['.castep', '.res', '.history', '.history.gz']
                    ext = [_ext for _ext in structure_exts if _file.endswith(_ext)]
                    assert len(ext) == 1
                    ext = ext[0]
                    structure_count = 0
                    for other_ext in structure_exts:
                        structure_count += self.repo.count_documents(
                            {'source': {'$in': [_file.replace(ext, other_ext)]}})

                    if structure_count > 1:
                        self.log.error('Found double in database, this needs to be dealt with manually: {}'.format(_file))

                    # if duplicate found, don't reimport
                    if structure_count >= 1:
                        # need to loop over possible file types to get correct fname
                        _add_to_delete_lists(_file, root, new_file_lists, delete_list)
                        if self.recent_only:
                            skipping = len(new_file_lists[root][structure_type][file_ind:])
                            skipped += skipping
                            for other_file in new_file_lists[root][structure_type][file_ind:]:
                                _add_to_delete_lists(other_file, root, new_file_lists, delete_list)
                            break

                        assert len(new_file_lists[root][structure_type]) == len(file_lists[root][structure_type]) + skipped

            for _type in types:
                new_file_lists[root][_type] = [_file for _file in new_file_lists[root][_type]
                                               if _file not in delete_list[_type]]

            for file_type in types:
                if len(new_file_lists[root][file_type]) != new_file_lists[root][file_type + '_count']:
                    raise RuntimeError('{} count {} does not match number of files {}'
                                       .format(file_type,
                                               len(new_file_lists[root][file_type]),
                                               new_file_lists[root][file_type + '_count']))

        self.log.info("Removed duplicates from import list...")
        return new_file_lists, skipped

    def _display_import(self):
        """ Display number of files to be imported and a breakdown of
        their types.

        """
        exts = ['res', 'cell', 'castep', 'param', 'synth', 'expt']
        counts = {ext: 0 for ext in exts}
        for root in self.file_lists:
            for ext in exts:
                counts[ext] += self.file_lists[root]['{}_count'.format(ext)]

        print('\n\n')
        if self.recent_only:
            print('\t\tSkipped {} files filtered by creation date (st_ctime).'.format(self.skipped))

        for ext in exts:
            print('\t\t{:8d}\t\t.{} files'.format(counts[ext], ext))


def _add_to_delete_lists(_file, root, new_file_lists, delete_list):
    """ Add files to the delete list, with the correct file types. """
    structure_exts = {'castep': 'castep', 'res': 'res', 'history': 'castep', 'history.gz': 'castep'}
    for _type in structure_exts:
        fname_trial = _file.replace(_file.split('.')[-1], _type)
        list_type = structure_exts[_type]
        if fname_trial in new_file_lists[root][list_type] and fname_trial not in delete_list[list_type]:
            delete_list[list_type].add(fname_trial)
            new_file_lists[root][list_type + '_count'] -= 1
