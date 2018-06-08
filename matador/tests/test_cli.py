#!/usr/bin/env python

""" This file implements several tests that directly act on
a database instance running locally, more integration tests
than unittests per se.

"""


import unittest
import os
import sys
import glob
import pymongo as pm

import matador.cli.cli
from matador.config import load_custom_settings
from matador.query import DBQuery
from matador.scrapers.castep_scrapers import cell2dict, res2dict

REAL_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
if os.uname()[1] == 'cluster2':
    CONFIG_FNAME = None
else:
    CONFIG_FNAME = REAL_PATH + 'data/matador_pipelines_conf.yml'
DB_NAME = 'ci_test'
ROOT_DIR = os.getcwd()
SETTINGS = load_custom_settings(config_fname=CONFIG_FNAME)

MONGO_PRESENT = True
try:
    MONGO_CLIENT = pm.MongoClient(SETTINGS['mongo']['host'], serverSelectionTimeoutMS=1000)
    MONGO_DB_NAMES = MONGO_CLIENT.database_names()
except pm.errors.ServerSelectionTimeoutError:
    MONGO_PRESENT = False

MONGO_CLIENT.crystals[DB_NAME].drop()
MONGO_CLIENT.crystals['__changelog_' + DB_NAME].drop()

for _file in glob.glob('*spatula*'):
    os.remove(_file)


@unittest.skipIf(not MONGO_PRESENT, 'MongoDB instance not found, skipping tests...')
class IntegrationTest(unittest.TestCase):
    """ Test  functionality acting on local database. """
    def testIntegration(self):
        """ Test import and query. """
        query, files_to_delete, err_flines, manifest_flines = test_import_castep()
        self.assertEqual(len(files_to_delete), 2, msg='Failed to write spatula files')
        print(err_flines)
        print(manifest_flines)
        self.assertEqual(len(err_flines), 3, msg='Failed to report errors correctly')
        self.assertEqual(len(manifest_flines), 3, msg='Failed to report successes correctly')
        self.assertEqual(len(query.cursor), 3, msg='Failed to import structures correctly')

        query_1, query_2, files_to_delete = test_import_res()
        self.assertEqual(len(query_1.cursor), 7, msg='Failed to import res files')
        self.assertEqual(len(query_2.cursor), 4, msg='Failed to import res files')
        self.assertEqual(len(files_to_delete), 2, msg='Failed to write spatula files')

        output_folder_exists, successes, elem_successes = test_swaps()
        self.assertTrue(output_folder_exists, msg='No folder created')
        self.assertTrue(all(successes), msg='Failed to even read files')
        self.assertFalse(all(elem_successes), msg='Swaps had wrong elements')

        query_1, query_2, changes_count = test_changes()
        self.assertEqual(len(query_1.cursor), 3, msg='matador changes did not remove files')
        self.assertEqual(len(query_2.cursor), 0, msg='matador changes did not remove files')
        self.assertEqual(changes_count, 1, msg='matador changes did not changelog')


def test_import_castep():
    """ Import from castep files, returning data to be checked. """
    # import from CASTEP files only
    os.chdir(REAL_PATH + '/data/castep_files')
    sys.argv = ['matador', 'import', '--db', DB_NAME]

    if CONFIG_FNAME is not None:
        sys.argv += ['--config', CONFIG_FNAME]

    matador.cli.cli.main(override=True)

    query = DBQuery(db=DB_NAME, config=CONFIG_FNAME)

    with open('spatula.err', 'r') as err_file:
        err_flines = err_file.readlines()
    with open('spatula.manifest', 'r') as manifest_file:
        manifest_flines = manifest_file.readlines()

    files_to_delete = glob.glob('*spatula*')
    for f in files_to_delete:
        os.remove(f)
    os.chdir(REAL_PATH)

    return query, files_to_delete, err_flines, manifest_flines


def test_import_res():
    """ Import from res files, returning data to be checked. """
    # import from combined res/cell/param files
    os.chdir(REAL_PATH + '/data/res_files')

    sys.argv = ['matador', 'import', '--db', DB_NAME]

    if CONFIG_FNAME is not None:
        sys.argv += ['--config', CONFIG_FNAME]

    matador.cli.cli.main(override=True)

    query_1 = DBQuery(db=DB_NAME, config=CONFIG_FNAME)
    query_2 = DBQuery(db=DB_NAME, composition='KSnP', config=CONFIG_FNAME)

    files_to_delete = glob.glob('*spatula*')
    for f in files_to_delete:
        os.remove(f)

    os.chdir(REAL_PATH)

    return query_1, query_2, files_to_delete


def test_swaps():
    """ Run swaps, returning data to be checked. """
    sys.argv = ['matador', 'swaps', '--db', DB_NAME, '--res', '--cell',
                '-c', 'NaPZn', '-int', '-sw', 'NaLi:PSi:ZnSn']
    if CONFIG_FNAME is not None:
        sys.argv += ['--config', CONFIG_FNAME]
    matador.cli.cli.main()
    expected_dir = 'swaps-NaPZn-ci_test-NaLi:PSi:ZnSn'
    output_folder_exists = os.path.isdir(expected_dir)
    successes = []
    elem_successes = []
    if output_folder_exists:
        os.chdir(expected_dir)
        expected_files = ['LiSi-swap-NaP_intermediates', 'LiSn-swap-Na3Zn4-OQMD_759599', 'Li-swap-Na-edgecase-CollCode10101']
        target_elems = ['Li', 'Si', 'Sn']

        for files in expected_files:
            cell, s = cell2dict(files, db=False, positions=True)
            successes.append(s)
            if s:
                elems = set(cell['atom_types'])
                elem_successes.append(any([elem not in target_elems for elem in elems]))
            else:
                elem_successes.append(False)
            res, s = res2dict(files, db=False)
            successes.append(s)
            if s:
                elems = set(res['atom_types'])
                elem_successes.append(any([elem not in target_elems for elem in elems]))
            else:
                elem_successes.append(False)

            os.remove(files + '.cell')
            os.remove(files + '.res')

        os.chdir(REAL_PATH)
        os.rmdir(expected_dir)

    return output_folder_exists, successes, elem_successes


def test_changes():
    """ Test matador changes functionality by undoing the second change (i.e. the res import). """
    sys.argv = ['matador', 'changes', '--db', DB_NAME, '-c', '2', '--undo']

    if CONFIG_FNAME is not None:
        sys.argv += ['--config', CONFIG_FNAME]

    matador.cli.cli.main(override=True)

    query_1 = DBQuery(db=DB_NAME, config=CONFIG_FNAME)
    query_2 = DBQuery(db=DB_NAME, composition='KSnP', config=CONFIG_FNAME)
    changes_count = MONGO_CLIENT.crystals['__changelog_' + DB_NAME].count()

    return query_1, query_2, changes_count


if __name__ == '__main__':
    unittest.main()
