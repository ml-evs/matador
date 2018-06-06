#!/usr/bin/env python
import unittest
import os
import sys
import glob
import pymongo as pm

import matador.cli.cli
from matador.config import load_custom_settings

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
    cli = pm.MongoClient(SETTINGS['mongo']['host'], serverSelectionTimeoutMS=1000)
    dbs = cli.database_names()
except pm.errors.ServerSelectionTimeoutError:
    MONGO_PRESENT = False

try:
    cli.crystals[DB_NAME].drop()
except:
    pass

for _file in glob.glob('*spatula*'):
    os.remove(_file)


@unittest.skipIf(not MONGO_PRESENT, 'MongoDB instance not found, skipping tests...')
class CLIIntegrationTest(unittest.TestCase):
    """ Test CLI functionality. """
    def testImporterAndMore(self):
        """ Test import and query. """
        from matador.query import DBQuery
        from matador.scrapers.castep_scrapers import cell2dict, res2dict

        # import from CASTEP files only
        os.chdir(REAL_PATH + '/data/castep_files')
        sys.argv = ['matador', 'import', '--db', DB_NAME]

        if CONFIG_FNAME is not None:
            sys.argv += ['--config', CONFIG_FNAME]

        matador.cli.cli.main(override=True)

        query = DBQuery(db='ci_test', config=CONFIG_FNAME)
        self.assertEqual(len(query.cursor), 3)

        with open('spatula.err', 'r') as err_file:
            err_flines = err_file.readlines()
        with open('spatula.manifest', 'r') as manifest_file:
            manifest_flines = manifest_file.readlines()

        files_to_delete = glob.glob('*spatula*')
        self.assertEqual(len(files_to_delete), 2)
        self.assertEqual(len(err_flines), 3)
        self.assertEqual(len(manifest_flines), 3)
        for f in files_to_delete:
            os.remove(f)

        # import from combined res/cell/param files
        os.chdir(REAL_PATH + '/data/res_files')

        sys.argv = ['matador', 'import', '--db', DB_NAME]

        if CONFIG_FNAME is not None:
            sys.argv += ['--config', CONFIG_FNAME]

        matador.cli.cli.main(override=True)

        query = DBQuery(db='ci_test', config=CONFIG_FNAME)
        self.assertEqual(len(query.cursor), 7)

        query = DBQuery(db='ci_test', composition='KSnP', config=CONFIG_FNAME)
        self.assertEqual(len(query.cursor), 4)

        files_to_delete = glob.glob('*spatula*')
        self.assertEqual(len(files_to_delete), 2)
        for f in files_to_delete:
            os.remove(f)

        os.chdir(REAL_PATH)

        sys.argv = ['matador', 'swaps', '--db', DB_NAME, '--res', '--cell',
                    '-c', 'NaPZn', '-int', '-sw', 'NaLi:PSi:ZnSn']

        if CONFIG_FNAME is not None:
            sys.argv += ['--config', CONFIG_FNAME]

        matador.cli.cli.main()

        expected_dir = 'swaps-NaPZn-ci_test-NaLi:PSi:ZnSn'
        output_folder_exists = os.path.isdir(expected_dir)
        print(output_folder_exists)
        successes = []
        elem_successes = []

        if output_folder_exists:
            os.chdir(expected_dir)
            expected_files = ['LiSi-swap-NaP_intermediates', 'LiSn-swap-Na3Zn4-OQMD_759599', 'Li-swap-Na-edgecase']
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

        self.assertTrue(output_folder_exists, msg='No folder created')
        self.assertTrue(all(successes), msg='Failed to even read files')
        self.assertFalse(all(elem_successes), msg='Swaps had wrong elements')


if __name__ == '__main__':
    unittest.main()
