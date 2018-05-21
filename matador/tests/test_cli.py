#!/usr/bin/env python
import unittest
import subprocess as sp
import os
import sys
import glob
import pymongo as pm

import matador.cli.cli
from matador.utils.db_utils import load_custom_settings

REAL_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
if os.uname()[1] == 'cluster2':
    CONFIG_FNAME = REAL_PATH + '../config/matador_conf.json'
else:
    CONFIG_FNAME = REAL_PATH + 'data/matador_pipelines_conf.json'
DB_NAME = 'ci_test'
ROOT_DIR = os.getcwd()
SETTINGS = load_custom_settings(CONFIG_FNAME)

MONGO_PRESENT = True
print(SETTINGS)
try:
    cli = pm.MongoClient(SETTINGS['mongo']['host'], serverSelectionTimeoutMS=1000)
    dbs = cli.database_names()
except pm.errors.ServerSelectionTimeoutError:
    MONGO_PRESENT = False


@unittest.skipIf(not MONGO_PRESENT, 'MongoDB instance not found, skipping tests...')
class CLITest(unittest.TestCase):
    """ Test CLI functionality. """
    def testImporter(self):
        """ Test import and query. """
        from matador.query import DBQuery
        os.chdir(REAL_PATH + '/data/castep_files')
        sys.argv = ['matador', 'import', '--db', DB_NAME, '--config', CONFIG_FNAME]

        matador.cli.cli.main(testing=True)

        query = DBQuery(db='ci_test', config=CONFIG_FNAME)
        self.assertEqual(len(query.cursor), 4)

        files_to_delete = glob.glob('*spatula*')
        self.assertEqual(len(files_to_delete), 2)
        for f in files_to_delete:
            os.remove(f)


if __name__ == '__main__':
    unittest.main()
