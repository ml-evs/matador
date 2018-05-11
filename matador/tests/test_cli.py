#!/usr/bin/env python
import unittest
import subprocess as sp
import os
import sys

import matador.cli.matador_cli

REAL_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
CONFIG_FNAME = REAL_PATH + 'data/matador_pipelines_conf.json'
DB_NAME = 'ci_test'
ROOT_DIR = os.getcwd()

MONGO_PRESENT = True
import pymongo as pm
try:
    cli = pm.MongoClient(serverSelectionTimeoutMS=100)
    dbs = cli.database_names()
except pm.errors.ServerSelectionTimeoutError:
    MONGO_PRESENT = False


@unittest.skipIf(not MONGO_PRESENT , 'MongoDB instance not found, skipping tests...')
class CLITest(unittest.TestCase):
    """ Test CLI functionality. """
    def testImporter(self):
        from matador.query import DBQuery
        os.chdir(REAL_PATH + '/data/castep_files')
        sys.argv = ['/home/matthew/.local/conda/bin/matador', 'import', '--db', DB_NAME, '--config', CONFIG_FNAME]

        matador.cli.matador_cli.main(testing=True)

        query = DBQuery(db='ci_test', config=CONFIG_FNAME)
        self.assertEqual(len(query.cursor), 4)


if __name__ == '__main__':
    unittest.main()
