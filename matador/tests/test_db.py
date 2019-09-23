""" Run some importer tests directly, rather than from CLI. """

import unittest
import os
import glob
import getpass
import pymongo as pm

from matador.config import load_custom_settings
from matador.db.importer import Spatula

REAL_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
CONFIG_FNAME = None
DB_NAME = 'ci_test'
ROOT_DIR = os.getcwd()
SETTINGS = load_custom_settings(config_fname=CONFIG_FNAME, override=True)
SETTINGS['mongo']['default_collection'] = DB_NAME
SETTINGS['mongo']['default_collection_file_path'] = '/data/'

DEBUG = False
MONGO_PRESENT = True
try:
    MONGO_CLIENT = pm.MongoClient(SETTINGS['mongo']['host'], serverSelectionTimeoutMS=1000)
    MONGO_DB_NAMES = MONGO_CLIENT.list_database_names()
    MONGO_CLIENT.crystals[DB_NAME].drop()
    MONGO_CLIENT.crystals[getpass.getuser() + '_' + DB_NAME].drop()
except pm.errors.ServerSelectionTimeoutError:
    MONGO_PRESENT = False


@unittest.skipIf(not MONGO_PRESENT, 'MongoDB not found, skipping')
class TestDatabase(unittest.TestCase):
    """ Tests the Spatula class. """
    def test_failed_import(self):
        """ Try to import to the default collection from some
        random folder.

        """
        os.chdir(REAL_PATH + '/data/castep_files')

        args = {'override': True}
        with self.assertRaises(RuntimeError):
            importer = Spatula(args, settings=SETTINGS)
        for _file in glob.glob('spatula.*'):
            os.remove(_file)

        args = {'db': DB_NAME, 'override': True}
        with self.assertRaises(RuntimeError):
            importer = Spatula(args, settings=SETTINGS)
        for _file in glob.glob('spatula.*'):
            os.remove(_file)

        args = {'dryrun': True, 'override': True}
        importer = Spatula(args, settings=SETTINGS)
        self.assertEqual(importer.import_count, 0)
        self.assertEqual(importer.skipped, 0)
        self.assertEqual(importer.errors, 3)
        for _file in glob.glob('spatula.*'):
            os.remove(_file)

        args = {'force': True, 'override': True}
        importer = Spatula(args, settings=SETTINGS)
        self.assertEqual(importer.import_count, 3)
        self.assertEqual(importer.skipped, 0)
        self.assertEqual(importer.errors, 5)
        for _file in glob.glob('spatula.*'):
            os.remove(_file)

        args = {'force': True, 'recent_only': True, 'override': True}
        importer = Spatula(args, settings=SETTINGS)
        self.assertEqual(importer.skipped, 5)
        self.assertEqual(importer.errors, 3)
        self.assertEqual(importer.import_count, 0)
        for _file in glob.glob('spatula.*'):
            os.remove(_file)
