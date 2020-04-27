""" Run some importer tests directly, rather than from CLI. """

import unittest
import os
import glob
import mongomock

from matador.db.importer import Spatula
from matador.query import DBQuery

REAL_PATH = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"
CONFIG_FNAME = None
DB_NAME = "ci_test"
ROOT_DIR = os.getcwd()
DEBUG = False
MONGO_CLIENT = mongomock.MongoClient()


@mongomock.patch(servers=(("mongo_test.com:27017"),))
class TestDatabaseImport(unittest.TestCase):
    """ Tests the Spatula class. """

    def setUp(self):
        from matador.config import load_custom_settings, SETTINGS

        SETTINGS.reset()
        SETTINGS = load_custom_settings(config_fname=CONFIG_FNAME, debug=True)
        SETTINGS["mongo"]["default_collection"] = DB_NAME
        SETTINGS["mongo"]["default_collection_file_path"] = "/data/"
        SETTINGS["mongo"]["host"] = "mongo_test.com"
        SETTINGS["mongo"]["port"] = 27017
        self.settings = SETTINGS

    def tearDown(self):
        from matador.config import SETTINGS

        SETTINGS.reset()

    def test_failed_import(self):
        """ Try to import to the default collection from some
        random folder.

        """
        os.chdir(REAL_PATH + "/data/castep_files")

        args = {"no_quickstart": True}
        with self.assertRaises(RuntimeError):
            importer = Spatula(args, settings=self.settings)
        for _file in glob.glob("spatula.*"):
            os.remove(_file)

        args = {"db": DB_NAME, "no_quickstart": True}
        with self.assertRaises(RuntimeError):
            importer = Spatula(args, settings=self.settings)
        for _file in glob.glob("spatula.*"):
            os.remove(_file)

        args = {"dryrun": True, "no_quickstart": True}
        importer = Spatula(args, settings=self.settings)
        print(importer.errors)
        self.assertEqual(importer.import_count, 0)
        self.assertEqual(importer.skipped, 0)
        self.assertEqual(importer.errors, 3)
        for _file in glob.glob("spatula.*"):
            os.remove(_file)

        args = {"force": True, "no_quickstart": True}
        importer = Spatula(args, settings=self.settings)
        self.assertEqual(importer.import_count, 3)
        self.assertEqual(importer.skipped, 0)
        self.assertEqual(importer.errors, 5)
        for _file in glob.glob("spatula.*"):
            os.remove(_file)

        args = {"force": True, "recent_only": True, "no_quickstart": True}
        importer = Spatula(args, settings=self.settings)
        # this can be fiddly and a bit system dependent, just check errors+skips are constant
        self.assertEqual(importer.skipped + importer.errors, 8)
        self.assertEqual(importer.import_count, 0)
        for _file in glob.glob("spatula.*"):
            os.remove(_file)

        query = DBQuery(db=DB_NAME, mongo_settings=self.settings)
        self.assertEqual(len(query.cursor), 3)

        query = DBQuery(
            db=DB_NAME, mongo_settings=self.settings, available_values="root_source"
        )
        self.assertEqual(len(query.cursor), 3)

        query = DBQuery(
            db=DB_NAME,
            mongo_settings=self.settings,
            root_src="Na3Zn4-swap-ReOs-OQMD_759599",
        )
        self.assertEqual(len(query.cursor), 1)

        query = DBQuery(db=DB_NAME, mongo_settings=self.settings, id="no chance")
        self.assertEqual(len(query.cursor), 0)
