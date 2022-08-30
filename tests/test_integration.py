#!/usr/bin/env python

""" This file implements several tests that directly act on
a database instance running locally, more integration tests
than unittests per se.

"""


import unittest
import os
import sys
import glob

import mongomock

import matador.cli.cli
from matador.query import DBQuery
from matador.hull import QueryConvexHull
from matador.scrapers.castep_scrapers import cell2dict, res2dict
from matador.utils.chem_utils import get_root_source

REAL_PATH = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"
CONFIG_FNAME = None
DB_NAME = "ci_test"
ROOT_DIR = os.getcwd()
DEBUG = False
MONGO_CLIENT = mongomock.MongoClient()
OUTPUT_DIR = ROOT_DIR + "/integration_test_results"


@mongomock.patch(servers=("mongo_test.com:27017",))
class IntegrationTest(unittest.TestCase):
    """Test functionality acting on local database."""

    def tearDown(self):
        import shutil

        self.settings.reset()
        os.chdir(ROOT_DIR)
        shutil.rmtree(OUTPUT_DIR)

    def setUp(self):
        from matador.config import load_custom_settings, SETTINGS

        SETTINGS.reset()
        _ = load_custom_settings(config_fname=CONFIG_FNAME)
        SETTINGS["mongo"]["default_collection"] = DB_NAME
        SETTINGS["mongo"]["default_collection_file_path"] = "/data/"
        SETTINGS["mongo"]["host"] = "mongo_test.com"
        SETTINGS["mongo"]["port"] = 27017
        self.settings = SETTINGS

        os.makedirs(OUTPUT_DIR, exist_ok=False)
        os.chdir(OUTPUT_DIR)

    def test_integration(self):
        """Test import and query."""
        print("IMPORT CASTEP 1")
        query = import_castep()
        self.assertEqual(
            len(query.cursor), 3, msg="Failed to import structures correctly"
        )

        # run again and hopefully nothing will change, i.e. no duplication
        print("IMPORT CASTEP 2")
        query = import_castep()
        self.assertEqual(
            len(query.cursor), 3, msg="Failed to import structures correctly"
        )

        print("IMPORT CASTEP 3")
        query = import_castep(extra_flags="--recent_only")
        self.assertEqual(
            len(query.cursor), 3, msg="Failed to import structures correctly"
        )

        print("IMPORT RES")
        query_1, query_2 = import_res()
        self.assertEqual(len(query_1.cursor), 7, msg="Failed to import res files")
        self.assertEqual(len(query_2.cursor), 4, msg="Failed to import res files")
        self.assertEqual(
            query_2.cursor[0]["species_pot"]["K"],
            "2|1.5|9|10|11|30U:40:31(qc=6)",
            msg="Failed to scrape OTF with weird name",
        )
        self.assertEqual(
            query_2.cursor[0]["species_pot"]["Sn"],
            "2|2|2|1.6|9.6|10.8|11.7|50U=-0.395U=+0.25:51U=-0.14U=+0.25",
            msg="Failed to scrape OTF with linebreak",
        )
        self.assertFalse(
            any(["Sb" in doc["species_pot"] for doc in query_2.cursor]),
            msg="pspots over-scraped!",
        )

        print("SWAP")
        output_folder_exists, successes, elem_successes = swaps()
        self.assertTrue(output_folder_exists, msg="No folder created")
        self.assertTrue(all(successes), msg="Failed to even read files")
        self.assertFalse(all(elem_successes), msg="Swaps had wrong elements")

        print("CHANGES")
        query_1, query_2, changes_count = changes()
        self.assertEqual(
            len(query_1.cursor), 3, msg="matador changes did not remove files"
        )
        self.assertEqual(
            len(query_2.cursor), 0, msg="matador changes did not remove files"
        )
        # unclear that mongomock can handle this; just check the queries instead
        # self.assertEqual(changes_count, 1, msg='matador changes did not changelog')

        expected_dir = "query-ci_test"
        if os.path.isdir(expected_dir):
            for _file in glob.glob(expected_dir + "/*"):
                os.remove(_file)
            os.removedirs(expected_dir)

        print("EXPORT")
        export()
        expected_files = [
            "query-ci_test/query-ci_test.md",
            "query-ci_test/query-ci_test.tex",
            "query-ci_test/Na3Zn4-swap-ReOs-OQMD_759599.pdb",
            "query-ci_test/Na3Zn4-swap-ReOs-OQMD_759599.xsf",
            "query-ci_test/Na3Zn4-swap-ReOs-OQMD_759599.json",
            "query-ci_test/Na-edgecase-CollCode10101.pdb",
            "query-ci_test/Na-edgecase-CollCode10101.xsf",
            "query-ci_test/Na-edgecase-CollCode10101.json",
            "query-ci_test/NaP_intermediates.pdb",
            "query-ci_test/NaP_intermediates.xsf",
            "query-ci_test/NaP_intermediates.json",
        ]
        dir_exists = os.path.isdir(expected_dir)
        files_exist = all([os.path.isfile(_file) for _file in expected_files])

        for _file in expected_files:
            if os.path.isfile(_file):
                os.remove(_file)

        if os.path.isdir(expected_dir):
            for _file in glob.glob(expected_dir + "/*"):
                os.remove(_file)
            os.removedirs(expected_dir)

        self.assertTrue(dir_exists, msg="Failed to create output directory")
        self.assertTrue(files_exist, msg="Some files missing from export")

        print("PSEUDOTERNARY HULL")
        query, hull = pseudoternary_hull()
        self.assertTrue(query.args.get("intersection"))
        self.assertTrue(query._non_elemental)
        self.assertTrue(query._create_hull)
        self.assertEqual(len(query.cursor), 7)
        self.assertEqual(len(hull.cursor), 7)
        self.assertEqual(len(hull.hull_cursor), 5)

        expected_dir = "query-LaLiZrO-ci_test"
        if os.path.isdir(expected_dir):
            for _file in glob.glob(expected_dir + "/*"):
                os.remove(_file)
            os.removedirs(expected_dir)

        print("PSEUDOTERNARY HULL 2")
        hull = pseudoternary_hull_no_query()
        self.assertTrue(hull._query.args.get("intersection"))
        self.assertTrue(hull._query._non_elemental)
        self.assertTrue(hull._query._create_hull)
        self.assertEqual(len(query.cursor), 7)
        self.assertEqual(len(hull.cursor), 7)
        self.assertEqual(len(hull.hull_cursor), 5)

        expected_dir = "query-LaLiZrO-ci_test"
        if os.path.isdir(expected_dir):
            for _file in glob.glob(expected_dir + "/*"):
                os.remove(_file)
            os.removedirs(expected_dir)

        print("UNIQ")
        uniq()
        expected_files = [expected_dir + "/cubic-LLZO-CollCode999999.res"]

        dir_exists = os.path.isdir(expected_dir)
        files_exist = all([os.path.isfile(_file) for _file in expected_files])
        correct_num = len(glob.glob(expected_dir + "/*.res"))

        if os.path.isdir(expected_dir):
            for _file in glob.glob(expected_dir + "/*"):
                os.remove(_file)
            os.removedirs(expected_dir)

        self.assertTrue(dir_exists, msg="Failed to create output directory, uniq")
        self.assertTrue(files_exist, msg="Some files missing from export, uniq")
        self.assertTrue(correct_num == len(expected_files), msg="Incorrect filter")

        print("STATS")
        try:
            stats()
        except TypeError:
            print("Unable to test stats module due to mongomock limitations.")
            pass

        print("ID QUERY")
        query = id_query()
        self.assertEqual(len(query.cursor), 0)

        print("REFINE")
        cursor = refine().cursor
        self.assertTrue(all([doc["doi"] == ["10/12345"] for doc in cursor]))
        self.assertTrue(all([doc["tags"] == ["integration_test"] for doc in cursor]))
        self.assertTrue(
            all([doc["root_source"] == get_root_source(doc) for doc in cursor])
        )
        self.assertTrue(all([isinstance(doc["_raw"], list) for doc in cursor]))


def import_castep(extra_flags=None):
    """Import from castep files, returning data to be checked."""
    # import from CASTEP files only
    os.chdir(REAL_PATH + "/data/castep_files")
    sys.argv = ["matador", "import", "--force", "--db", DB_NAME]

    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]

    if DEBUG:
        sys.argv += ["--debug"]

    matador.cli.cli.main(no_quickstart=True)

    query = DBQuery(db=DB_NAME, config=CONFIG_FNAME, details=True, source=True)

    files_to_delete = glob.glob("*spatula*")
    for f in files_to_delete:
        os.remove(f)
    os.chdir(OUTPUT_DIR)

    return query


def import_res():
    """Import from res files, returning data to be checked."""
    # import from combined res/cell/param files
    os.chdir(REAL_PATH + "/data/res_files")

    sys.argv = ["matador", "import", "--force", "--db", DB_NAME]

    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]

    if DEBUG:
        sys.argv += ["--debug"]

    matador.cli.cli.main(no_quickstart=True)

    query_1 = DBQuery(db=DB_NAME, config=CONFIG_FNAME, details=True, source=True)
    query_2 = DBQuery(
        db=DB_NAME, composition="KSnP", config=CONFIG_FNAME, details=True, source=True
    )

    files_to_delete = glob.glob("*spatula*")
    for f in files_to_delete:
        os.remove(f)

    os.chdir(OUTPUT_DIR)

    return query_1, query_2


def pseudoternary_hull():
    """Import some other res files ready to make a hull."""
    os.chdir(REAL_PATH + "/data/hull-LLZO")
    sys.argv = ["matador", "import", "--force", "--db", DB_NAME]

    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]

    if DEBUG:
        sys.argv += ["--debug"]

    matador.cli.cli.main(no_quickstart=True)

    query = DBQuery(
        db=DB_NAME,
        subcmd="hull",
        composition="La2O3:Li2O:ZrO2",
        config=CONFIG_FNAME,
        details=True,
        source=True,
        no_plot=True,
    )
    hull = QueryConvexHull(query=query)

    os.chdir(OUTPUT_DIR)

    return query, hull


def pseudoternary_hull_no_query():
    """Import some other res files ready to make a hull."""
    hull = QueryConvexHull(
        db=DB_NAME,
        composition="La2O3:Li2O:ZrO2",
        config=CONFIG_FNAME,
        details=True,
        source=True,
        no_plot=True,
    )

    os.chdir(OUTPUT_DIR)

    return hull


def stats():
    """Run the stats command just to check it works at all."""
    sys.argv = ["matador", "stats", "--db", DB_NAME]
    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]
    matador.cli.cli.main()


def swaps():
    """Run swaps, returning data to be checked."""
    sys.argv = [
        "matador",
        "swaps",
        "--db",
        DB_NAME,
        "--res",
        "--cell",
        "-c",
        "NaPZn",
        "-int",
        "-sw",
        "NaLi:PSi:ZnSn",
    ]
    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]
    matador.cli.cli.main()
    expected_dir = "swaps-NaPZn-ci_test-NaLi:PSi:ZnSn"
    output_folder_exists = os.path.isdir(expected_dir)
    successes = []
    elem_successes = []
    if output_folder_exists:
        os.chdir(expected_dir)
        expected_files = [
            "Li3Si2-swap-NaP_intermediates",
            "Li3Sn4-swap-Na3Zn4-ReOs-OQMD_759599",
            "Li-swap-Na-edgecase-CollCode10101",
        ]
        target_elems = ["Li", "Si", "Sn"]

        for files in expected_files:
            cell, s = cell2dict(files, db=False, positions=True)
            successes.append(s)
            if s:
                elems = set(cell["atom_types"])
                elem_successes.append(any([elem not in target_elems for elem in elems]))
            else:
                elem_successes.append(False)
            res, s = res2dict(files, db=False)
            successes.append(s)
            if s:
                elems = set(res["atom_types"])
                elem_successes.append(any([elem not in target_elems for elem in elems]))
            else:
                elem_successes.append(False)

    os.chdir(OUTPUT_DIR)

    return output_folder_exists, successes, elem_successes


def changes():
    """Test matador changes functionality by undoing the second
    change (i.e. the res import).
    """
    sys.argv = ["matador", "changes", "--db", DB_NAME, "-c", "2", "--undo"]

    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]

    matador.cli.cli.main(no_quickstart=True)

    query_1 = DBQuery(db=DB_NAME, config=CONFIG_FNAME)
    query_2 = DBQuery(db=DB_NAME, composition="KSnP", config=CONFIG_FNAME)
    changes_count = MONGO_CLIENT.crystals["__changelog_" + DB_NAME].count_documents({})

    return query_1, query_2, changes_count


def refine():
    """Test various matador refine tasks."""
    sys.argv = [
        "matador",
        "refine",
        "--db",
        DB_NAME,
        "--task",
        "sym",
        "--mode",
        "overwrite",
    ]
    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]
    matador.cli.cli.main(no_quickstart=True)

    sys.argv = [
        "matador",
        "refine",
        "--db",
        DB_NAME,
        "--task",
        "source",
        "--mode",
        "set",
    ]
    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]
    matador.cli.cli.main(no_quickstart=True)

    sys.argv = [
        "matador",
        "refine",
        "--db",
        DB_NAME,
        "--task",
        "doi",
        "--mode",
        "set",
        "--new_doi",
        "10/12345",
    ]
    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]
    matador.cli.cli.main(no_quickstart=True)

    sys.argv = [
        "matador",
        "refine",
        "--db",
        DB_NAME,
        "--task",
        "tag",
        "--mode",
        "overwrite",
        "--new_tag",
        "integration_test",
    ]
    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]
    matador.cli.cli.main(no_quickstart=True)

    sys.argv = [
        "matador",
        "refine",
        "--db",
        DB_NAME,
        "--task",
        "pspot",
        "--mode",
        "overwrite",
    ]
    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]
    matador.cli.cli.main(no_quickstart=True)

    sys.argv = [
        "matador",
        "refine",
        "--db",
        DB_NAME,
        "--task",
        "raw",
        "--mode",
        "overwrite",
    ]
    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]
    matador.cli.cli.main(no_quickstart=True)

    query = DBQuery(db=DB_NAME, config=CONFIG_FNAME)

    return query


def export():
    """Test exporting to some random file types. Don't worry too much about
    the contents of the files yet, just that they exist with non-zero size.
    """
    sys.argv = [
        "matador",
        "query",
        "--db",
        DB_NAME,
        "--xsf",
        "--pdb",
        "--latex",
        "--json",
    ]

    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]

    matador.cli.cli.main(no_quickstart=True)


def uniq():
    """Test filtering of structures by PDF."""
    sys.argv = ["matador", "query", "--db", DB_NAME, "-u 0.1", "--res", "-c", "LaLiZrO"]

    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]

    matador.cli.cli.main(no_quickstart=True)


def id_query():
    """Test a simple ID query, and that nothing is found..."""
    sys.argv = ["matador", "query", "--db", DB_NAME, "-i", "testing testing"]

    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]

    matador.cli.cli.main(no_quickstart=True)

    query = DBQuery(db=DB_NAME, config=CONFIG_FNAME, id="testing testing")
    return query


def drop(collname):
    """Drop collection."""
    sys.argv = ["matador", "stats", "--delete", "--db", collname]

    if CONFIG_FNAME is not None:
        sys.argv += ["--config", CONFIG_FNAME]

    matador.cli.cli.main(no_quickstart=True)
    print("Dropped ci_test!")
