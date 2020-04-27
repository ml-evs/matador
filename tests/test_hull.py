""" Test hull routines. """

import unittest
import copy
from os.path import realpath
from glob import glob

import numpy as np

from matador.hull import QueryConvexHull
from matador.query import DBQuery
from matador.scrapers.castep_scrapers import res2dict, castep2dict
from matador.utils.chem_utils import get_concentration

# grab abs path for accessing test data
REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"


class HullTest(unittest.TestCase):
    """ Test Convex hull functionality. """

    def test_hull_from_file(self):
        """ Loading hull structures from files. """
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(
            len(res_list),
            87,
            "Could not find test res files, please check installation...",
        )
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["K", "Sn", "P"], no_plot=True)
        self.assertEqual(len(hull.hull_cursor), 16)

    def test_hull_from_custom_chempots(self):
        """ Loading hull structures from files. """
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(
            len(res_list),
            87,
            "Could not find test res files, please check installation...",
        )
        cursor = [res2dict(res)[0] for res in res_list]
        cursor = [doc for doc in cursor if len(doc["stoichiometry"]) != 1]
        hull = QueryConvexHull(
            cursor=cursor,
            elements=["K", "Sn", "P"],
            no_plot=True,
            chempots=[-791.456765, -928.045026 / 2.0, 878.326441 / 4.0],
        )
        self.assertEqual(len(hull.hull_cursor), 16)

    def test_hull_with_crystal_models(self):
        """ Loading hull structures from files. """
        cursor, failures = res2dict(
            REAL_PATH + "data/hull-KPSn-KP/*.res", as_model=True
        )
        self.assertEqual(
            len(cursor),
            87,
            "Could not find all test res files, please check installation...",
        )
        cursor = [doc for doc in cursor if doc.num_elements != 1]
        hull = QueryConvexHull(
            cursor=cursor,
            elements=["K", "Sn", "P"],
            no_plot=True,
            chempots=[-791.456765, -928.045026 / 2.0, 878.326441 / 4.0],
        )
        self.assertEqual(len(hull.hull_cursor), 16)

    def test_hull_from_file_with_extraneous_elements(self):
        """ Loading hull structures from files with too many elements. """
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(
            cursor=cursor, elements=["K", "Sn"], no_plot=True, debug=True
        )
        self.assertEqual(len(hull.hull_cursor), 5)

    def test_binary_hull_distances(self):
        """ Test computing binary hull distances. """
        res_list = glob(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")
        self.assertEqual(
            len(res_list),
            295,
            "Could not find test res files, please check installation...",
        )
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["K", "P"], no_plot=True)

        hull_dist_test = np.loadtxt(REAL_PATH + "data/test_KP_hull_dist.dat")
        np.testing.assert_array_almost_equal(
            np.sort(hull_dist_test), np.sort(hull.hull_dist), decimal=3
        )

    def test_ternary_hull_distances(self):
        """ Test computing ternary hull distances. """
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(
            len(res_list),
            87,
            "Could not find test res files, please check installation...",
        )
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["K", "Sn", "P"], no_plot=True)
        self.assertEqual(len(hull.hull_cursor), 16)
        self.assertEqual(len(hull.cursor), 87)
        for ind, doc in enumerate(hull.cursor):
            hull.cursor[ind]["filename"] = doc["source"][0].split("/")[-1]

        structures = np.loadtxt(REAL_PATH + "data/test_KSnP.dat")
        hull_dist_test = np.loadtxt(REAL_PATH + "data/test_KSnP_hull_dist.dat")
        precomp_hull_dist = hull.get_hull_distances(structures, precompute=True)
        no_precomp_hull_dist = hull.get_hull_distances(structures, precompute=False)
        np.testing.assert_array_almost_equal(
            np.sort(hull_dist_test), np.sort(hull.hull_dist), decimal=3
        )
        np.testing.assert_array_almost_equal(
            np.sort(hull.hull_dist), np.sort(precomp_hull_dist), decimal=3
        )
        np.testing.assert_array_almost_equal(
            no_precomp_hull_dist, precomp_hull_dist, decimal=5
        )
        self.assertFalse(np.isnan(hull.hull_dist).any())

    def test_toy_ternary(self):
        cursor = [
            {
                "stoichiometry": [["K", 1.0]],
                "enthalpy_per_atom": 0,
                "cell_volume": 100,
                "num_atoms": 10,
            },
            {
                "stoichiometry": [["Sn", 1.0]],
                "enthalpy_per_atom": 0,
                "cell_volume": 100,
                "num_atoms": 10,
            },
            {
                "stoichiometry": [["P", 1.0]],
                "enthalpy_per_atom": 0,
                "cell_volume": 100,
                "num_atoms": 10,
            },
            {
                "stoichiometry": [["K", 1.0], ["Sn", 1.0], ["P", 1.0]],
                "enthalpy_per_atom": -1,
                "cell_volume": 100,
                "num_atoms": 30,
            },
            {
                "stoichiometry": [["K", 1.0], ["Sn", 1.0], ["P", 1.0]],
                "enthalpy_per_atom": -0.5,
                "cell_volume": 100,
                "num_atoms": 30,
            },
            {
                "stoichiometry": [["Sn", 1.0], ["P", 1.0]],
                "enthalpy_per_atom": 0.1,
                "cell_volume": 100,
                "num_atoms": 30,
            },
            {
                "stoichiometry": [["K", 4.0], ["Sn", 2.0], ["P", 2.0]],
                "enthalpy_per_atom": -0.74,
                "cell_volume": 100,
                "num_atoms": 30,
            },
            {
                "stoichiometry": [["K", 2.0], ["Sn", 2.0], ["P", 4.0]],
                "enthalpy_per_atom": -0.74,
                "cell_volume": 100,
                "num_atoms": 30,
            },
            {
                "stoichiometry": [["K", 2.0], ["Sn", 4.0], ["P", 2.0]],
                "enthalpy_per_atom": -0.74,
                "cell_volume": 100,
                "num_atoms": 30,
            },
            {
                "stoichiometry": [["Sn", 1.0], ["P", 1.0]],
                "enthalpy_per_atom": 0,
                "cell_volume": 100,
                "num_atoms": 30,
            },
        ]
        for ind, doc in enumerate(cursor):
            cursor[ind]["concentration"] = get_concentration(doc, ["K", "Sn", "P"])
            cursor[ind]["source"] = ["abcde"]
            cursor[ind]["num_fu"] = 10
            cursor[ind]["enthalpy"] = (
                cursor[ind]["enthalpy_per_atom"] * doc["num_atoms"]
            )

        hull = QueryConvexHull(cursor=cursor, elements=["K", "Sn", "P"], no_plot=True)

        hull_dists = [doc["hull_distance"] for doc in hull.cursor]
        np.testing.assert_array_almost_equal(
            hull_dists, [0, 0, 0, 0, 0.5, 0.1, 0.01, 0.01, 0.01, 0]
        )

    def test_pseudoternary_hull(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        print()
        print(80 * "-")
        self.assertEqual(
            len(cursor), 12, "Error with test res files, please check installation..."
        )
        hull = QueryConvexHull(
            cursor=cursor, elements=["La2O3", "ZrO2", "Li2O"], no_plot=True
        )
        self.assertEqual(len(hull.cursor), 7)
        self.assertEqual(len(hull.hull_cursor), 5)
        hull = QueryConvexHull(
            cursor=cursor,
            elements=["La2O3", "ZrO2", "Li2O"],
            hull_cutoff=0.01,
            no_plot=True,
        )
        self.assertEqual(len(hull.cursor), 7)
        self.assertEqual(len(hull.hull_cursor), 6)
        self.assertEqual(len(hull.convex_hull.vertices), 5)
        for doc in hull.hull_cursor:
            if "cubic-LLZO" in doc["source"][0]:
                self.assertAlmostEqual(
                    doc["formation_enthalpy_per_atom"], -0.05758265622
                )
                self.assertAlmostEqual(doc["hull_distance"], 0.008746875)
                break
        else:
            raise RuntimeError("Did not find cubic-LLZO in cursor")

    def test_pseudoternary_hull_NaVSO4(self):
        cursor, s = castep2dict(REAL_PATH + "data/hull-NaVSO4/*.castep")
        print(80 * "-")
        self.assertEqual(
            len(cursor), 9, "Error with test castep files, please check installation..."
        )
        hull = QueryConvexHull(
            cursor=cursor,
            elements=["Na", "V", "SO4"],
            chempots=[-999, -999, 100],
            no_plot=True,
        )
        self.assertEqual(len(hull.cursor), 8)
        self.assertEqual(len(hull.hull_cursor), 6)

    def test_pseudoternary_hull_with_custom_chempots(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        print()
        print(80 * "-")
        self.assertEqual(
            len(cursor), 12, "Error with test res files, please check installation..."
        )
        hull = QueryConvexHull(
            cursor=cursor,
            elements=["La2O3", "ZrO2", "Li2O"],
            chempots=[-3275.040 / 5, -2178.987 / 3, -848.148 / 3],
            no_plot=True,
        )
        self.assertEqual(len(hull.cursor), 10)
        self.assertEqual(len(hull.hull_cursor), 5)

        # cursor, s = res2dict(REAL_PATH + 'data/hull-LLZO/*.res')
        hull = QueryConvexHull(
            cursor=cursor,
            elements=["La2O3", "ZrO2", "Li2O"],
            hull_cutoff=0.20,
            chempots=[26200.3194 / 40, -8715.94784 / 12, -3392.59361 / 12],
            no_plot=True,
        )
        self.assertEqual(len(hull.cursor), 10)
        self.assertEqual(len(hull.hull_cursor), 10)

        for doc in hull.hull_cursor:
            if "cubic-LLZO" in doc["source"][0]:
                self.assertAlmostEqual(
                    doc["formation_enthalpy_per_atom"], -0.05758265622
                )
                self.assertAlmostEqual(doc["hull_distance"], 0.008746875)
                break
        else:
            raise RuntimeError("Did not find cubic-LLZO in cursor")

    def test_pseudoternary_hull_with_custom_chempots_below_hull(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        hull = QueryConvexHull(
            cursor=cursor,
            elements=["La2O3", "ZrO2", "Li2O"],
            hull_cutoff=0.02,
            chempots=[26200.3194 / 40, -8715.94784 / 12, -3392.59361 / 12],
            no_plot=True,
        )
        self.assertEqual(len(hull.cursor), 10)
        self.assertEqual(len(hull.hull_cursor), 9)
        for doc in cursor:
            if "La2O3" in doc["source"][0]:
                new_doc = copy.deepcopy(doc)
                new_doc["enthalpy_per_atom"] -= 1
                new_doc["source"][0] = "IMAGINARY.res"
        hull_dists = hull.get_hull_distances([[1, 0, -1]])
        self.assertAlmostEqual(hull_dists[0], -1)

    def test_hull_pseudobinary_hull(self):
        cursor, s = res2dict(REAL_PATH + "data/query-SrTiO-oqmd_1.1/*.res")
        print()
        print(80 * "-")
        self.assertEqual(
            len(cursor), 202, "Error with test res files, please check installation..."
        )
        hull = QueryConvexHull(cursor=cursor, species=["Sr", "TiO3"], no_plot=True)

        self.assertEqual(len(hull.cursor), 37)
        self.assertEqual(len(hull.hull_cursor), 3)

        hull = QueryConvexHull(
            cursor=cursor, species=["Sr", "TiO3"], hull_cutoff=0.05, no_plot=True
        )
        self.assertEqual(len(hull.cursor), 37)
        self.assertEqual(len(hull.hull_cursor), 10)

    def test_pseudoternary_from_fake_query(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        print()
        print(80 * "-")
        self.assertEqual(
            len(cursor), 12, "Error with test res files, please check installation..."
        )
        hull = QueryConvexHull(
            cursor=cursor, elements=["La2O3", "ZrO2", "Li2O"], no_plot=True
        )
        self.assertEqual(len(hull.cursor), 7)

        fake_query = DBQuery.__new__(DBQuery)
        fake_query.cursor = hull.cursor
        for ind, doc in enumerate(fake_query.cursor):
            fake_query.cursor[ind]["_id"] = None
            fake_query.cursor[ind]["text_id"] = [doc["source"][0], "."]
        fake_query._non_elemental = True
        fake_query._create_hull = True
        fake_query.args = dict()
        fake_query.args["intersection"] = True
        fake_query.args["subcmd"] = "hull"
        fake_query.args["composition"] = ["La2O3:ZrO2:Li2O"]
        hull = QueryConvexHull(
            query=fake_query,
            hull_cutoff=0.01,
            chempots=[26200.3194 / 40, -8715.94784 / 12, -3392.59361 / 12],
            no_plot=True,
        )
        self.assertEqual(len(hull.cursor), 10)
        self.assertEqual(len(hull.hull_cursor), 9)

    def test_variations_of_element_parsing(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")
        self.assertEqual(
            len(cursor), 295, "Error with test res files, please check installation..."
        )
        msg = "failed to parse elements=['KP']"
        hull = QueryConvexHull(cursor=cursor, elements=["KP"], no_plot=True)
        self.assertEqual(len(hull.cursor), 295, msg=msg)
        self.assertEqual(len(hull.hull_cursor), 7, msg=msg)

        msg = "failed to parse species=['KP']"
        hull = QueryConvexHull(cursor=cursor, species=["KP"], no_plot=True)
        self.assertEqual(len(hull.cursor), 295, msg=msg)
        self.assertEqual(len(hull.hull_cursor), 7, msg=msg)

        msg = "failed to parse elements=['K', 'P']"
        hull = QueryConvexHull(cursor=cursor, elements=["K", "P"], no_plot=True)
        self.assertEqual(len(hull.cursor), 295, msg=msg)
        self.assertEqual(len(hull.hull_cursor), 7, msg=msg)

        msg = "failed to parse species=['K', 'P']"
        hull = QueryConvexHull(cursor=cursor, species=["K", "P"], no_plot=True)
        self.assertEqual(len(hull.cursor), 295, msg=msg)
        self.assertEqual(len(hull.hull_cursor), 7, msg=msg)

        msg = "failed to parse species='KP'"
        hull = QueryConvexHull(cursor=cursor, species="KP", no_plot=True)
        self.assertEqual(len(hull.cursor), 295, msg=msg)
        self.assertEqual(len(hull.hull_cursor), 7, msg=msg)

        msg = "failed to parse elements='KP'"
        hull = QueryConvexHull(cursor=cursor, elements="KP", no_plot=True)
        self.assertEqual(len(hull.cursor), 295, msg=msg)
        self.assertEqual(len(hull.hull_cursor), 7, msg=msg)

        msg = "failed to parse elements='K:P'"
        hull = QueryConvexHull(cursor=cursor, elements="K:P", no_plot=True)
        self.assertEqual(len(hull.cursor), 295, msg=msg)
        self.assertEqual(len(hull.hull_cursor), 7, msg=msg)

        msg = "failed to parse species='K:P'"
        hull = QueryConvexHull(cursor=cursor, species="K:P", no_plot=True)
        self.assertEqual(len(hull.cursor), 295, msg=msg)
        self.assertEqual(len(hull.hull_cursor), 7, msg=msg)

    def test_binary_from_fake_query(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")
        print()
        print(80 * "-")
        self.assertEqual(
            len(cursor), 295, "Error with test res files, please check installation..."
        )
        hull = QueryConvexHull(cursor=cursor, elements=["KP"], no_plot=True)
        self.assertEqual(len(hull.cursor), 295)
        self.assertEqual(len(hull.hull_cursor), 7)

        fake_query = DBQuery.__new__(DBQuery)
        fake_query.cursor = hull.cursor
        for ind, doc in enumerate(fake_query.cursor):
            fake_query.cursor[ind]["_id"] = None
            fake_query.cursor[ind]["text_id"] = [doc["source"][0], "."]
            del fake_query.cursor[ind]["hull_distance"]
            del fake_query.cursor[ind]["concentration"]

        fake_query._non_elemental = False
        fake_query._create_hull = True
        fake_query.args = dict()
        fake_query.args["intersection"] = False
        fake_query.args["subcmd"] = "hull"
        fake_query.args["composition"] = ["KP"]
        hull = QueryConvexHull(
            query=fake_query, chempots=[-791.456765, -219.58161025], no_plot=True
        )
        # now need to include custom chempots in counts
        self.assertEqual(len(hull.cursor), 297)
        self.assertEqual(len(hull.hull_cursor), 9)

    def test_filter_cursor(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        species = ["Li2O", "La2O3", "ZrO2"]
        new_cursor = QueryConvexHull.filter_cursor_by_chempots(species, cursor)
        self.assertEqual(len(cursor), 12)
        self.assertEqual(len(new_cursor), 7)
        cursor = [{"stoichiometry": [["Li", 7], ["O", 12], ["Zr", 2], ["La", 3]]}]
        new_cursor = QueryConvexHull.filter_cursor_by_chempots(species, cursor)
        self.assertEqual(len(new_cursor), 1)
        self.assertAlmostEqual(
            new_cursor[0]["concentration"][0], 0.5, msg="Concentrations do not match"
        )
        self.assertAlmostEqual(
            new_cursor[0]["concentration"][1],
            1.5 / 7.0,
            msg="Concentrations do not match",
        )

    def test_pseudoternary_hull_failure(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        print()
        print(80 * "-")
        self.assertEqual(
            len(cursor), 12, "Error with test res files, please check installation..."
        )
        with self.assertRaises(RuntimeError):
            QueryConvexHull(cursor=cursor, no_plot=True)


class EnsembleHullTest(unittest.TestCase):
    """ Test of Ensemble Hulls for BEEF/temperature. """

    def test_beef_hull(self):
        from matador.hull import EnsembleHull
        from matador.scrapers import castep2dict

        cursor, s = castep2dict(REAL_PATH + "data/beef_files/*.castep", db=False)

        beef_hull = EnsembleHull(
            cursor, "_beef", energy_key="total_energy_per_atom", parameter_key="thetas"
        )

        self.assertEqual(len(beef_hull.phase_diagrams), 5000)
        self.assertEqual(len(beef_hull.cursor[0]["_beef"]["hull_distance"]), 5000)
        self.assertEqual(len(beef_hull.cursor[1]["_beef"]["hull_distance"]), 5000)
