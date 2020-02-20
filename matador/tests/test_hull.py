#!/usr/bin/env python

""" Test hull routines. """

import unittest
import json
import copy
from os.path import realpath
from glob import glob

import numpy as np

from matador.hull import QueryConvexHull
from matador.query import DBQuery
from matador.scrapers.castep_scrapers import res2dict, castep2dict
from matador.export import generate_hash

# grab abs path for accessing test data
REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"


class HullTest(unittest.TestCase):
    """ Test Convex hull functionality. """

    def test_hull_from_file(self):
        """ Loading hull structures from files. """
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(len(res_list), 87, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["K", "Sn", "P"], no_plot=True)
        self.assertEqual(len(hull.hull_cursor), 16)

    def test_hull_from_custom_chempots(self):
        """ Loading hull structures from files. """
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(len(res_list), 87, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        cursor = [doc for doc in cursor if len(doc["stoichiometry"]) != 1]
        hull = QueryConvexHull(
            cursor=cursor, elements=["K", "Sn", "P"], no_plot=True, chempots=[-791.456765, -928.045026 / 2.0, 878.326441 / 4.0]
        )
        self.assertEqual(len(hull.hull_cursor), 16)

    def test_hull_from_file_with_extraneous_elements(self):
        """ Loading hull structures from files with too many elements. """
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["K", "Sn"], no_plot=True, debug=True)
        self.assertEqual(len(hull.hull_cursor), 5)

    def test_binary_hull_distances(self):
        """ Test computing binary hull distances. """
        res_list = glob(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")
        self.assertEqual(len(res_list), 295, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["K", "P"], no_plot=True)

        hull_dist_test = np.loadtxt(REAL_PATH + "data/test_KP_hull_dist.dat")
        np.testing.assert_array_almost_equal(np.sort(hull_dist_test), np.sort(hull.hull_dist), decimal=3)

    def test_ternary_hull_distances(self):
        """ Test computing ternary hull distances. """
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(len(res_list), 87, "Could not find test res files, please check installation...")
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
        np.testing.assert_array_almost_equal(np.sort(hull_dist_test), np.sort(hull.hull_dist), decimal=3)
        np.testing.assert_array_almost_equal(np.sort(hull.hull_dist), np.sort(precomp_hull_dist), decimal=3)
        np.testing.assert_array_almost_equal(no_precomp_hull_dist, precomp_hull_dist, decimal=5)

    def test_pseudoternary_hull(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        print()
        print(80 * "-")
        self.assertEqual(len(cursor), 12, "Error with test res files, please check installation...")
        hull = QueryConvexHull(cursor=cursor, elements=["La2O3", "ZrO2", "Li2O"], no_plot=True)
        self.assertEqual(len(hull.cursor), 7)
        self.assertEqual(len(hull.hull_cursor), 5)
        hull = QueryConvexHull(cursor=cursor, elements=["La2O3", "ZrO2", "Li2O"], hull_cutoff=0.01, no_plot=True)
        self.assertEqual(len(hull.cursor), 7)
        self.assertEqual(len(hull.hull_cursor), 6)
        self.assertEqual(len(hull.convex_hull.vertices), 5)
        for doc in hull.hull_cursor:
            if "cubic-LLZO" in doc["source"][0]:
                self.assertAlmostEqual(doc["formation_enthalpy_per_atom"], -0.05758265622)
                self.assertAlmostEqual(doc["hull_distance"], 0.008746875)
                break
        else:
            raise RuntimeError("Did not find cubic-LLZO in cursor")

    def test_pseudoternary_hull_NaVSO4(self):
        cursor, s = castep2dict(REAL_PATH + "data/hull-NaVSO4/*.castep")
        print(80 * "-")
        self.assertEqual(len(cursor), 9, "Error with test castep files, please check installation...")
        hull = QueryConvexHull(cursor=cursor, elements=["Na", "V", "SO4"], chempots=[-999, -999, 100], no_plot=True)
        self.assertEqual(len(hull.cursor), 8)
        self.assertEqual(len(hull.hull_cursor), 6)

    def test_pseudoternary_hull_with_custom_chempots(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        print()
        print(80 * "-")
        self.assertEqual(len(cursor), 12, "Error with test res files, please check installation...")
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
                self.assertAlmostEqual(doc["formation_enthalpy_per_atom"], -0.05758265622)
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
        self.assertEqual(len(cursor), 202, "Error with test res files, please check installation...")
        hull = QueryConvexHull(cursor=cursor, species=["Sr", "TiO3"], no_plot=True)

        self.assertEqual(len(hull.cursor), 37)
        self.assertEqual(len(hull.hull_cursor), 3)

        hull = QueryConvexHull(cursor=cursor, species=["Sr", "TiO3"], hull_cutoff=0.05, no_plot=True)
        self.assertEqual(len(hull.cursor), 37)
        self.assertEqual(len(hull.hull_cursor), 10)

    def test_pseudoternary_from_fake_query(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        print()
        print(80 * "-")
        self.assertEqual(len(cursor), 12, "Error with test res files, please check installation...")
        hull = QueryConvexHull(cursor=cursor, elements=["La2O3", "ZrO2", "Li2O"], no_plot=True)
        self.assertEqual(len(hull.cursor), 7)

        fake_query = DBQuery.__new__(DBQuery)
        fake_query.cursor = hull.cursor
        for ind, doc in enumerate(fake_query.cursor):
            fake_query.cursor[ind]["_id"] = generate_hash(hash_len=20)
            fake_query.cursor[ind]["text_id"] = [doc["source"][0], "."]
        fake_query._non_elemental = True
        fake_query._create_hull = True
        fake_query.args = dict()
        fake_query.args["intersection"] = True
        fake_query.args["subcmd"] = "hull"
        fake_query.args["composition"] = ["La2O3:ZrO2:Li2O"]
        hull = QueryConvexHull(
            query=fake_query, hull_cutoff=0.01, chempots=[26200.3194 / 40, -8715.94784 / 12, -3392.59361 / 12], no_plot=True
        )
        self.assertEqual(len(hull.cursor), 10)
        self.assertEqual(len(hull.hull_cursor), 9)

    def test_variations_of_element_parsing(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")
        self.assertEqual(len(cursor), 295, "Error with test res files, please check installation...")
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
        self.assertEqual(len(cursor), 295, "Error with test res files, please check installation...")
        hull = QueryConvexHull(cursor=cursor, elements=["KP"], no_plot=True)
        self.assertEqual(len(hull.cursor), 295)
        self.assertEqual(len(hull.hull_cursor), 7)

        fake_query = DBQuery.__new__(DBQuery)
        fake_query.cursor = hull.cursor
        for ind, doc in enumerate(fake_query.cursor):
            fake_query.cursor[ind]["_id"] = generate_hash(hash_len=20)
            fake_query.cursor[ind]["text_id"] = [doc["source"][0], "."]
            del fake_query.cursor[ind]["hull_distance"]
            del fake_query.cursor[ind]["concentration"]

        fake_query._non_elemental = False
        fake_query._create_hull = True
        fake_query.args = dict()
        fake_query.args["intersection"] = False
        fake_query.args["subcmd"] = "hull"
        fake_query.args["composition"] = ["KP"]
        hull = QueryConvexHull(query=fake_query, chempots=[-791.456765, -219.58161025], no_plot=True)
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
        self.assertAlmostEqual(new_cursor[0]["concentration"][0], 0.5, msg="Concentrations do not match")
        self.assertAlmostEqual(new_cursor[0]["concentration"][1], 1.5 / 7.0, msg="Concentrations do not match")

    def test_pseudoternary_hull_failure(self):
        cursor, s = res2dict(REAL_PATH + "data/hull-LLZO/*.res")
        print()
        print(80 * "-")
        self.assertEqual(len(cursor), 12, "Error with test res files, please check installation...")
        with self.assertRaises(RuntimeError):
            QueryConvexHull(cursor=cursor, no_plot=True)


class EnsembleHullTest(unittest.TestCase):
    """ Test of Ensemble Hulls for BEEF/temperature. """

    def test_beef_hull(self):
        from matador.hull import EnsembleHull
        from matador.scrapers import castep2dict

        cursor, s = castep2dict(REAL_PATH + "data/beef_files/*.castep", db=False)

        beef_hull = EnsembleHull(cursor, "_beef", energy_key="total_energy_per_atom", parameter_key="thetas")

        self.assertEqual(len(beef_hull.phase_diagrams), 5000)
        self.assertEqual(len(beef_hull.cursor[0]["_beef"]["hull_distance"]), 5000)
        self.assertEqual(len(beef_hull.cursor[1]["_beef"]["hull_distance"]), 5000)


class VoltageTest(unittest.TestCase):
    """ Test voltage curve functionality. """

    def test_binary_voltage(self):
        """ Test simple binary voltage curve. """
        hull_cursor = []
        test_x = np.loadtxt(REAL_PATH + "data/voltage_data/LiAs_x.dat")
        test_Q = np.loadtxt(REAL_PATH + "data/voltage_data/LiAs_Q.dat")
        test_V = np.loadtxt(REAL_PATH + "data/voltage_data/LiAs_V.dat")
        for i in range(5):
            with open(REAL_PATH + "data/voltage_data/hull_data" + str(i) + ".json") as f:
                hull_cursor.append(json.load(f))
        bare_hull = QueryConvexHull(cursor=hull_cursor, species=["Li", "As"], subcmd="voltage", no_plot=True)
        self.assertTrue(len(bare_hull.voltage_data["voltages"]) == 1)
        np.testing.assert_array_almost_equal(bare_hull.voltage_data["voltages"][0], test_V, decimal=5, verbose=True)
        np.testing.assert_array_almost_equal(bare_hull.voltage_data["x"][0], test_x, decimal=5)
        np.testing.assert_array_almost_equal(bare_hull.voltage_data["Q"][0], test_Q, decimal=5)
        for ind in range(len(bare_hull.voltage_data["voltages"])):
            self.assertTrue(np.isnan(bare_hull.voltage_data["Q"][ind][-1]))
            self.assertTrue(bare_hull.voltage_data["voltages"][ind][-1] == 0)

    def test_binary_voltage_mayo(self):
        """ Test binary voltages from cursor for Mayo et al,
        DOI: 10.1021/acs.chemmater.5b04208.

        """
        res_list = glob(REAL_PATH + "data/hull-LiP-mdm_chem_mater/*.res")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["Li", "P"], no_plot=True, subcmd="voltage")
        self.assertEqual(len(hull.voltage_data["voltages"]), len(hull.voltage_data["Q"]))
        self.assertEqual(len(hull.voltage_data["voltages"]), len(hull.voltage_data["Q"]))
        LiP_voltage_curve = np.loadtxt(REAL_PATH + "data/LiP_voltage.csv", delimiter=",")
        self.assertTrue(len(hull.voltage_data["voltages"]) == 1)
        np.testing.assert_allclose(hull.voltage_data["voltages"][0], LiP_voltage_curve[:, 1], verbose=True, rtol=1e-4)
        np.testing.assert_allclose(hull.voltage_data["Q"][0], LiP_voltage_curve[:, 0], verbose=True, rtol=1e-4)

    def test_ternary_voltage(self):
        """ Test ternary voltages from cursor. """
        # test data from LiSnS
        res_list = glob(REAL_PATH + "data/hull-LiSnS/*.res")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(
            cursor=cursor, elements=["Li", "Sn", "S"], no_plot=True, pathways=True, subcmd="voltage", debug=True
        )
        pin = np.array(
            [
                [2, 0, 0, -380.071],
                [0, 2, 4, -1305.0911],
                [2, 0, 1, -661.985],
                [6, 2, 0, -1333.940],
                [16, 4, 16, -7906.417],
                [4, 4, 0, -1144.827],
                [0, 4, 4, -1497.881],
                [0, 1, 0, -95.532],
                [0, 0, 48, -13343.805],
            ]
        )
        tot = pin[:, 0] + pin[:, 1] + pin[:, 2]
        points = pin / tot[:, None]

        voltage_data = [
            np.asarray(
                [
                    1.9415250000000697,
                    1.9415250000000697,
                    1.8750000000001705,
                    1.4878749999999741,
                    0.63925000000000409,
                    0.34612500000000068,
                    0.0,
                ]
            ),
            np.asarray([1.4878749999999741, 1.4878749999999741, 0.63925000000000409, 0.34612500000000068, 0.0]),
        ]

        Q_data = [np.array([0, 195, 293, 586, 733, 1026, np.NaN]), np.array([0, 356, 533, 889, np.NaN])]

        points = np.delete(points, 2, axis=1)

        self.assertEqual(len(hull.voltage_data["Q"]), len(Q_data))
        for i in range(len(hull.voltage_data["voltages"])):
            np.testing.assert_array_almost_equal(hull.voltage_data["voltages"][i], voltage_data[i], decimal=3)
            np.testing.assert_array_almost_equal(hull.voltage_data["Q"][i], Q_data[i], decimal=0)
        for ind in range(len(hull.voltage_data["voltages"])):
            self.assertEqual(len(hull.voltage_data["Q"][ind]), len(hull.voltage_data["voltages"][ind]))
            self.assertTrue(np.isnan(hull.voltage_data["Q"][ind][-1]))
            self.assertTrue(hull.voltage_data["voltages"][ind][-1] == 0)

    def test_ternary_voltage_with_one_two_phase_region(self):
        """ Test ternary voltages with awkward two-phase region. """
        # load old hull then rejig it to go through a ternary phase
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(len(res_list), 87, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        cursor = [
            doc
            for doc in cursor
            if (
                doc["stoichiometry"] != [["P", 3], ["Sn", 4]]
                and doc["stoichiometry"] != [["P", 3], ["Sn", 1]]
                and doc["stoichiometry"] != [["K", 3], ["P", 7]]
                and doc["stoichiometry"] != [["K", 1], ["P", 7]]
                and doc["stoichiometry"] != [["K", 2], ["P", 3]]
                and doc["stoichiometry"] != [["K", 8], ["P", 4], ["Sn", 1]]
                and doc["stoichiometry"] != [["K", 1], ["P", 2], ["Sn", 2]]
                and doc["stoichiometry"] != [["K", 1], ["Sn", 1]]
                and doc["stoichiometry"] != [["K", 4], ["Sn", 9]]
                and doc["stoichiometry"] != [["K", 5], ["P", 4]]
                and doc["stoichiometry"] != [["P", 2], ["Sn", 1]]
            )
        ]
        hull = QueryConvexHull(cursor=cursor, elements=["K", "Sn", "P"], no_plot=True, pathways=True, subcmd="voltage")
        self.assertEqual(len(hull.voltage_data["voltages"]), len(hull.voltage_data["Q"]))
        np.testing.assert_array_almost_equal(
            np.asarray(hull.voltage_data["voltages"]), np.asarray([[1.0229, 1.0229, 0.2676, 0.000]]), decimal=3
        )
        for ind in range(len(hull.voltage_data["voltages"])):
            self.assertEqual(len(hull.voltage_data["Q"][ind]), len(hull.voltage_data["voltages"][ind]))
            self.assertTrue(np.isnan(hull.voltage_data["Q"][ind][-1]))
            self.assertTrue(hull.voltage_data["voltages"][ind][-1] == 0)

    def test_ternary_voltage_with_two_two_phase_regions(self):
        """ Test ternary voltages with two awkward two-phase regions. """
        # load old hull then rejig it to go through a ternary phase
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(len(res_list), 87, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        cursor = [
            doc
            for doc in cursor
            if (
                doc["stoichiometry"] != [["P", 3], ["Sn", 4]]
                and doc["stoichiometry"] != [["P", 3], ["Sn", 1]]
                and doc["stoichiometry"] != [["K", 3], ["P", 7]]
                and doc["stoichiometry"] != [["K", 1], ["P", 7]]
                and doc["stoichiometry"] != [["K", 2], ["P", 3]]
                and doc["stoichiometry"] != [["K", 8], ["P", 4], ["Sn", 1]]
                and doc["stoichiometry"] != [["K", 1], ["Sn", 1]]
                and doc["stoichiometry"] != [["K", 4], ["Sn", 9]]
                and doc["stoichiometry"] != [["K", 5], ["P", 4]]
                and doc["stoichiometry"] != [["P", 2], ["Sn", 1]]
            )
        ]
        hull = QueryConvexHull(cursor=cursor, elements=["K", "Sn", "P"], no_plot=True, pathways=True, subcmd="voltage")
        self.assertEqual(len(hull.voltage_data["voltages"]), len(hull.voltage_data["Q"]))
        np.testing.assert_array_almost_equal(
            np.asarray(hull.voltage_data["voltages"][0]), np.asarray([1.1845, 1.1845, 0.8612, 0.2676, 0.000]), decimal=3
        )
        self.assertAlmostEqual(hull.voltage_data["Q"][0][-2], 425.7847612, places=5)
        for ind in range(len(hull.voltage_data["voltages"])):
            self.assertEqual(len(hull.voltage_data["Q"][ind]), len(hull.voltage_data["voltages"][ind]))
            self.assertTrue(np.isnan(hull.voltage_data["Q"][ind][-1]))
            self.assertTrue(hull.voltage_data["voltages"][ind][-1] == 0)

    def test_ternary_voltage_with_exclusively_two_phase_regions(self):
        """ Test ternary voltages exclusively awkward two-phase regions. """
        # load old hull then rejig it to go through a ternary phase
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(len(res_list), 87, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        cursor = [
            doc
            for doc in cursor
            if (
                doc["stoichiometry"] == [["P", 1], ["Sn", 1]]
                or doc["stoichiometry"] == [["K", 1], ["P", 1], ["Sn", 1]]
                or doc["stoichiometry"] == [["P", 1]]
                or doc["stoichiometry"] == [["K", 1]]
                or doc["stoichiometry"] == [["Sn", 1]]
            )
        ]

        hull = QueryConvexHull(cursor=cursor, elements=["K", "Sn", "P"], no_plot=True, pathways=True, subcmd="voltage")
        self.assertEqual(len(hull.voltage_data["voltages"]), len(hull.voltage_data["Q"]))
        self.assertEqual(len(hull.voltage_data["voltages"]), 1)
        for ind in range(len(hull.voltage_data["voltages"])):
            self.assertEqual(len(hull.voltage_data["Q"][ind]), len(hull.voltage_data["voltages"][ind]))
            self.assertTrue(np.isnan(hull.voltage_data["Q"][ind][-1]))

    def test_ternary_voltage_with_single_phase_region(self):
        """ Test ternary voltages with single-phase regions. """
        # load old hull then rejig it to go through a ternary phase
        res_list = glob(REAL_PATH + "data/hull-LiSiP/*.res")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["Li", "Si", "P"], pathways=True, no_plot=True, subcmd="voltage")
        self.assertEqual(len(hull.voltage_data["voltages"]), len(hull.voltage_data["Q"]))
        np.testing.assert_array_almost_equal(
            np.asarray(hull.voltage_data["voltages"][0]),
            np.asarray([1.1683, 1.1683, 1.0759, 0.7983, 0.6447, 0.3726, 0.3394, 0.1995, 0.1570, 0.1113, 0.1041, 0.0000]),
            decimal=3,
        )
        for ind in range(len(hull.voltage_data["voltages"])):
            self.assertEqual(len(hull.voltage_data["voltages"][ind]), len(hull.voltage_data["Q"][ind]))
            self.assertTrue(np.isnan(hull.voltage_data["Q"][ind][-1]))
            self.assertTrue(hull.voltage_data["voltages"][ind][-1] == 0)

    def test_angelas_awkward_voltage(self):
        """ Test a particular example of Angela's awkward ternary voltages. """
        # test data from NaFeP
        res_list = glob(REAL_PATH + "data/hull-NaFeP-afh41_new_Na+Fe+P/*.res")
        self.assertEqual(len(res_list), 16, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["Na", "Fe", "P"], no_plot=True, subcmd="voltage")
        self.assertEqual(len(hull.voltage_data["voltages"][0]), 8)
        self.assertEqual(len(hull.voltage_data["voltages"][1]), 5)
        self.assertEqual(len(hull.voltage_data["voltages"][2]), 3)


class VolumeTest(unittest.TestCase):
    """ Test simple binary volume curve. """

    def test_binary_volume_curve(self):
        """ Test simple binary volume curve. """
        res_list = glob(REAL_PATH + "data/hull-LiP-mdm_chem_mater/*.res")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, elements=["Li", "P"], no_plot=True)
        hull.volume_curve()
        self.assertEqual(len(hull.volume_data["x"]), len(hull.hull_cursor) - 1)
        self.assertEqual(len(hull.volume_data["vol_per_y"]), len(hull.hull_cursor) - 1)


if __name__ == "__main__":
    unittest.main(buffer=True, verbosity=2)
