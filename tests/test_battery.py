""" Test battery routines. """

import unittest

import json
from os.path import realpath
from glob import glob

import numpy as np

from matador.battery import Electrode
from matador.hull import QueryConvexHull
from matador.scrapers.castep_scrapers import res2dict

# grab abs path for accessing test data
REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"


class VoltageTest(unittest.TestCase):
    """ Test Electrode voltage functionality. """

    def test_average_voltage(self):
        capacities = [0, 10, 100]
        voltages = [1, 1, 1]
        self.assertEqual(Electrode.calculate_average_voltage(capacities, voltages), 1.0)

        capacities = [0, 10, 100, np.nan]
        voltages = [1, 1, 1, 0]
        self.assertEqual(Electrode.calculate_average_voltage(capacities, voltages), 1.0)

        capacities = [0, 100, 1000]
        voltages = [1, 1, 0.1]
        self.assertEqual(Electrode.calculate_average_voltage(capacities, voltages), 0.19)

        capacities = [0, 100, 1000, np.nan]
        voltages = [1, 1, 0.1, 0]
        self.assertEqual(Electrode.calculate_average_voltage(capacities, voltages), 0.19)

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
        self.assertAlmostEqual(bare_hull.voltage_data["average_voltage"][0], 0.949184, places=3)
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
        self.assertTrue(np.isnan(hull.voltage_data["Q"][0][-1]))
        self.assertEqual(0.0, hull.voltage_data["voltages"][0][-1])

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
        self.assertAlmostEqual(hull.voltage_data["average_voltage"][0], 0.58523, places=4)
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
        cursor = res2dict(REAL_PATH + "data/hull-LiSiP/*.res")[0]
        hull = QueryConvexHull(cursor=cursor, species="LiSiP", no_plot=True, subcmd="voltage")
        self.assertEqual(len(hull.voltage_data["voltages"]), len(hull.voltage_data["Q"]))
        np.testing.assert_array_almost_equal(
            np.asarray(hull.voltage_data["voltages"][0]),
            np.asarray([1.1683, 1.1683, 1.0759, 0.7983, 0.6447, 0.3726, 0.3394, 0.1995, 0.1570, 0.1113, 0.1041, 0.0000]),
            decimal=3,
        )
        self.assertEqual(len(hull.voltage_data["Q"][0]), 12)
        self.assertEqual(len(hull.voltage_data["Q"][1]), 11)
        for ind in range(len(hull.voltage_data["voltages"])):
            self.assertEqual(len(hull.voltage_data["voltages"][ind]), len(hull.voltage_data["Q"][ind]))
            self.assertTrue(np.isnan(hull.voltage_data["Q"][ind][-1]))
            self.assertTrue(hull.voltage_data["voltages"][ind][-1] == 0)

    def test_ternary_voltage_problematic(self):
        """ Test for NaSnP voltages which triggered a bug in the capacity
        calculation.

        """
        cursor = res2dict(REAL_PATH + "data/voltage_data/voltage-NaSnP/*.res")[0]
        hull = QueryConvexHull(cursor=cursor, species="NaSnP", no_plot=True, subcmd="voltage")
        self.assertEqual(len(hull.voltage_data["voltages"]), 2)
        self.assertEqual(len(hull.voltage_data["Q"]), 2)
        self.assertEqual(len(hull.voltage_data["Q"][0]), 12)
        self.assertEqual(len(hull.voltage_data["Q"][1]), 13)

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

        for ind in range(len(hull.voltage_data["voltages"])):
            self.assertEqual(len(hull.voltage_data["voltages"][ind]), len(hull.voltage_data["Q"][ind]))
            self.assertTrue(np.isnan(hull.voltage_data["Q"][ind][-1]))
            self.assertTrue(hull.voltage_data["voltages"][ind][-1] == 0)



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
