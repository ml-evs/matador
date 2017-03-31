#!/usr/bin/env python

# matador modules
from matador.hull import QueryConvexHull
from matador.utils.chem_utils import get_generic_grav_capacity
# external libraries
from scipy.spatial import ConvexHull
import numpy as np
# standard library
import sys
import os
import json
from os.path import realpath
import unittest

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class VoltageTest(unittest.TestCase):
    """ Test voltage curve functionality. """
    def testBinaryVoltage(self):
        match, hull_cursor = [], []
        test_x = np.loadtxt(REAL_PATH + 'data/x.dat')
        test_Q = np.loadtxt(REAL_PATH + 'data/Q.dat')
        test_V = np.loadtxt(REAL_PATH + 'data/V.dat')
        for i in range(5):
            with open(REAL_PATH + 'data/hull_data' + str(i) + '.json') as f:
                hull_cursor.append(json.load(f))
        for i in range(2):
            with open(REAL_PATH + 'data/mu' + str(i) + '.json') as f:
                match.append(json.load(f))
        with open(REAL_PATH + 'data/elements.json') as f:
            elements = json.load(f)
        bare_hull = QueryConvexHull.__new__(QueryConvexHull)
        bare_hull.cursor = list(hull_cursor)
        bare_hull.ternary = False
        bare_hull.elements = list(elements)
        bare_hull.hull_cursor = list(hull_cursor)
        bare_hull.match = list(match)
        with open(os.devnull, 'w') as sys.stdout:
            bare_hull.voltage_curve(bare_hull.hull_cursor)
        sys.stdout = sys.__stdout__
        self.assertTrue(len(bare_hull.voltages) == 1)
        np.testing.assert_array_equal(bare_hull.voltages[0], test_V, verbose=True)
        np.testing.assert_array_equal(bare_hull.x[0], test_x)
        np.testing.assert_array_equal(bare_hull.Q[0], test_Q)

    def testTernaryVoltage(self):
        # test data from LiSnS
        pin = np.array([[2, 0, 0, -380.071],
                        [0, 2, 4, -1305.0911],
                        [2, 0, 1, -661.985],
                        [6, 2, 0, -1333.940],
                        [16, 4, 16, -7906.417],
                        [4, 4, 0, -1144.827],
                        [0, 4, 4, -1497.881],
                        [0, 1, 0, -95.532],
                        [0, 0, 48, -13343.805]])
        tot = pin[:, 0] + pin[:, 1] + pin[:, 2]
        points = pin/tot[:, None]
        hull_cursor = []
        for ind, point in enumerate(points):
            hull_cursor.append(dict())
            hull_cursor[-1]['gravimetric_capacity'] = get_generic_grav_capacity(point[0:3],
                                                                                ['Li', 'Sn', 'S'])
            hull_cursor[-1]['stoichiometry'] = [['Li', int(pin[ind][0])],
                                                ['Sn', int(pin[ind][1])],
                                                ['S', int(pin[ind][2])]]
            hull_cursor[-1]['concentration'] = point[0:2]
            hull_cursor[-1]['enthalpy_per_atom'] = point[-1]

        voltage_data = np.array([[-0.0, 1.9415250000000697],
                                 [0.30769230769230765, 1.9415250000000697],
                                 [0.30769230769230765, 1.8750000000001705],
                                 [0.39999999999999997, 1.8750000000001705],
                                 [0.39999999999999997, 1.487874999999974],
                                 [0.5714285714285714, 1.487874999999974],
                                 [0.5714285714285714, 0.6392500000000041],
                                 [0.625, 0.6392500000000041],
                                 [0.625, 0.3461250000000007],
                                 [0.7, 0.3461250000000007],
                                 [0.7, 0.0],
                                 [1.0, 0.0]])

        points = np.delete(points, 2, axis=1)

        bare_hull = QueryConvexHull.__new__(QueryConvexHull)
        bare_hull.hull = ConvexHull(points)
        bare_hull.hull_cursor = hull_cursor
        bare_hull.match = [{'enthalpy_per_atom': -380.071/2.0}]
        bare_hull.ternary = True
        bare_hull.voltage_curve(bare_hull.hull_cursor)
        try:
            np.testing.assert_array_almost_equal(bare_hull.Vpoints[0], voltage_data)
        except:
            print('calculated: ', np.shape(bare_hull.Vpoints[0]))
            print(bare_hull.Vpoints[0])
            print('data: ', np.shape(voltage_data))
            print(voltage_data)
            np.testing.assert_array_almost_equal(bare_hull.Vpoints[0], voltage_data, decimal=8)


if __name__ == '__main__':
    unittest.main()
