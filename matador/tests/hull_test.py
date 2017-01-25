#!/usr/bin/env python

import unittest
import json
import numpy as np
from matador.hull import QueryConvexHull
from os.path import realpath
import sys
import os

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class VoltageTest(unittest.TestCase):
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
        bare_hull.elements = list(elements)
        bare_hull.hull_cursor = list(hull_cursor)
        bare_hull.match = list(match)
        with open(os.devnull, 'w') as sys.stdout:
            bare_hull.voltage_curve()
        sys.stdout = sys.__stdout__
        np.testing.assert_array_equal(bare_hull.voltages, test_V, verbose=True)
        np.testing.assert_array_equal(bare_hull.x, test_x)
        np.testing.assert_array_equal(bare_hull.Q, test_Q)

if __name__ == '__main__':
    unittest.main()
