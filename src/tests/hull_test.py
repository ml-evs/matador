#!/usr/bin/env python

import unittest
import json
import numpy as np
from hull import QueryConvexHull
from traceback import print_exc

class VoltageTest(unittest.TestCase):
    """ Test voltage_curve(). """

    def test(self):
        print('Testing voltage_curve()...')
        match, hull_cursor = [], []
        test_x = np.loadtxt('data/x.dat')
        test_Q = np.loadtxt('data/Q.dat')
        test_V = np.loadtxt('data/V.dat')
        for i in range(5):
            with open('data/hull_data' + str(i) +'.json') as f:
                hull_cursor.append(json.load(f))
        for i in range(2):
            with open('data/mu' + str(i) +'.json') as f:
                match.append(json.load(f))
        with open('data/elements.json') as f:
            elements = json.load(f)
        bare_hull = QueryConvexHull.__new__(QueryConvexHull)
        bare_hull.cursor = list(hull_cursor)
        bare_hull.elements = list(elements)
        bare_hull.hull_cursor = list(hull_cursor)
        bare_hull.match = list(match)
        bare_hull.voltage_curve()
        np.testing.assert_array_equal(bare_hull.voltages, test_V, verbose=True)
        np.testing.assert_array_equal(bare_hull.x, test_x)
        np.testing.assert_array_equal(bare_hull.Q, test_Q)

if __name__ == '__main__':
    unittest.main()
