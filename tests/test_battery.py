""" Test battery routines. """

import unittest
import numpy as np
from matador.battery import Electrode


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
