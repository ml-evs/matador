import unittest
import numpy as np


class VizUtilTest(unittest.TestCase):
    """Tests viz util functions."""

    def test_formula_to_colour(self):
        from matador.utils.viz_utils import formula_to_colour

        self.assertListEqual(formula_to_colour("K"), [0.4, 0.14, 0.43, 1])
        self.assertListEqual(formula_to_colour("KP"), [0.62, 0.275, 0.255, 1])
        self.assertListEqual(formula_to_colour("P"), [0.84, 0.41, 0.08, 1])
        np.testing.assert_array_almost_equal(
            formula_to_colour("KSnP"), [0.60666, 0.366666, 0.3999999, 1], decimal=3
        )
        with self.assertRaises(RuntimeError):
            formula_to_colour("KSnPSb")
