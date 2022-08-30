# coding: utf-8
# Distributed under the terms of the MIT license.

import unittest
import numpy as np
import pytest


class TestMagresReferencer(unittest.TestCase):
    def test_fit_only(self):
        from matador.magres import MagresReferencer
        from matador.crystal import Crystal

        SiO2 = Crystal(
            dict(
                lattice_cart=[[10, 0, 0], [0, 10, 0], [0, 0, 10]],
                positions_frac=[[0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [0.2, 0.2, 0.2]],
                atom_types=["Si", "O", "O"],
                chemical_shielding_isos=[-100, -200, -300],
            )
        )

        expt_shifts = {"SiO2": {"O": [100, 200]}}

        ref = MagresReferencer({"SiO2": SiO2}, expt_shifts, species="O")

        self.assertAlmostEqual(ref.fit_gradient, -1)
        self.assertAlmostEqual(ref.fit_intercept, -100)
        self.assertAlmostEqual(ref.fit_rsquared, 1)

    def test_assign_shifts(self):
        from matador.magres import MagresReferencer
        from matador.crystal import Crystal

        SiO2 = Crystal(
            dict(
                lattice_cart=[[10, 0, 0], [0, 10, 0], [0, 0, 10]],
                positions_frac=[[0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [0.2, 0.2, 0.2]],
                atom_types=["Si", "O", "O"],
                chemical_shielding_isos=[-100, -200, -300],
            )
        )

        LiCoO3_supercell = Crystal(
            dict(
                lattice_cart=[[10, 0, 0], [0, 10, 0], [0, 0, 10]],
                positions_frac=[
                    [0.0, 0.0, 0.0],
                    [0.1, 0.1, 0.1],
                    [0.2, 0.2, 0.2],
                    [0.3, 0.3, 0.3],
                    [0.4, 0.4, 0.4],
                    [0.5, 0.5, 0.5],
                    [0.6, 0.6, 0.6],
                    [0.7, 0.7, 0.7],
                    [0.8, 0.8, 0.8],
                    [0.9, 0.9, 0.9],
                ],
                atom_types=["Li", "Li", "Co", "Co", "O", "O", "O", "O", "O", "O"],
                chemical_shielding_isos=[
                    -100,
                    -120,
                    -321,
                    -992,
                    -1500,
                    -1600,
                    -1700,
                    -1800,
                    -1900,
                    -2000,
                ],
            )
        )

        LiCoO3_theory = Crystal(
            dict(
                lattice_cart=[[10, 0, 0], [0, 10, 0], [0, 0, 10]],
                positions_frac=[
                    [0.0, 0.0, 0.0],
                    [0.1, 0.1, 0.1],
                    [0.2, 0.2, 0.2],
                    [0.3, 0.3, 0.3],
                    [0.4, 0.4, 0.4],
                    [0.5, 0.5, 0.5],
                    [0.6, 0.6, 0.6],
                    [0.7, 0.7, 0.7],
                    [0.8, 0.8, 0.8],
                    [0.9, 0.9, 0.9],
                ],
                atom_types=["Li", "Li", "Co", "Co", "O", "O", "O", "O", "O", "O"],
                chemical_shielding_isos=[-100, -120, -321, -992, 1, 2, 3, 4, 5, 6],
            )
        )

        expt_shifts = {
            "SiO2": {"O": [100, 200]},
            "LiCoO3": {"O": [1400, 1500, 1600, 1700, 1800, 1900]},
        }

        ref = MagresReferencer(
            {"SiO2": SiO2, "LiCoO3": LiCoO3_supercell},
            expt_shifts,
            species="O",
            structures=[LiCoO3_theory],
        )

        self.assertAlmostEqual(ref.fit_gradient, -1)
        self.assertAlmostEqual(ref.fit_intercept, -100)
        self.assertAlmostEqual(ref.fit_rsquared, 1)

        LiCoO3_theory = ref.structures[0]
        shifts = [
            site["chemical_shift_iso"] for site in LiCoO3_theory if site.species == "O"
        ]

        np.testing.assert_array_almost_equal(
            shifts, [-101, -102, -103, -104, -105, -106]
        )

        np.testing.assert_array_almost_equal(
            LiCoO3_theory["chemical_shift_isos"][4:],
            [-101, -102, -103, -104, -105, -106],
        )

        with pytest.raises(AttributeError):
            LiCoO3_theory[-1]["chemical_shift_iso"] = 123

        with pytest.raises(ValueError):
            LiCoO3_theory[0][-1]["chemical_shift_iso"] = 123
