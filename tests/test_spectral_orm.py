#!/usr/bin/env python

import unittest
import numpy as np

from matador.scrapers import castep2dict
from matador.orm.spectral import (
    VibrationalDOS,
    VibrationalDispersion,
    ElectronicDispersion,
    ElectronicDOS,
)
from .utils import REAL_PATH


class SpectralOrmTest(unittest.TestCase):
    def test_failed_make_from_castep(self):
        """ Tests that an error is raised when trying to create
        a Spectral object from a castep file that contains no
        spectral data.

        """
        fname = REAL_PATH + "data/castep_files/KP-castep17.castep"
        doc, s = castep2dict(fname, db=False)

        with self.assertRaises(RuntimeError):
            ElectronicDispersion(doc)

        with self.assertRaises(RuntimeError):
            ElectronicDispersion(**doc)

        with self.assertRaises(RuntimeError):
            ElectronicDOS(doc)

        with self.assertRaises(RuntimeError):
            ElectronicDOS(**doc)

        with self.assertRaises(RuntimeError):
            VibrationalDOS(doc)

        with self.assertRaises(RuntimeError):
            VibrationalDOS(**doc)

        with self.assertRaises(RuntimeError):
            VibrationalDispersion(doc)

        with self.assertRaises(RuntimeError):
            VibrationalDispersion(**doc)

    def test_dos_construction_from_frequencies(self):
        """ Test that we can make sensible phonon DOS's and free
        energy calculations from frequency and qpoint data, and that
        it agrees with CASTEP's calculated values.

        """

        fname = REAL_PATH + "data/castep_files/CuP-thermo-test.castep"
        doc, _ = castep2dict(fname, db=False)
        dos = VibrationalDOS(doc)
        backup_eigs = np.copy(dos.eigs)
        self.assertAlmostEqual(
            dos.zpe, doc["thermo_zero_point_energy"] / doc["num_atoms"], places=4
        )
        np.testing.assert_array_equal(
            dos.eigs, backup_eigs, err_msg="Eigenvalues were changed by some function"
        )
        self.assertAlmostEqual(
            dos.vibrational_free_energy(1000),
            doc["thermo_free_energy"][1000.0] / doc["num_atoms"],
            places=4,
        )
        np.testing.assert_array_equal(
            dos.eigs, backup_eigs, err_msg="Eigenvalues were changed by some function"
        )
        self.assertEqual(dos.compute_free_energy(0.0), dos.zpe)
