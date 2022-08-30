#!/usr/bin/env python

import unittest

try:
    import ase  # noqa
    import ase.build

    ASE_IMPORTED = True
except ImportError:
    ASE_IMPORTED = False

import numpy as np

from matador.utils.ase_utils import ase2dict
from .utils import REAL_PATH


@unittest.skipIf(not ASE_IMPORTED, "ASE not found")
class AseUtilTest(unittest.TestCase):
    """Tests cursor util functions."""

    def setUp(self):
        # construct ASE atoms object manually
        self.ase_atoms = ase.build.bulk(
            "SiO", crystalstructure="rocksalt", a=12, cubic=True
        )
        self.ase_atoms.info = {"test_info": "dictionary"}

    def test_ase2dict(self):
        doc = ase2dict(self.ase_atoms)
        self.assertEqual(
            doc["atom_types"], sorted(["Si", "O", "Si", "O", "Si", "O", "Si", "O"])
        )
        self.assertListEqual(doc["lattice_cart"][0], [12, 0, 0])
        self.assertListEqual(doc["lattice_cart"][1], [0, 12, 0])
        self.assertListEqual(doc["lattice_cart"][2], [0, 0, 12])
        self.assertListEqual(doc["lattice_abc"][0], [12, 12, 12])
        self.assertListEqual(doc["lattice_abc"][1], [90, 90, 90])
        self.assertAlmostEqual(doc["cell_volume"], 12**3)
        self.assertEqual(len(doc["positions_frac"]), len(doc["atom_types"]))
        self.assertListEqual(doc["positions_frac"][0], [0.5, 0.0, 0.0])
        self.assertListEqual(doc["positions_frac"][1], [0.5, 0.5, 0.5])
        self.assertListEqual(doc["positions_frac"][3], [0.0, 0.5, 0.0])
        self.assertEqual(doc["num_atoms"], 8)
        self.assertEqual(doc["num_fu"], 4)
        self.assertListEqual(doc["stoichiometry"][0], ["O", 1.0])
        self.assertListEqual(doc["stoichiometry"][1], ["Si", 1.0])
        self.assertListEqual(sorted(list(doc["elems"])), ["O", "Si"])
        self.assertEqual(doc["space_group"], "Fm-3m")
        self.assertDictEqual(doc["ase_info"], self.ase_atoms.info)

    def test_doc2ase(self):
        from matador.scrapers import castep2dict

        doc, s = castep2dict(
            REAL_PATH + "data/castep_files/KP-castep17.castep", as_model=True
        )

        ase_atoms = doc.ase_atoms
        np.testing.assert_array_almost_equal(
            ase_atoms.cell.array, doc.lattice_cart, decimal=16
        )

        np.testing.assert_array_almost_equal(
            doc.positions_frac, ase_atoms.get_scaled_positions(), decimal=12
        )

        for ind, site in enumerate(doc):
            np.testing.assert_array_almost_equal(
                site.coords_cartesian, ase_atoms[ind].position, decimal=12
            )
            self.assertEqual(site.species, ase_atoms[ind].symbol)

        self.assertDictEqual(ase_atoms.info["matador"], doc._data)
