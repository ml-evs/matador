#!/usr/bin/env python

import unittest

try:
    import ase  # noqa
    import ase.build

    ASE_IMPORTED = True
except ImportError:
    ASE_IMPORTED = False

from matador.utils.ase_utils import ase2dict


@unittest.skipIf(not ASE_IMPORTED, "ASE not found")
class AseUtilTest(unittest.TestCase):
    """ Tests cursor util functions. """

    def setUp(self):
        # construct ASE atoms object manually
        self.ase_atoms = ase.build.bulk(
            "SiO", crystalstructure="rocksalt", a=12, cubic=True
        )

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
        self.assertAlmostEqual(doc["cell_volume"], 12 ** 3)
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

    def test_doc2ase(self):
        pass


if __name__ == "__main__":
    unittest.main()
