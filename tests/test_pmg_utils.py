import unittest

try:
    import ase  # noqa

    ASE_IMPORTED = True
except ImportError:
    ASE_IMPORTED = False

try:
    import pymatgen  # noqa

    PMG_IMPORTED = True
except ImportError:
    PMG_IMPORTED = False

import numpy as np


from .utils import REAL_PATH


@unittest.skipIf((not (PMG_IMPORTED and ASE_IMPORTED)), "Unable to import ASE/pymtagen")
class PMGUtilTest(unittest.TestCase):
    """Tests cursor util functions."""

    def test_doc2pmg(self):
        from matador.scrapers import castep2dict

        doc, s = castep2dict(
            REAL_PATH + "data/castep_files/KP-castep17.castep", as_model=True
        )

        pmg_structure = doc.pmg_structure
        self.assertEqual(pmg_structure.formula, "K7 P2")
        np.testing.assert_array_almost_equal(
            pmg_structure.lattice.matrix, doc.lattice_cart, decimal=16
        )
        for ind, site in enumerate(doc):
            np.testing.assert_array_almost_equal(
                site.coords, pmg_structure[ind].frac_coords, decimal=16
            )
            self.assertEqual(site.species, pmg_structure[ind].species_string)

        self.assertDictEqual(pmg_structure.info["matador"], doc._data)

    def test_pmg2dict(self):
        import ase.build
        from pymatgen.io.ase import AseAtomsAdaptor
        from matador.utils.pmg_utils import pmg2dict

        # construct ASE atoms object manually
        ase_atoms = ase.build.bulk("SiO", crystalstructure="rocksalt", a=12, cubic=True)

        pmg_structure = AseAtomsAdaptor.get_structure(ase_atoms)
        doc = pmg2dict(pmg_structure, as_model=True)

        self.assertEqual(
            doc["atom_types"], sorted(["Si", "O", "Si", "O", "Si", "O", "Si", "O"])
        )

        np.testing.assert_array_almost_equal(
            pmg_structure.lattice.matrix, doc.lattice_cart, decimal=16
        )
        self.assertAlmostEqual(doc["cell_volume"], 12**3)
        self.assertEqual(len(doc["positions_frac"]), len(doc["atom_types"]))
        self.assertEqual(doc["num_atoms"], 8)
        self.assertEqual(doc["num_fu"], 4)
        self.assertListEqual(doc["stoichiometry"][0], ["O", 1.0])
        self.assertListEqual(doc["stoichiometry"][1], ["Si", 1.0])
        self.assertListEqual(sorted(list(doc["elems"])), ["O", "Si"])
        self.assertEqual(doc["space_group"], "Fm-3m")
        self.assertEqual(pmg_structure.formula, "Si4 O4")
