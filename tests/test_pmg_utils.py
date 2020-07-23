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
    """ Tests cursor util functions. """

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
