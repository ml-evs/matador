import unittest
import copy
from os.path import realpath

import numpy as np

from matador.crystal.crystal import Crystal, UnitCell
from matador.crystal.crystal_site import Site
from matador.scrapers.castep_scrapers import castep2dict, res2dict
from matador.utils.cell_utils import frac2cart
from matador.scrapers.magres_scrapers import magres2dict

# grab abs path for accessing test data
REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"

try:
    import networkx  # noqa

    imported_networkx = True
except ImportError:
    imported_networkx = False

imported_vornet = False


class UnitCellTest(unittest.TestCase):
    def test_cart_init(self):
        lattice_cart = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        lat_tup = tuple(tuple(vec) for vec in lattice_cart)
        cell = UnitCell(lattice_cart)
        self.assertEqual(cell.lattice_cart, lat_tup)
        self.assertEqual(cell.lattice_abc, ((3, 3, 3), (90, 90, 90)))
        self.assertEqual(cell.volume, 27)
        self.assertEqual(cell.lengths, (3, 3, 3))
        self.assertEqual(cell.angles, (90, 90, 90))

        lattice_cart = np.asarray([[3, 0, 0], [0, 3, 0], [0, 0, 3]])
        cell_2 = UnitCell(lattice_cart)
        self.assertAlmostEqual(cell_2.lattice_cart, lat_tup)
        self.assertAlmostEqual(cell_2.lattice_abc, ((3, 3, 3), (90, 90, 90)))
        self.assertEqual(cell_2.volume, 27)
        self.assertAlmostEqual(cell_2.lengths, (3, 3, 3))
        self.assertAlmostEqual(cell_2.angles, (90, 90, 90))
        self.assertEqual(cell.lattice_cart, lat_tup)
        self.assertEqual(cell.lattice_abc, ((3, 3, 3), (90, 90, 90)))
        self.assertEqual(cell.volume, 27)
        self.assertEqual(cell.lengths, (3, 3, 3))

        lattice_cart = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        cell.lattice_cart = lattice_cart
        lat_tup = tuple(tuple(vec) for vec in lattice_cart)
        self.assertEqual(cell.lattice_cart, lat_tup)
        lattice_cart = "aadsfadsf"
        self.assertEqual(cell.lattice_cart, lat_tup)
        self.assertEqual(cell.lattice_abc, ((10, 10, 10), (90, 90, 90)))
        self.assertEqual(cell.volume, 1000)

    def test_abc_init(self):
        lattice_abc = [[2, 3, 4], [60, 60, 60]]
        lat_tup = tuple(tuple(elem) for elem in lattice_abc)
        cell = UnitCell(lattice_abc)
        self.assertAlmostEqual(cell.lattice_abc, lat_tup)
        cell.lengths = [10, 10, 10]
        self.assertEqual(cell.lattice_abc, ((10, 10, 10), (60, 60, 60)))
        cell.angles = [90, 90, 90]
        self.assertEqual(cell.lattice_abc, ((10, 10, 10), (90, 90, 90)))
        lattice_cart = ((10, 0, 0), (0, 10, 0), (0, 0, 10))
        self.assertEqual(cell.lattice_cart, lattice_cart)


class CrystalTest(unittest.TestCase):
    def test_getters_setters(self):
        doc, s = castep2dict(REAL_PATH + "data/Na3Zn4-swap-ReOs-OQMD_759599.castep")
        crystal = Crystal(doc)
        self.assertEqual(
            list(crystal.lattice_cart[0]), [9.0397727, 0.0081202, 0.0000000]
        )
        self.assertEqual(crystal.num_atoms, 14)
        with self.assertRaises(AttributeError):
            crystal["positions_frac"] = [[0, 1, 2]]

        # check we can set fields to the same value
        crystal["new_field"] = [1, 2, 3]
        crystal["new_field"] = [1, 2, 3]

        crystal["new_field_2"] = np.nan
        crystal["new_field_2"] = np.nan

        crystal["new_field_3"] = [1, 2, 4]
        with self.assertRaises(AttributeError):
            crystal["new_field_3"] = [1, 2, 5]

        crystal["new_field_4"] = [1, 2, np.nan]
        crystal["new_field_4"] = [1, 2, np.nan]

        crystal["new_field_5"] = [1, np.nan, 2]
        with self.assertRaises(AttributeError):
            crystal["new_field_5"] = [1, 2, np.nan]

        crystal["new_field_6"] = np.linspace(0, 1, 1000).tolist()
        crystal["new_field_6"] = np.array(crystal["new_field_6"], copy=True).tolist()

    def test_set_positions(self):
        doc, s = castep2dict(REAL_PATH + "data/Na3Zn4-swap-ReOs-OQMD_759599.castep")
        doc = Crystal(doc)

        copydoc = copy.deepcopy(doc)
        old_pos = np.asarray(doc.positions_frac)
        copydoc.set_positions(np.zeros_like(old_pos), fractional=True)

        np.testing.assert_array_almost_equal(
            np.asarray(copydoc.positions_frac), np.zeros_like(old_pos)
        )
        np.testing.assert_array_almost_equal(
            np.asarray(copydoc.positions_abs), np.zeros_like(old_pos)
        )

        self.assertNotAlmostEqual(doc.positions_frac[-1][0], 0.0)

    def test_convert_positions(self):
        doc = res2dict(REAL_PATH + "data/structures/Li7Sn-Fmmm.res")[0]
        crystal = res2dict(REAL_PATH + "data/structures/Li7Sn-Fmmm.res", as_model=True)[
            0
        ]

        doc["positions_abs"] = frac2cart(doc["lattice_cart"], doc["positions_frac"])

        np.testing.assert_array_almost_equal(
            doc["positions_abs"], crystal.positions_abs
        )
        for ind, site in enumerate(crystal):
            np.testing.assert_array_almost_equal(
                doc["positions_abs"][ind], site.coords_cartesian
            )

        crystal.cell.lengths = np.asarray(crystal.cell.lengths) * 10

        rescaled_pos = frac2cart(
            np.asarray(doc["lattice_cart"]) * 10, doc["positions_frac"]
        )

        for ind, site in enumerate(crystal):
            np.testing.assert_array_almost_equal(
                doc["positions_frac"][ind], site.coords
            )
            np.testing.assert_array_almost_equal(
                rescaled_pos[ind], site.coords_cartesian
            )

    def test_minimal_init(self):
        doc = Crystal(
            dict(
                lattice_abc=np.asarray([[3, 3, 3], [90, 90, 90]]),
                atom_types=["Na", "Cl"],
                positions_frac=[[0, 0, 0], [0.5, 0.5, 0.5]],
            )
        )
        self.assertEqual(doc.stoichiometry, [["Cl", 1.0], ["Na", 1.0]])
        self.assertEqual(doc.lattice_abc, ((3.0, 3.0, 3.0), (90.0, 90.0, 90.0)))
        self.assertEqual(
            doc.lattice_cart, ((3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))
        )
        self.assertEqual(len(doc.sites), 2)
        self.assertEqual(doc.num_atoms, 2)
        self.assertEqual(doc.get_concentration(), [0.5, 0.5])
        self.assertEqual(doc.concentration, [0.5, 0.5])
        self.assertEqual(
            doc.get_concentration(elements=["Na", "Cl", "Y"]), [0.5, 0.5, 0.0]
        )
        self.assertEqual(doc.concentration, [0.5, 0.5, 0.0])
        self.assertEqual(doc.positions_abs, [[0, 0, 0], [1.5, 1.5, 1.5]])
        self.assertEqual(doc.positions_frac, [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.assertEqual(doc.formula, "NaCl")
        self.assertEqual(doc.cell_volume, 27.0)
        self.assertEqual(doc.space_group, "Pm-3m")
        self.assertEqual(doc.space_group_tex, "$Pm\\bar{3}m$")

        doc = Crystal(
            dict(
                lattice_cart=((3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0)),
                atom_types=["Na", "Cl"],
                positions_abs=[[0, 0, 0], [1.5, 1.5, 1.5]],
            )
        )
        self.assertEqual(doc.lattice_abc, ((3.0, 3.0, 3.0), (90.0, 90.0, 90.0)))
        self.assertEqual(
            doc.lattice_cart, ((3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0))
        )
        self.assertEqual(doc.stoichiometry, [["Cl", 1.0], ["Na", 1.0]])
        self.assertEqual(len(doc.sites), 2)
        self.assertEqual(doc.num_atoms, 2)
        self.assertEqual(doc.get_concentration(), [0.5, 0.5])
        self.assertEqual(doc.positions_abs, [[0.0, 0.0, 0.0], [1.5, 1.5, 1.5]])
        self.assertEqual(doc.positions_frac, [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.assertEqual(doc.formula, "NaCl")
        self.assertEqual(doc.cell_volume, 27.0)
        self.assertEqual(doc.space_group, "Pm-3m")

    def testSites(self):
        doc, s = castep2dict(REAL_PATH + "data/Na3Zn4-swap-ReOs-OQMD_759599.castep")
        del doc["lattice_cart"]
        crystal = Crystal(doc)
        np.testing.assert_array_almost_equal(
            crystal[0].coords, [0.776467, 0.466319, 0.0]
        )

        with self.assertRaises(RuntimeError):
            crystal[0].set_position([0.5, 0.6, 0.7, 0.8], "fractional")
        with self.assertRaises(RuntimeError):
            crystal[0].set_position([[1, 2, 3], [4, 5, 6], [7, 8, 9]], "fractional")
        self.assertEqual(
            [atom for atom in crystal], [atom[1] for atom in enumerate(crystal)]
        )

        atom = Site(
            species="Cl",
            position=[0.2, 0.5, 0.2],
            lattice=[[10, 0, 0], [0, 10, 0], [0, 0, 10]],
        )
        atom2 = copy.deepcopy(atom)
        atom2.species = "Br"

        self.assertEqual(atom.species, "Cl")
        self.assertEqual(atom2.species, "Br")

        atom2.set_position([1.2, -0.5, 0.2], "fractional")
        np.testing.assert_array_almost_equal(
            atom2.displacement_between_sites(atom), [0.0, 0.0, 0.0], decimal=10
        )
        self.assertAlmostEqual(atom2.distance_between_sites(atom), 0.0, places=10)
        atom2.set_position([1.3, -0.5, 0.2], "fractional")
        np.testing.assert_array_almost_equal(
            atom2.displacement_between_sites(atom), [1.0, 0.0, 0.0], decimal=10
        )
        self.assertAlmostEqual(atom2.distance_between_sites(atom), 1.0, places=10)
        atom2.set_position([1.3, -0.5, 0.3], "fractional")
        np.testing.assert_array_almost_equal(
            atom2.displacement_between_sites(atom), [1.0, 0.0, 1.0], decimal=10
        )
        self.assertAlmostEqual(
            atom2.distance_between_sites(atom), np.sqrt(2), places=10
        )

    def testSpg(self):
        doc, s = castep2dict(REAL_PATH + "data/Na3Zn4-swap-ReOs-OQMD_759599.castep")
        crystal = Crystal(doc)
        print(crystal.get_space_group(symprec=0.01))
        print(crystal.get_space_group(symprec=0.001))
        self.assertEqual(crystal.get_space_group(symprec=0.0000001), "Pm")

    def testFromMagres(self):
        doc, s = magres2dict(REAL_PATH + "data/magres_files/NaP_QE6.magres")
        crystal = Crystal(doc)
        for atom in crystal:
            print(
                atom, atom["chemical_shielding_iso"], atom["chemical_shift_asymmetry"]
            )

    @unittest.skipIf(not imported_vornet, "Voronoi code not found in this distribution")
    def testCoordination(self):
        doc, s = magres2dict(REAL_PATH + "data/magres_files/NaP_QE6.magres")
        crystal = Crystal(doc, voronoi=True)
        for atom in crystal:
            print(atom, atom.coordination)
        print(crystal.coordination_lists)
        print(crystal.coordination_stats)

    @unittest.skipIf(not imported_vornet, "Voronoi code not found in this distribution")
    def testVoronoi(self):
        doc, s = magres2dict(REAL_PATH + "data/magres_files/NaP_QE6.magres")
        crystal = Crystal(doc)
        print(crystal.unique_sites)

    @unittest.skipIf(not imported_networkx, "NetworkX missing")
    def testBondLengths(self):
        doc, s = magres2dict(REAL_PATH + "data/magres_files/NaP_QE6.magres")
        crystal = Crystal(doc)
        print(crystal.bond_lengths)

    @unittest.skipIf(not imported_networkx, "NetworkX missing")
    def testBondStats(self):
        doc, s = magres2dict(REAL_PATH + "data/magres_files/NaP_QE6.magres")
        crystal = Crystal(doc)
        print(crystal.bonding_stats)


class ElasticCrystalTest(unittest.TestCase):
    """Test the elastic functionality of the Crystal module."""

    def testKBulkModulus(self):
        from matador.crystal.elastic import get_equation_of_state

        results = get_equation_of_state(
            REAL_PATH + "/data/bulk_modulus/K-bulk_modulus", plot=False
        )
        self.assertTrue("eos" in results)
        self.assertEqual(len(results["eos"]), 3)
        self.assertAlmostEqual(results["eos"][0].bulk_modulus, 3.696117355)
        self.assertAlmostEqual(results["eos"][1].bulk_modulus, 3.699072676)
        self.assertAlmostEqual(results["eos"][2].bulk_modulus, 3.691406442)
        self.assertAlmostEqual(results["eos"][0].bulk_modulus_err, 3e-6, places=1)
        self.assertAlmostEqual(results["eos"][1].bulk_modulus_err, 2e-6, places=1)
        self.assertAlmostEqual(results["eos"][2].bulk_modulus_err, 2e-6, places=1)


if __name__ == "__main__":
    unittest.main(buffer=False, verbosity=2)
