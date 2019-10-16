#!/usr/bin/env python
# standard library
import unittest
from os.path import realpath

import numpy as np

# matador modules
from matador.crystal.crystal import Crystal, UnitCell
from matador.scrapers.castep_scrapers import castep2dict
from matador.scrapers.magres_scrapers import magres2dict

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'

try:
    import networkx
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
        cell = UnitCell(lattice_cart)
        self.assertAlmostEqual(cell.lattice_cart, lat_tup)
        self.assertAlmostEqual(cell.lattice_abc,((3, 3, 3), (90, 90, 90)))
        self.assertEqual(cell.volume, 27)
        self.assertAlmostEqual(cell.lengths, (3, 3, 3))
        self.assertAlmostEqual(cell.angles, (90, 90, 90))

        lattice_cart = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        cell.lattice_cart = lattice_cart
        lat_tup = tuple(tuple(vec) for vec in lattice_cart)
        self.assertEqual(cell.lattice_cart, lat_tup)
        lattice_cart = 'aadsfadsf'
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
        doc, s = castep2dict(REAL_PATH + 'data/Na3Zn4-swap-ReOs-OQMD_759599.castep')
        crystal = Crystal(doc)
        self.assertEqual(list(crystal.lattice_cart[0]), [9.0397727, 0.0081202, 0.0000000])
        self.assertEqual(crystal.num_atoms, 14)
        with self.assertRaises(AttributeError):
            crystal['positions_frac'] = [[0, 1, 2]]

    def testSites(self):
        doc, s = castep2dict(REAL_PATH + 'data/Na3Zn4-swap-ReOs-OQMD_759599.castep')
        del doc['lattice_cart']
        crystal = Crystal(doc)
        print(crystal)
        self.assertEqual(crystal[0].coords, [0.776467, 0.466319, 0.0])
        with self.assertRaises(RuntimeError):
            crystal[0].set_position([0.5, 0.6, 0.7, 0.8], 'fractional')
        with self.assertRaises(RuntimeError):
            crystal[0].set_position([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'fractional')
        self.assertEqual([atom for atom in crystal], [atom[1] for atom in enumerate(crystal)])

    def testSpg(self):
        doc, s = castep2dict(REAL_PATH + 'data/Na3Zn4-swap-ReOs-OQMD_759599.castep')
        crystal = Crystal(doc)
        print(crystal.get_space_group(symprec=0.01))
        print(crystal.get_space_group(symprec=0.001))
        self.assertEqual(crystal.get_space_group(symprec=0.0000001), 'Pm')

    def testFromMagres(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP_QE6.magres')
        crystal = Crystal(doc)
        for atom in crystal:
            print(atom, atom.magres_shift)

    @unittest.skipIf(not imported_vornet, 'Voronoi code not found in this distribution')
    def testCoordination(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP_QE6.magres')
        crystal = Crystal(doc, voronoi=True)
        for atom in crystal:
            print(atom, atom.coordination)
        print(crystal.coordination_lists)
        print(crystal.coordination_stats)

    @unittest.skipIf(not imported_vornet, 'Voronoi code not found in this distribution')
    def testVoronoi(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP_QE6.magres')
        crystal = Crystal(doc)
        print(crystal.unique_sites)

    @unittest.skipIf(not imported_networkx, 'NetworkX missing')
    def testBondLengths(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP_QE6.magres')
        crystal = Crystal(doc)
        print(crystal.bond_lengths)

    @unittest.skipIf(not imported_networkx, 'NetworkX missing')
    def testBondStats(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP_QE6.magres')
        crystal = Crystal(doc)
        print(crystal.bonding_stats)


class ElasticCrystalTest(unittest.TestCase):
    """ Test the elastic functionality of the Crystal module. """
    def testKBulkModulus(self):
        from matador.crystal.elastic import get_equation_of_state
        results = get_equation_of_state(REAL_PATH + '/data/bulk_modulus/K-bulk_modulus', plot=False)
        self.assertTrue('eos' in results)
        self.assertEqual(len(results['eos']), 3)
        self.assertAlmostEqual(results['eos'][0].bulk_modulus, 3.696117355)
        self.assertAlmostEqual(results['eos'][1].bulk_modulus, 3.699072676)
        self.assertAlmostEqual(results['eos'][2].bulk_modulus, 3.691406442)
        self.assertAlmostEqual(results['eos'][0].bulk_modulus_err, 3e-6, places=1)
        self.assertAlmostEqual(results['eos'][1].bulk_modulus_err, 2e-6, places=1)
        self.assertAlmostEqual(results['eos'][2].bulk_modulus_err, 2e-6, places=1)


if __name__ == '__main__':
    unittest.main(buffer=False, verbosity=2)
