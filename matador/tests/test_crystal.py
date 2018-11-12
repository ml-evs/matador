#!/usr/bin/env python
# standard library
import unittest
from os.path import realpath
# matador modules
from matador.crystal import Crystal
from matador.scrapers.castep_scrapers import castep2dict
from matador.scrapers.magres_scrapers import magres2dict

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'

try:
    from Vornetclass import VoronoiNetwork
    imported_vornet = True
except:
    imported_vornet = False


class CrystalTest(unittest.TestCase):
    def testSites(self):
        doc, s = castep2dict(REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep')
        del doc['lattice_cart']
        crystal = Crystal(doc)
        assert [atom for atom in crystal] == [atom[1] for atom in enumerate(crystal)]

    def testSpg(self):
        doc, s = castep2dict(REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep')
        crystal = Crystal(doc)
        print(crystal.get_space_group(symprec=0.01))
        print(crystal.get_space_group(symprec=0.001))

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

    def testBondLengths(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP_QE6.magres')
        crystal = Crystal(doc)
        print(crystal.bond_lengths)


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
