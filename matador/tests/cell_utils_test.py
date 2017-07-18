#!/usr/bin/env python
import unittest
from matador.utils.cell_utils import abc2cart, cart2abc, cart2volume, create_simple_supercell, standardize_doc_cell
from matador.scrapers.castep_scrapers import castep2dict
import numpy as np
from os.path import realpath

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class CellUtilTest(unittest.TestCase):
    """ Tests cell util functions. """
    def testCart2AbcConversions(self):
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', castep_fname, '- please check installation.')
        if not failed_open:
            f.close()
            test_doc, s = castep2dict(castep_fname, db=True, verbosity=5)
            try:
                self.assertTrue(np.allclose(test_doc['lattice_abc'], cart2abc(test_doc['lattice_cart'])),
                                msg='Conversion cart2abc failed.')
                self.assertTrue(np.allclose(cart2abc(test_doc['lattice_cart']), cart2abc(abc2cart(test_doc['lattice_abc']))),
                                msg='Conversion abc2cart failed.')
                self.assertAlmostEqual(test_doc['cell_volume'], cart2volume(test_doc['lattice_cart']),
                                       msg='Failed to calculate volume from lattice vectors.', places=5)
                self.assertIsInstance(test_doc['lattice_abc'], list, msg='Failed abc numpy cast to list')
                self.assertIsInstance(test_doc['lattice_cart'], list, msg='Failed cartesian numpy cast to list')
            except(AssertionError):
                print('cart:', test_doc['lattice_cart'], abc2cart(test_doc['lattice_abc']))
                print('abc:', test_doc['lattice_abc'], cart2abc(test_doc['lattice_cart']))
                print('volume:', test_doc['cell_volume'], cart2volume(test_doc['lattice_cart']))
                raise AssertionError

    def testSupercellCreator(self):
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', castep_fname, '- please check installation.')
        if not failed_open:
            f.close()
            test_doc, s = castep2dict(castep_fname, db=True, verbosity=5)
            # test simple 2x2x2
            supercell = create_simple_supercell(test_doc, (2, 2, 2))
            self.assertEqual(supercell['num_atoms'], 8*test_doc['num_atoms'])
            self.assertEqual(len(supercell['positions_frac']), 8*len(test_doc['positions_frac']))
            np.testing.assert_array_equal(np.asarray(supercell['lattice_cart']), 2*np.asarray(test_doc['lattice_cart']))

            # test error for 1x1x1
            try:
                supercell = create_simple_supercell(test_doc, (1, 1, 1))
                error = False
            except:
                error = True
            self.assertTrue(error)

            # test non-diagonal
            supercell = create_simple_supercell(test_doc, (2, 1, 2))
            self.assertEqual(supercell['num_atoms'], 4*test_doc['num_atoms'])
            self.assertEqual(len(supercell['positions_frac']), 4*len(test_doc['positions_frac']))
            np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][0]), 2*np.asarray(test_doc['lattice_cart'][0]))
            np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][1]), np.asarray(test_doc['lattice_cart'][1]))
            np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][2]), 2*np.asarray(test_doc['lattice_cart'][2]))

            # test cell standardization
            supercell = create_simple_supercell(test_doc, (2, 1, 2), standardize=True)
            test_doc = standardize_doc_cell(test_doc)
            self.assertEqual(supercell['num_atoms'], 4*test_doc['num_atoms'])
            self.assertEqual(len(supercell['positions_frac']), 4*len(test_doc['positions_frac']))
            np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][0]), 2*np.asarray(test_doc['lattice_cart'][0]))
            np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][1]), np.asarray(test_doc['lattice_cart'][1]))
            np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][2]), 2*np.asarray(test_doc['lattice_cart'][2]))


if __name__ == '__main__':
    unittest.main()
