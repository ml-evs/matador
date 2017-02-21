#!/usr/bin/env python
import unittest
from matador.utils.cell_utils import abc2cart, cart2abc, cart2volume
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
            except(AssertionError):
                print('cart:', test_doc['lattice_cart'], abc2cart(test_doc['lattice_abc']))
                print('abc:', test_doc['lattice_abc'], cart2abc(test_doc['lattice_cart']))
                print('volume:', test_doc['cell_volume'], cart2volume(test_doc['lattice_cart']))
                raise AssertionError


if __name__ == '__main__':
    unittest.main()
