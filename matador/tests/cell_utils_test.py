#!/usr/bin/env python
import unittest
from matador.utils.cell_utils import abc2cart, cart2abc, cart2volume, create_simple_supercell, doc2spg
from matador.scrapers.castep_scrapers import castep2dict, res2dict
from matador.export import doc2res
from matador.similarity.pdf_similarity import PDF, PDFOverlap
from functools import reduce
from spglib import find_primitive
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
        num_tests = 3
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', castep_fname, '- please check installation.')
        if not failed_open:
            f.close()
            test_doc, s = castep2dict(castep_fname, db=True, verbosity=5)
            _iter = 0
            while _iter < num_tests:
                extension = np.random.randint(low=1, high=5, size=(3)).tolist()
                if extension == [1, 1, 1]:
                    extension[np.random.randint(low=0, high=2)] += 1
                num_images = reduce(lambda x, y: x*y, extension)

                supercell = create_simple_supercell(test_doc, tuple(extension))
                self.assertEqual(supercell['num_atoms'], num_images*test_doc['num_atoms'])
                self.assertAlmostEqual(supercell['cell_volume'], num_images*test_doc['cell_volume'], places=3)
                self.assertEqual(len(supercell['positions_frac']), num_images*len(test_doc['positions_frac']))
                for i in range(3):
                    np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][i]), extension[i]*np.asarray(test_doc['lattice_cart'][i]))
                self.assertLess(pdf_sim_dist(test_doc, supercell), 1e-3)
                _iter += 1

            # test error for 1x1x1
            try:
                supercell = create_simple_supercell(test_doc, (1, 1, 1))
                error = False
            except:
                error = True
            self.assertTrue(error)

        res_fname = REAL_PATH + 'data/parent2.res'
        failed_open = False
        try:
            f = open(res_fname, 'r')
        except:
            failed_open = True
        if not failed_open:
            f.close()
            test_doc, s = res2dict(res_fname, db=False, verbosity=0)
            while _iter < num_tests:
                extension = np.random.randint(low=1, high=5, size=(3, 1)).tolist()
                num_images = reduce(map(lambda x, y: x*y, extension))

                supercell = create_simple_supercell(test_doc, tuple(extension))
                self.assertEqual(supercell['num_atoms'], num_images*test_doc['num_atoms'])
                self.assertAlmostEqual(supercell['cell_volume'], num_images*test_doc['cell_volume'], places=3)
                self.assertEqual(len(supercell['positions_frac']), num_images*len(test_doc['positions_frac']))
                for i in range(3):
                    np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][i]), extension[i]*np.asarray(test_doc['lattice_cart'][i]))
                self.assertLess(pdf_sim_dist(test_doc, supercell), 1e-3)
                _iter += 1

            # test error for 1x1x1
            try:
                supercell = create_simple_supercell(test_doc, (1, 1, 1))
                error = False
            except:
                error = True
            self.assertTrue(error)

            # spg_test_doc = doc2spg(test_doc)
            # spg_supercell = doc2spg(supercell)

            # if spg_supercell is None:
                # raise RuntimeError('Unable to convert supercell to spg.')

            # doc2res(supercell, '/home/matthew/super_test.res', info=False)
            # print(supercell['positions_frac'])
            # print(spg_supercell)

            # primitive_test_doc = find_primitive(spg_test_doc)
            # primitive_supercell = find_primitive(spg_supercell)

            # if primitive_supercell is None:
                # raise RuntimeError('Unable to convert supercell to spg.')

            # print('og')
            # print(primitive_test_doc)
            # print('sup')
            # print(primitive_supercell)


def pdf_sim_dist(doc_test, doc_supercell):
    doc_test['text_id'] = ['test', 'cell']
    doc_supercell['text_id'] = ['super', 'cell']
    pdf_test = PDF(doc_test, low_mem=True)
    pdf_supercell = PDF(doc_supercell, low_mem=True)
    overlap = PDFOverlap(pdf_test, pdf_supercell)
    return overlap.similarity_distance

if __name__ == '__main__':
    unittest.main()
