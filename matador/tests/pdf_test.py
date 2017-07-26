#!/usr/bin/env python
import unittest
from matador.similarity.pdf_similarity import PDF, PDFOverlap
from matador.scrapers.castep_scrapers import res2dict
from matador.utils.cell_utils import abc2cart, cart2volume
import numpy as np
from os.path import realpath

REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class PDFCalculatorTest(unittest.TestCase):
    """ Test PDF calculator. """
    def testIdealGasPDF(self, retry=0):
        # create fake matador doc
        doc = dict()
        max_retries = 3
        self.assertLess(retry, max_retries, msg='After {} attempts, PDF still failed.'.format(retry))
        num_atoms = 400
        box_size = 40
        num_samples = 10
        rmax = 15
        dr = 0.01
        i = 0
        doc['atom_types'] = num_atoms * ['C']
        doc['lattice_cart'] = np.asarray([[box_size, 0, 0], [0, box_size, 0], [0, 0, box_size]])
        doc['cell_volume'] = box_size**3
        doc['text_id'] = ['ideal', 'gas']
        doc['Gr_smear'] = np.zeros_like(np.arange(0, rmax+dr, dr))
        doc['Gr_hist'] = np.zeros_like(doc['Gr_smear'])
        while i < num_samples:
            doc['positions_frac'] = np.random.rand(num_atoms, 3)
            doc['pdf'] = PDF(doc, num_images=1, dr=dr, rmax=rmax, lazy=True, style='histogram')
            doc['pdf']._calc_pdf()
            doc['pdf_smear'] = PDF(doc, num_images=1, gaussian_width=0.01, dr=dr, rmax=rmax, lazy=True, style='smear')
            doc['pdf_smear']._calc_pdf()
            doc['Gr_smear'] += doc['pdf_smear'].Gr / num_samples
            doc['Gr_hist'] += doc['pdf'].Gr / num_samples
            i += 1
        try:
            self.assertAlmostEqual(np.mean(doc['Gr_smear']), np.mean(doc['Gr_hist']), places=1)
            self.assertAlmostEqual(np.mean(doc['Gr_smear']), 1.0, places=1)
            self.assertAlmostEqual(np.mean(doc['Gr_hist']), 1.0, places=1)
        except:
            self.testIdealGasPDF(retry=retry+1)

    def testPDFAutoImageNumber(self):
        doc, success = res2dict(REAL_PATH + 'data/LiPZn-r57des.res')
        doc['lattice_cart'] = abc2cart(doc['lattice_abc'])
        doc['text_id'] = ['pdf', 'test']
        doc['pdf_num_images'] = PDF(doc, num_images=5, **{'debug': True})
        doc['pdf_auto_images'] = PDF(doc, num_images='auto', **{'debug': True})
        np.testing.assert_array_almost_equal(doc['pdf_num_images'].Gr, doc['pdf_auto_images'].Gr)

    def testOverlapPDFSameStructure(self):
        doc, success = res2dict(REAL_PATH + 'data/LiPZn-r57des.res')
        doc['lattice_cart'] = abc2cart(doc['lattice_abc'])
        doc['text_id'] = ['pdf', 'test']
        doc['pdf_smear'] = PDF(doc, num_images=3, dr=0.001, gaussian_width=0.01, style='smear', debug=True, low_mem=True)
        overlap = PDFOverlap(doc['pdf_smear'], doc['pdf_smear'])
        self.assertEqual(overlap.similarity_distance, 0.0)

    def testOverlapPDFHistVSSmear(self):
        doc, success = res2dict(REAL_PATH + 'data/LiPZn-r57des.res')
        doc['lattice_cart'] = abc2cart(doc['lattice_abc'])
        doc['text_id'] = ['pdf', 'test']
        doc['pdf_smear'] = PDF(doc, num_images=3, dr=0.001, gaussian_width=0.01, style='smear', low_mem=True)
        doc['pdf_hist'] = PDF(doc, num_images=3, dr=0.1, style='histogram')
        overlap = PDFOverlap(doc['pdf_smear'], doc['pdf_hist'])
        self.assertLessEqual(overlap.similarity_distance, 0.02)
        self.assertGreater(overlap.similarity_distance, 0.0)

    def testPDFPrimitiveVsSupercell(self):
        test_doc, success = res2dict(REAL_PATH + 'data/KP_primitive.res', db=False)
        test_doc['text_id'] = ['primitive', 'cell']
        test_doc['lattice_cart'] = abc2cart(test_doc['lattice_abc'])
        test_doc['cell_volume'] = cart2volume(test_doc['lattice_cart'])
        supercell_doc, success = res2dict(REAL_PATH + 'data/KP_supercell.res', db=False)
        supercell_doc['text_id'] = ['supercell', 'cell']
        supercell_doc['lattice_cart'] = abc2cart(supercell_doc['lattice_abc'])
        supercell_doc['cell_volume'] = cart2volume(supercell_doc['lattice_cart'])
        test_doc['pdf'] = PDF(test_doc, dr=0.01, low_mem=True, rmax=10, num_images='auto', debug=True)
        supercell_doc['pdf'] = PDF(supercell_doc, dr=0.01, low_mem=True, rmax=10, num_images='auto', debug=True)
        overlap = PDFOverlap(test_doc['pdf'], supercell_doc['pdf'])
        self.assertLessEqual(overlap.similarity_distance, 1e-3)
        self.assertGreater(overlap.similarity_distance, 0.0)


if __name__ == '__main__':
    unittest.main()
