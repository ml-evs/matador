#!/usr/bin/env python
import unittest
from matador.similarity.pdf_similarity import PDF
import numpy as np


class PDFCalculatorTest(unittest.TestCase):
    """ Test PDF functionality on ideal gas. """
    def testIdealGasPDF(self):
        # create fake matador doc
        doc = dict()
        num_atoms = 200
        box_size = 40
        num_samples = 10
        rmax = 15
        dr = 0.01
        i = 0
        doc['atom_types'] = num_atoms * ['C']
        doc['lattice_cart'] = np.asarray([[box_size, 0, 0], [0, box_size, 0], [0, 0, box_size]])
        doc['cell_volume'] = box_size**3
        doc['text_id'] = ['ideal', 'gas']
        doc['Gr_smear'] = np.zeros_like(np.arange(0, rmax, dr))
        doc['Gr_hist'] = np.zeros_like(doc['Gr_smear'])
        while i < num_samples:
            doc['positions_frac'] = np.random.rand(num_atoms, 3)
            doc['pdf'] = PDF(doc, num_images=1, dr=dr, rmax=rmax, lazy=True)
            doc['pdf']._calc_pdf(style='histogram')
            doc['pdf_smear'] = PDF(doc, num_images=1, gaussian_width=0.01, dr=dr, rmax=rmax, lazy=True)
            doc['pdf_smear']._calc_pdf(style='smear')
            doc['Gr_smear'] += doc['pdf_smear'].Gr / num_samples
            doc['Gr_hist'] += doc['pdf'].Gr / num_samples
            i += 1
        print('Smear: {}, Hist: {}'.format(np.mean(doc['Gr_smear']), np.mean(doc['Gr_hist'])))
        self.assertAlmostEqual(np.mean(doc['Gr_smear']), np.mean(doc['Gr_hist']), places=1)
        self.assertAlmostEqual(np.mean(doc['Gr_smear']), 1.0, places=1)
        self.assertAlmostEqual(np.mean(doc['Gr_hist']), 1.0, places=1)


if __name__ == '__main__':
    unittest.main()
