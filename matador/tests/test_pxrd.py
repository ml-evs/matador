#!/usr/bin/env python
import unittest
from matador.fingerprints.pxrd import PXRD, PXRDFactory
from matador.scrapers.castep_scrapers import res2dict
import numpy as np
from os.path import realpath

REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'
DEBUG = True
PEAK_POS_TOL = 2


class PDFCalculatorTest(unittest.TestCase):
    """ Test PDF calculator. """

    def test_simple_pxrd(self):
        doc, success = res2dict(REAL_PATH + 'data/structures/Li.res')
        doc['pxrd-Cu'] = PXRD(doc, two_theta_resolution=0.0001, gaussian_width=0, plot=False)
        peaks = doc['pxrd-Cu'].peak_positions
        # peaks grabbed from Materials Project for Cu K-alpha
        Li_peaks = [36.01, 41.822, 60.63, 72.581, 76.369, 91.095, 102.134, 105.896, 121.916, 136.034]
        for peak in Li_peaks:
            # dodgy test for now: check every peak is present to within 2 degrees
            self.assertTrue(len(peaks[np.abs(peaks - peak) < 2]) > 2,
                            msg='Missing peak at {}'.format(peak))
        # peaks grabbed from Materials Project for Ag K-alpha
        doc['pxrd-Ag'] = PXRD(doc, two_theta_resolution=0.0001, gaussian_width=0, wavelength=0.56)
        peaks = doc['pxrd-Ag'].peak_positions
        Li_peaks = [12.912, 14.92, 21.161, 24.868, 25.992, 30.102, 32.876]
        for peak in Li_peaks:
            # dodgy test for now: check every peak is present to within 1 degree
            self.assertTrue(len(peaks[np.abs(peaks - peak) < 2]) > 2,
                            msg='Missing peak at {}'.format(peak))

        # TODO: need to check amplitudes too and factory
