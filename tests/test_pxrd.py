#!/usr/bin/env python
import unittest
from matador.fingerprints.pxrd import PXRD
from matador.scrapers.castep_scrapers import res2dict
import numpy as np
from os.path import realpath

REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"
DEBUG = True
PEAK_POS_TOL = 0.02


def _match_peaks(res_file, gsas_reflections, **kwargs):
    """ Check if peaks match with HKL between matador and GSAS. """
    doc, s = res2dict(res_file, as_model=True)

    pxrd = PXRD(doc, **kwargs)
    gsas_peaks = np.loadtxt(gsas_reflections)
    hkl_array = {}
    for peak in gsas_peaks:
        hkl_array[(int(peak[0]), int(peak[1]), int(peak[2]))] = peak[3]

    for ind, hkl in enumerate(pxrd.hkls):
        hkl_tuple = (int(hkl[0]), int(hkl[1]), int(hkl[2]))
        if hkl_tuple in hkl_array:
            assert (
                np.abs(hkl_array[hkl_tuple] - pxrd.peak_positions[ind]) < PEAK_POS_TOL
            )

    return pxrd


class PXRDCalculatorTest(unittest.TestCase):
    """ Test PXRD calculator. """

    def test_simple_pxrd(self):
        """ Test Li PXRD vs GSAS. """
        _ = _match_peaks(
            REAL_PATH + "data/structures/Li.res",
            REAL_PATH + "data/pxrd_files/Li_reflections_Cu.txt",
        )
        pxrd_args = {"wavelength": 0.559363}
        _ = _match_peaks(
            REAL_PATH + "data/structures/Li.res",
            REAL_PATH + "data/pxrd_files/Li_reflections_Ag.txt",
            **pxrd_args
        )

    def test_CuP2_vs_GSAS(self):
        """ Test CuP2 peak positions vs GSAS. """
        pxrd = _match_peaks(
            REAL_PATH + "data/pxrd_files/CuP2.res",
            REAL_PATH + "data/pxrd_files/CuP2_GSASII_reflections.txt",
        )
        self.assertAlmostEqual(
            pxrd.two_thetas[np.argmax(pxrd.pattern)], 30.969, places=2
        )
