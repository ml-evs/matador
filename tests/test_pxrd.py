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
    """Check if peaks match with HKL between matador and GSAS."""
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
    """Test PXRD calculator."""

    def test_simple_pxrd(self):
        """Test Li PXRD vs GSAS."""
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

    def test_tchz_peak_shape(self):
        """Test the TCHZ pseudo-Voigt peak shape."""
        doc, s = res2dict(REAL_PATH + "data/structures/Li.res", as_model=True)
        pxrd_lor = PXRD(doc)
        pxrd = PXRD(doc, peak_shape="tchz")
        self.assertEqual(len(pxrd.pattern), len(pxrd.two_thetas))
        self.assertTrue(np.all(np.isfinite(pxrd.pattern)))
        self.assertAlmostEqual(np.max(pxrd.pattern), 1.0)
        # peak shape should not shift the strongest peak
        self.assertAlmostEqual(
            pxrd.two_thetas[np.argmax(pxrd.pattern)],
            pxrd_lor.two_thetas[np.argmax(pxrd_lor.pattern)],
            delta=0.05,
        )

        # adding a Lorentzian strain term should broaden the peaks
        pxrd_broad = PXRD(doc, peak_shape="tchz", tchz_params={"y": 0.05})
        self.assertGreater(
            np.sum(pxrd_broad.pattern > 0.01), np.sum(pxrd.pattern > 0.01)
        )

        with self.assertRaises(RuntimeError):
            PXRD(doc, peak_shape="voigt")

        with self.assertRaises(RuntimeError):
            PXRD(doc, peak_shape="tchz", tchz_params={"bad_param": 1.0})

    def test_CuP2_vs_GSAS(self):
        """Test CuP2 peak positions vs GSAS."""
        pxrd = _match_peaks(
            REAL_PATH + "data/pxrd_files/CuP2.res",
            REAL_PATH + "data/pxrd_files/CuP2_GSASII_reflections.txt",
        )
        self.assertAlmostEqual(
            pxrd.two_thetas[np.argmax(pxrd.pattern)], 30.969, places=2
        )
