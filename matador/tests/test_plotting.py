#!/usr/bin/env python

""" These tests only check whether plots are created,
not that they look correct!

"""

import unittest
import os
import sys
from glob import glob
import matador.cli.dispersion
from matador.scrapers import res2dict
from matador.hull import QueryConvexHull

REAL_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
ROOT_DIR = os.getcwd()


class SpectralPlotTests(unittest.TestCase):
    """ Test Dispersion script. """
    def test_pdis_plot(self):
        """ Test combined spectral plots. """
        os.chdir(REAL_PATH + '/data/dispersion')
        expected_file = 'K3P-OQMD_4786-CollCode25550_spectral.png'
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = ['dispersion', 'K3P-OQMD_4786-CollCode25550', '--png',
                    '-scale', '10', '-interp', '2', '-pw', '-5', '5', '--gap',
                    '--preserve_kspace_distance', '--figsize', '10', '10']
        error = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            print(exc)
            error = True
        file_exists = os.path.isfile(expected_file)
        os.remove(expected_file)
        os.chdir(ROOT_DIR)
        self.assertFalse(error)
        self.assertTrue(file_exists)

    def test_dos_only(self):
        """ Test combined spectral plots. """
        os.chdir(REAL_PATH + '/data/dispersion')
        expected_file = 'K3P-OQMD_4786-CollCode25550_spectral.png'
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = ['dispersion', 'K3P-OQMD_4786-CollCode25550', '--png',
                    '--dos_only',
                    '--figsize', '10', '10']
        error = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            print(exc)
            error = True
        file_exists = os.path.isfile(expected_file)
        os.remove(expected_file)
        os.chdir(ROOT_DIR)
        self.assertFalse(error)
        self.assertTrue(file_exists)

    def test_multiseed(self):
        """ Test plotting two seed bandstructures on top of each other. """
        os.chdir(REAL_PATH + '/data/bands_files')
        expected_file = 'KPSn_spectral.png'
        sys.argv = ['dispersion', 'KPSn', 'KPSn_2.bands',
                    '--dos_only', '--cmap', 'viridis',
                    '--png',
                    '--labels', 'PBE, LDA',
                    '--figsize', '10', '10']
        error = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            print(exc)
            error = True
        file_exists = os.path.isfile(expected_file)
        os.remove(expected_file)
        os.chdir(ROOT_DIR)
        self.assertFalse(error)
        self.assertTrue(file_exists)

    def test_x11_no_fail(self):
        """ Test combined spectral plots. """
        os.chdir(REAL_PATH + '/data/dispersion')
        sys.argv = ['dispersion', 'K3P-OQMD_4786-CollCode25550',
                    '--dos_only', '--cmap', 'viridis',
                    '--figsize', '10', '10']
        error = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            print(exc)
            error = True
        os.chdir(ROOT_DIR)
        self.assertFalse(error)

    def test_phonon(self):
        """ Test phonon dispersion plot. """
        os.chdir(REAL_PATH + '/data')
        expected_file = 'K8SnP4_spectral.png'
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = ['dispersion', 'K8SnP4', '--png', '-ph', '--band_colour', 'k',
                    '--figsize', '10', '10']
        error = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            print(exc)
            error = True
        file_exists = os.path.isfile(expected_file)
        os.remove(expected_file)
        os.chdir(ROOT_DIR)
        self.assertFalse(error)
        self.assertTrue(file_exists)


class HullPlotTests(unittest.TestCase):
    """ Tests for plotting phase diagrams. """
    def test_binary_hull_plot(self):
        """ Test plotting binary hull. """
        expected_files = ['KP_hull.png', 'KP_voltage.png', 'KP_volume.png']
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
        res_list = glob(REAL_PATH + 'data/hull-KP-KSnP_pub/*.res')
        self.assertEqual(len(res_list), 295, 'Could not find test res files, please check installation...')
        cursor = [res2dict(res)[0] for res in res_list]
        QueryConvexHull(cursor=cursor, elements=['K', 'P'], no_plot=False, png=True, quiet=False, subcmd='voltage',
                        labels=True, label_cutoff=0.05, colour_by_source=True, hull_cutoff=0.1, volume=True)
        for expected_file in expected_files:
            self.assertTrue(os.path.isfile(expected_file))
        for expected_file in expected_files:
            os.remove(expected_file)

    def test_ternary_hull_plot(self):
        """ Test plotting ternary hull. """
        expected_files = ['KSnP_hull.png', 'KSnP_voltage.png']
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
        res_list = glob(REAL_PATH + 'data/hull-KPSn-KP/*.res')
        self.assertEqual(len(res_list), 87, 'Could not find test res files, please check installation...')
        cursor = [res2dict(res)[0] for res in res_list]
        QueryConvexHull(cursor=cursor, elements=['K', 'Sn', 'P'], no_plot=False, png=True, quiet=False, subcmd='voltage',
                        labels=True, label_cutoff=0.05, hull_cutoff=0.1, capmap=True)
        self.assertTrue(os.path.isfile(expected_file))
        for expected_file in expected_files:
            os.remove(expected_file)
