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
        file_exists = os.path.isfile('K3P-OQMD_4786-CollCode25550_spectral.png')
        os.remove(expected_file)
        os.chdir(ROOT_DIR)
        self.assertFalse(error)
        self.assertTrue(file_exists)


class HullPlotTests(unittest.TestCase):
    """ Tests for plotting phase diagrams. """
    def test_binary_hull_plot(self):
        """ Test plotting binary hull. """
        expected_files = ['KP_hull.png', 'KP_voltage.png']
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
        res_list = glob(REAL_PATH + 'data/hull-KP-KSnP_pub/*.res')
        self.assertEqual(len(res_list), 295, 'Could not find test res files, please check installation...')
        cursor = [res2dict(res)[0] for res in res_list]
        QueryConvexHull(cursor=cursor, elements=['K', 'P'], no_plot=False, png=True, quiet=False, subcmd='voltage')
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
        QueryConvexHull(cursor=cursor, elements=['K', 'Sn', 'P'], no_plot=False, png=True, quiet=False, subcmd='voltage')
        self.assertTrue(os.path.isfile(expected_file))
        for expected_file in expected_files:
            os.remove(expected_file)
