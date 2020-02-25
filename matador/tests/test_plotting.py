#!/usr/bin/env python

""" These tests only check whether plots are created,
not that they look correct!

"""

import unittest
import os
import sys
from glob import glob
import numpy as np
import matador.cli.dispersion
from matador.scrapers import res2dict
from matador.hull import QueryConvexHull
from matador.plotting.battery_plotting import plot_voltage_curve
from matador.utils.chem_utils import get_formula_from_stoich

REAL_PATH = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"
ROOT_DIR = os.getcwd()

try:
    import matplotlib  # noqa

    MATPLOTLIB_PRESENT = True
except ImportError:
    MATPLOTLIB_PRESENT = False

try:
    import _tkinter  # noqa
except Exception:
    MATPLOTLIB_PRESENT = False


@unittest.skipIf(not MATPLOTLIB_PRESENT, "Skipping plotting tests.")
class SpectralPlotTests(unittest.TestCase):
    """ Test Dispersion script. """

    def test_pdis_plot(self):
        """ Test combined spectral plots. """
        os.chdir(REAL_PATH + "/data/dispersion")
        expected_file = "K3P-OQMD_4786-CollCode25550_spectral.png"
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = [
            "dispersion",
            "K3P-OQMD_4786-CollCode25550",
            "--png",
            "-scale",
            "10",
            "-interp",
            "2",
            "-pw",
            "-5",
            "5",
            "--gap",
            "--preserve_kspace_distance",
            "--figsize",
            "10",
            "10",
        ]
        errored = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            errored = True
            error = exc
        file_exists = os.path.isfile(expected_file)
        if file_exists:
            os.remove(expected_file)
        os.chdir(ROOT_DIR)
        if errored:
            raise error
        self.assertTrue(file_exists)

    def test_dos_only(self):
        """ Test combined spectral plots. """
        os.chdir(REAL_PATH + "/data/dispersion")
        expected_file = "K3P-OQMD_4786-CollCode25550_spectral.png"
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = ["dispersion", "K3P-OQMD_4786-CollCode25550", "--png", "--dos_only", "--figsize", "10", "10"]
        errored = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            errored = True
            error = exc
        file_exists = os.path.isfile(expected_file)
        if file_exists:
            os.remove(expected_file)
        os.chdir(ROOT_DIR)
        if errored:
            raise error
        self.assertTrue(file_exists)

    def test_multiseed(self):
        """ Test plotting two seed bandstructures on top of each other. """
        os.chdir(REAL_PATH + "/data/bands_files")
        expected_file = "KPSn_spectral.png"
        sys.argv = [
            "dispersion",
            "KPSn",
            "KPSn_2.bands",
            "--dos_only",
            "--cmap",
            "viridis",
            "--png",
            "--labels",
            "PBE, LDA",
            "--figsize",
            "10",
            "10",
        ]
        errored = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            errored = True
            error = exc
        file_exists = os.path.isfile(expected_file)
        if file_exists:
            os.remove(expected_file)
        os.chdir(ROOT_DIR)
        if errored:
            raise error
        self.assertTrue(file_exists)

    def test_x11_no_fail(self):
        """ Test combined spectral plots. """
        os.chdir(REAL_PATH + "/data/dispersion")
        sys.argv = ["dispersion", "K3P-OQMD_4786-CollCode25550", "--dos_only", "--cmap", "viridis", "--figsize", "10", "10"]
        errored = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            errored = True
            error = exc
        os.chdir(ROOT_DIR)
        if errored:
            raise error

    def test_phonon_dispersion(self):
        """ Test phonon dispersion plot. """
        os.chdir(REAL_PATH + "/data/phonon_dispersion")
        expected_file = "K3P_spectral.png"
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = ["dispersion", "K3P", "--png", "-ph", "--band_colour", "k", "--figsize", "10", "10"]
        errored = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            errored = True
            error = exc
        file_exists = os.path.isfile(expected_file)
        if file_exists:
            os.remove(expected_file)
        os.chdir(ROOT_DIR)
        if errored:
            raise error
        self.assertTrue(file_exists)

    def test_phonon_ir(self):
        """ Test phonon IR/Raman plot. """
        os.chdir(REAL_PATH + "/data/phonon_ir")
        expected_file = "h-BN_IRR_ir.svg"
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = ["dispersion", "h-BN_IRR", "--svg", "-ir", "--figsize", "5", "5"]
        errored = False
        try:
            matador.cli.dispersion.main()
        except Exception as exc:
            errored = True
            error = exc
        file_exists = os.path.isfile(expected_file)
        if file_exists:
            os.remove(expected_file)
        os.chdir(ROOT_DIR)
        if errored:
            raise error
        self.assertTrue(file_exists)


@unittest.skipIf(not MATPLOTLIB_PRESENT, "Skipping plotting tests.")
class HullPlotTests(unittest.TestCase):
    """ Tests for plotting phase diagrams. """

    def test_binary_hull_plot(self):
        """ Test plotting binary hull. """
        expected_files = ["KP_hull_simple.svg"]
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
        res_list = glob(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")
        self.assertEqual(len(res_list), 295, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        QueryConvexHull(
            cursor=cursor,
            elements=["K", "P"],
            subcmd="hull",
            svg=True,
            hull_cutoff=0.0,
            plot_kwargs={"plot_fname": "KP_hull_simple", "svg": True},
        )
        for expected_file in expected_files:
            self.assertTrue(os.path.isfile(expected_file))
        for expected_file in expected_files:
            os.remove(expected_file)

    def test_binary_battery_plots(self):
        """ Test plotting binary hull. """
        expected_files = ["KP_hull.png", "KP_voltage.png", "KP_volume.png"]
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
        res_list = glob(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")
        self.assertEqual(len(res_list), 295, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        QueryConvexHull(
            cursor=cursor,
            elements=["K", "P"],
            no_plot=False,
            png=True,
            quiet=False,
            subcmd="voltage",
            labels=True,
            label_cutoff=0.05,
            hull_cutoff=0.1,
            volume=True,
            plot_kwargs={"colour_by_source": True},
        )
        for expected_file in expected_files:
            self.assertTrue(os.path.isfile(expected_file))
        for expected_file in expected_files:
            os.remove(expected_file)

    def test_voltage_labels(self):
        expected_files = ["KP_voltage.png"]
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
        res_list = glob(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")
        self.assertEqual(len(res_list), 295, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        hull = QueryConvexHull(cursor=cursor, species="KP", no_plot=True, subcmd="voltage", labels=True)

        label_cursor = []
        plot_voltage_curve(hull, label_cursor=label_cursor, png=True)
        labels = [get_formula_from_stoich(doc["stoichiometry"], elements=hull.elements, tex=False) for doc in label_cursor]

        self.assertEqual(labels, ["KP7", "K3P7", "K2P3", "KP", "K5P4"])
        for expected_file in expected_files:
            self.assertTrue(os.path.isfile(expected_file))
        for expected_file in expected_files:
            os.remove(expected_file)

    def test_ternary_hull_plot(self):
        """ Test plotting ternary hull. """
        expected_files = ["KSnP_hull.png", "KSnP_voltage.png"]
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(len(res_list), 87, "Could not find test res files, please check installation...")
        cursor = [res2dict(res)[0] for res in res_list]
        QueryConvexHull(
            cursor=cursor,
            elements=["K", "Sn", "P"],
            no_plot=False,
            png=True,
            quiet=False,
            subcmd="voltage",
            labels=True,
            label_cutoff=0.05,
            hull_cutoff=0.1,
            capmap=True,
        )
        self.assertTrue(os.path.isfile(expected_file))
        for expected_file in expected_files:
            os.remove(expected_file)

    def test_beef_hull_plot(self):
        """ Test plotting BEEF hull. """
        from matador.hull import EnsembleHull
        from matador.scrapers import castep2dict

        expected_file = "KP_beef_hull.svg"
        if os.path.isfile(expected_file):
            os.remove(expected_file)

        cursor, s = castep2dict(REAL_PATH + "data/beef_files/*.castep", db=False)
        self.assertEqual(len(s), 0)

        beef_hull = EnsembleHull(
            cursor, "_beef", elements=["K", "P"], num_samples=10, energy_key="total_energy_per_atom", parameter_key="thetas"
        )

        beef_hull.plot_hull(svg=True)
        self.assertTrue(os.path.isfile(expected_file))
        os.remove(expected_file)


@unittest.skipIf(not MATPLOTLIB_PRESENT, "Skipping plotting tests.")
class ConvergencePlotTest(unittest.TestCase):
    """ Test the ability to read convergence data and make plots. """

    def setUp(self):
        os.chdir(REAL_PATH + "/data/convergence/")

    def tearDown(self):
        os.chdir(ROOT_DIR)

    def test_scraping_and_plotting(self):
        from matador.plotting.convergence_plotting import (
            get_convergence_data,
            get_convergence_files,
            combine_convergence_data,
            get_convergence_values,
        )
        from matador.plotting.convergence_plotting import plot_cutoff_kpt_grid

        kpt_files = get_convergence_files("completed_kpts")
        cutoff_files = get_convergence_files("completed_cutoff")
        kpt_data = get_convergence_data(kpt_files, conv_parameter="kpoints_mp_spacing", species=["Li"])
        cutoff_data = get_convergence_data(cutoff_files, conv_parameter="cut_off_energy", species=["Li"])
        data = combine_convergence_data(kpt_data, cutoff_data)
        self.assertEqual(data["Li-bcc"]["kpoints_mp_spacing"]["kpoints_mp_spacing"], [0.1, 0.07])
        self.assertEqual(data["Li-bcc"]["cut_off_energy"]["cut_off_energy"], [300, 400])
        values, parameters = get_convergence_values(data["Li-bcc"], "cut_off_energy", "formation_energy_per_atom", log=True)
        self.assertEqual(parameters.tolist(), [300.0, 400.0])
        self.assertAlmostEqual(values.tolist()[0], 0.7291198427497395, places=6)
        self.assertEqual(values.tolist()[1], -np.inf)
        self.data = data

        expected_files = ["conv.svg"]
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)

        plot_cutoff_kpt_grid(self.data, svg=True)

        for file in expected_files:
            self.assertTrue(os.path.isfile(file))

        for expected_file in expected_files:
            os.remove(expected_file)
