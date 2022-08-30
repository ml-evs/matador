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
from matador.scrapers import res2dict, magres2dict
from matador.hull import QueryConvexHull
from matador.plotting.battery_plotting import plot_voltage_curve
from matador.plotting.pdf_plotting import plot_pdf
from matador.plotting.pxrd_plotting import plot_pxrd
from matador.plotting.magres_plotting import plot_magres
from .utils import MatadorUnitTest

REAL_PATH = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"
ROOT_DIR = os.getcwd()

try:
    import matplotlib  # noqa

    matplotlib.use("Agg")

    MATPLOTLIB_PRESENT = True
except ImportError:
    MATPLOTLIB_PRESENT = False

try:
    import ternary  # noqa

    TERNARY_PRESENT = True
except ImportError:
    TERNARY_PRESENT = False

try:
    import _tkinter  # noqa
except Exception:
    MATPLOTLIB_PRESENT = False


@unittest.skipIf(not MATPLOTLIB_PRESENT, "Skipping plotting tests.")
class SpectralPlotTests(unittest.TestCase):
    """Test Dispersion script."""

    def test_pdis_plot(self):
        """Test combined spectral plots."""
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
        """Test combined spectral plots."""
        os.chdir(REAL_PATH + "/data/dispersion")
        expected_file = "K3P-OQMD_4786-CollCode25550_spectral.png"
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = [
            "dispersion",
            "K3P-OQMD_4786-CollCode25550",
            "--png",
            "--dos_only",
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

    def test_multiseed(self):
        """Test plotting two seed bandstructures on top of each other."""
        os.chdir(REAL_PATH + "/data/bands_files")
        expected_file = "KPSn_spectral.png"
        sys.argv = [
            "dispersion",
            "KPSn",
            "KPSn_2",
            "--dos_only",
            "--cmap",
            "viridis",
            "--png",
            "--band_reorder",
            "--labels",
            "PBE, LDA",
            "--figsize",
            "10",
            "10",
            "--colours",
            "green",
            "red",
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
        """Test combined spectral plots."""
        os.chdir(REAL_PATH + "/data/dispersion")
        sys.argv = [
            "dispersion",
            "K3P-OQMD_4786-CollCode25550",
            "--dos_only",
            "--cmap",
            "viridis",
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
        os.chdir(ROOT_DIR)
        if errored:
            raise error

    def test_phonon_dispersion(self):
        """Test phonon dispersion plot."""
        os.chdir(REAL_PATH + "/data/phonon_dispersion")
        expected_file = "K3P_spectral.png"
        if os.path.isfile(expected_file):
            os.remove(expected_file)
        sys.argv = [
            "dispersion",
            "K3P",
            "--png",
            "-ph",
            "--colours",
            "grey",
            "green",
            "blue",
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

    def test_phonon_ir(self):
        """Test phonon IR/Raman plot."""
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

    def test_projector_scraping(self):
        from matador.plotting.spectral_plotting import _parse_projectors_list

        self.assertEqual(
            _parse_projectors_list("K"),
            [
                ("K", "s", None),
                ("K", "p", None),
                ("K", "d", None),
                ("K", "f", None),
                ("K", None, None),
            ],
        )
        self.assertEqual(
            _parse_projectors_list("K,P"),
            [
                ("K", "s", None),
                ("K", "p", None),
                ("K", "d", None),
                ("K", "f", None),
                ("K", None, None),
                ("P", "s", None),
                ("P", "p", None),
                ("P", "d", None),
                ("P", "f", None),
                ("P", None, None),
            ],
        )
        self.assertEqual(
            _parse_projectors_list("K,P:s"),
            [
                ("K", "s", None),
                ("K", "p", None),
                ("K", "d", None),
                ("K", "f", None),
                ("K", None, None),
                ("P", "s", None),
            ],
        )
        self.assertEqual(
            _parse_projectors_list("123:x,P:s"), [("123", "x", None), ("P", "s", None)]
        )


@unittest.skipIf(not MATPLOTLIB_PRESENT, "Skipping plotting tests.")
class HullPlotTests(MatadorUnitTest):
    """Tests for plotting phase diagrams."""

    def test_binary_hull_plot(self):
        """Test plotting binary hull."""
        cursor = res2dict(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")[0]
        QueryConvexHull(
            cursor=cursor,
            elements=["K", "P"],
            hull_cutoff=0.0,
            plot_kwargs={"plot_fname": "KP_hull_simple", "svg": True},
        )

    def test_binary_battery_plots(self):
        """Test plotting binary hull."""
        cursor = res2dict(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")[0]
        QueryConvexHull(
            cursor=cursor,
            elements=["K", "P"],
            no_plot=False,
            quiet=False,
            voltage=True,
            labels=True,
            label_cutoff=0.05,
            hull_cutoff=0.1,
            volume=True,
            plot_kwargs={"colour_by_source": True},
        )

    def test_voltage_labels(self):
        cursor = res2dict(REAL_PATH + "data/hull-KP-KSnP_pub/*.res")[0]
        hull = QueryConvexHull(
            cursor=cursor, species="KP", no_plot=True, voltage=True, labels=True
        )

        plot_voltage_curve(hull.voltage_data, labels=True)

    @unittest.skipIf(not TERNARY_PRESENT, "Skipping as python-ternary not found")
    def test_ternary_hull_plot(self):
        """Test plotting ternary hull."""
        expected_files = ["KSnP_hull.png", "KSnP_voltage.png"]
        for expected_file in expected_files:
            if os.path.isfile(expected_file):
                os.remove(expected_file)
        res_list = glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        self.assertEqual(
            len(res_list),
            87,
            "Could not find test res files, please check installation...",
        )
        cursor = [res2dict(res)[0] for res in res_list]
        QueryConvexHull(
            cursor=cursor,
            elements=["K", "Sn", "P"],
            no_plot=False,
            png=True,
            quiet=False,
            voltage=True,
            labels=True,
            label_cutoff=0.05,
            hull_cutoff=0.1,
            capmap=True,
        )
        self.assertTrue(os.path.isfile(expected_file))
        for expected_file in expected_files:
            os.remove(expected_file)

    def test_beef_hull_plot(self):
        """Test plotting BEEF hull."""
        from matador.hull import EnsembleHull
        from matador.scrapers import castep2dict

        cursor, s = castep2dict(REAL_PATH + "data/beef_files/*.castep", db=False)
        self.assertEqual(len(s), 0)

        beef_hull = EnsembleHull(
            cursor,
            "_beef",
            elements=["K", "P"],
            num_samples=10,
            energy_key="total_energy_per_atom",
            parameter_key="thetas",
        )

        beef_hull.plot_hull(svg=True)

    def test_td_hull_plot(self):
        from matador.hull.hull_temperature import TemperatureDependentHull
        from matador.scrapers import castep2dict

        cursor, s = castep2dict(
            REAL_PATH + "data/castep_phonon_files/*.castep", db=False
        )
        td_hull = TemperatureDependentHull(
            cursor=cursor, energy_key="total_energy_per_atom"
        )
        td_hull.plot_hull()


@unittest.skipIf(not MATPLOTLIB_PRESENT, "Skipping plotting tests.")
class FingerprintPlotTests(MatadorUnitTest):
    """Test ability to plot PDF and PXRDs."""

    def test_pdf_plot(self):
        structure = res2dict(
            REAL_PATH + "data/res_files/KPSn-OQMD_123456.res", as_model=True
        )[0]
        plot_pdf(structure)
        plot_pdf([structure, structure])

    def test_pxrd_plot(self):
        structure = res2dict(
            REAL_PATH + "data/res_files/KPSn-OQMD_123456.res", as_model=True
        )[0]
        plot_pxrd(structure)
        plot_pdf([structure, structure])


@unittest.skipIf(not MATPLOTLIB_PRESENT, "Skipping plotting tests.")
class MagresPlotTests(MatadorUnitTest):
    """Test ability to plot magres data."""

    def test_magres_plot(self):
        magres, f = magres2dict(
            REAL_PATH + "data/magres_files/*P*.magres", as_model=True
        )
        plot_magres(
            magres,
            species="P",
            line_kwargs={"c": "green"},
        )
        plot_magres(
            magres,
            species="Li",
            broadening_width=0,
            magres_key="chemical_shift_aniso",
            signal_labels=["NaP", "LiP"],
            line_kwargs=[{"lw": 3}, {"ls": "--"}],
        )

        with self.assertRaises(RuntimeError):
            plot_magres(magres, species=None)

        with self.assertRaises(RuntimeError):
            plot_magres(magres, species="K")


@unittest.skipIf(not MATPLOTLIB_PRESENT, "Skipping plotting tests.")
class ConvergencePlotTest(unittest.TestCase):
    """Test the ability to read convergence data and make plots."""

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
        kpt_data = get_convergence_data(
            kpt_files, conv_parameter="kpoints_mp_spacing", species=["Li"]
        )
        cutoff_data = get_convergence_data(
            cutoff_files, conv_parameter="cut_off_energy", species=["Li"]
        )
        data = combine_convergence_data(kpt_data, cutoff_data)
        self.assertEqual(
            data["Li-bcc"]["kpoints_mp_spacing"]["kpoints_mp_spacing"], [0.1, 0.07]
        )
        self.assertEqual(data["Li-bcc"]["cut_off_energy"]["cut_off_energy"], [300, 400])
        values, parameters = get_convergence_values(
            data["Li-bcc"], "cut_off_energy", "formation_energy_per_atom", log=True
        )
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
