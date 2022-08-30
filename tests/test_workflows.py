#!/usr/bin/env python

""" Some tests for high-throughput calculations. """

import unittest
import os
import shutil
import glob

import psutil
import numpy as np

from .utils import MatadorUnitTest, REAL_PATH, detect_program
from matador.compute import ComputeTask
from matador.scrapers import cell2dict, param2dict, phonon2dict, magres2dict


HOSTNAME = os.uname()[1]
PATHS_TO_DEL = ["completed", "bad_castep", "input", "logs", HOSTNAME]
VERBOSITY = 2
EXECUTABLE = "castep"


CASTEP_PRESENT = detect_program(EXECUTABLE)
MPI_PRESENT = detect_program("mpirun")

if CASTEP_PRESENT and MPI_PRESENT:
    NCORES = max(psutil.cpu_count(logical=False) - 2, 1)
else:
    NCORES = 1


@unittest.skipIf(not CASTEP_PRESENT, "CASTEP not found.")
class ElasticWorkflowTest(MatadorUnitTest):
    """Run a elastic workflow calculation."""

    def test_bulk_mod(self):
        for _f in glob.glob(REAL_PATH + "data/elastic_workflow/*"):
            shutil.copy(_f, ".")

        cell_dict, _ = cell2dict("bulk_mod.cell", db=False)
        param_dict, _ = param2dict("bulk_mod.param", db=False)
        _ = ComputeTask(
            res="Si2.res",
            ncores=NCORES,
            nnodes=None,
            node=None,
            cell_dict=cell_dict,
            param_dict=param_dict,
            verbosity=VERBOSITY,
            compute_dir="/tmp/scratch_test",
            workflow_kwargs={"plot": False, "num_volumes": 5},
        )

        self.assertFalse(os.path.isfile("completed/Si2.bib"))
        self.assertTrue(os.path.isfile("completed/Si2.check"))

        self.assertTrue(os.path.isfile("completed/Si2.bulk_mod.results"))
        self.assertTrue(os.path.isfile("completed/Si2.bulk_mod.res"))
        self.assertTrue(os.path.isfile("completed/Si2.bulk_mod.castep"))

        with open("completed/Si2.bulk_mod.results", "r") as f:
            flines = f.readlines()

        B = []
        for line in flines:
            if "bulk modulus" in line:
                B.append(float(line.split()[3]))

        # check all computed bulk mods are between 88-92
        self.assertEqual(len(B), 3)
        self.assertTrue(all(abs(b - 90) < 2) for b in B)

        self.assertFalse(os.path.isfile("completed/Si2.bulk_mod.pdf"))

        self.assertTrue(os.path.isfile("completed/Si2.res"))
        self.assertTrue(os.path.isfile("completed/Si2.geom"))
        self.assertTrue(os.path.isfile("completed/Si2.castep"))

        self.assertTrue(os.path.isfile("completed/Si2.bulk_mod.cell"))
        self.assertFalse(os.path.exists("/tmp/scratch_test"))
        self.assertFalse(os.path.exists("scratch_test_link"))


@unittest.skipIf(not CASTEP_PRESENT, "CASTEP not found.")
class PhononWorkflowTest(MatadorUnitTest):
    """Run a phonon workflow calculation."""

    def test_phonon(self):
        for _f in glob.glob(REAL_PATH + "data/phonon_workflow/*"):
            shutil.copy(_f, ".")

        cell_dict, _ = cell2dict("Si.cell", db=False)
        param_dict, _ = param2dict("Si.param", db=False)
        _ = ComputeTask(
            res="Si2.res",
            ncores=NCORES,
            nnodes=None,
            node=None,
            cell_dict=cell_dict,
            param_dict=param_dict,
            verbosity=VERBOSITY,
            compute_dir="tmpier_tst",
        )

        self.assertFalse(os.path.isfile("completed/Si2.bib"))
        self.assertTrue(os.path.isfile("completed/Si2.check"))

        self.assertTrue(os.path.isfile("completed/Si2.bands"))
        self.assertTrue(os.path.isfile("completed/Si2.castep"))
        self.assertTrue(os.path.isfile("completed/Si2.phonon"))
        self.assertTrue(os.path.isfile("completed/Si2.phonon_dos"))

        phon, s = phonon2dict("completed/Si2.phonon")
        a = 2.7355124
        np.testing.assert_array_almost_equal(
            phon["lattice_cart"], np.array([[0, a, a], [a, 0, a], [a, a, 0]]), decimal=3
        )

        a = 3.869
        np.testing.assert_array_almost_equal(
            phon["lattice_abc"], np.array([[a, a, a], [60, 60, 60]]), decimal=3
        )
        self.assertTrue(s, msg="Failed to read phonon file")
        self.assertGreater(np.min(phon["eigenvalues_q"]), -0.05)

        self.assertTrue(os.path.isfile("completed/Si2.cell"))
        self.assertTrue(os.path.isfile("completed/Si2.res"))


@unittest.skipIf(not CASTEP_PRESENT, "CASTEP not found.")
class MagresWorkflowTest(MatadorUnitTest):
    """Run a magres workflow calculation."""

    def test_magres(self):
        for _f in glob.glob(REAL_PATH + "data/magres_workflow/*"):
            shutil.copy(_f, ".")

        cell_dict, _ = cell2dict("Si.cell", db=False)
        param_dict, _ = param2dict("Si.param", db=False)
        _ = ComputeTask(
            res="Si2.res",
            ncores=NCORES,
            nnodes=None,
            node=None,
            cell_dict=cell_dict,
            param_dict=param_dict,
            verbosity=VERBOSITY,
            compute_dir="tmpier_tst",
            workflow_kwargs={"final_elec_energy_tol": 1e-9},
        )

        self.assertTrue(os.path.isfile("completed/Si2.check"))

        self.assertTrue(os.path.isfile("completed/Si2.bands"))
        self.assertTrue(os.path.isfile("completed/Si2.castep"))
        self.assertTrue(os.path.isfile("completed/Si2.magres"))

        self.assertTrue(os.path.isfile("completed/Si2.cell_magres"))
        self.assertTrue(os.path.isfile("completed/Si2.param_magres"))

        self.assertTrue(os.path.isfile("completed/Si2.cell_scf"))
        self.assertTrue(os.path.isfile("completed/Si2.param_scf"))

        param, s = param2dict("completed/Si2.param_scf")
        self.assertTrue(s, msg="Failed to read param file")
        self.assertEqual(param["elec_energy_tol"], 1e-12)

        param, s = param2dict("completed/Si2.param_magres")
        self.assertEqual(param["elec_energy_tol"], 1e-12)
        self.assertTrue(s, msg="Failed to read param file")

        magres, s = magres2dict("completed/Si2.magres")
        self.assertTrue(s, msg="Failed to read magres file")

        a = 3.866895
        np.testing.assert_array_almost_equal(
            magres["lattice_abc"], np.array([[a, a, a], [60, 60, 60]]), decimal=3
        )

        np.testing.assert_array_almost_equal(
            magres["chemical_shielding_isos"], np.array([129.577, 129.577]), decimal=2
        )

        self.assertTrue(os.path.isfile("completed/Si2.cell"))
        self.assertTrue(os.path.isfile("completed/Si2.res"))


@unittest.skipIf(not CASTEP_PRESENT, "CASTEP not found.")
class SpectralWorkflowTest(MatadorUnitTest):
    """Run a spectral workflow calculation."""

    def test_full_spectral_in_compute_dir(self):
        for _f in glob.glob(REAL_PATH + "data/spectral_workflow/*"):
            shutil.copy(_f, ".")

        cell_dict, _ = cell2dict("Si.cell", db=False)
        param_dict, _ = param2dict("Si.param", db=False)
        _ = ComputeTask(
            res="Si2.res",
            ncores=NCORES,
            nnodes=None,
            node=None,
            cell_dict=cell_dict,
            param_dict=param_dict,
            verbosity=VERBOSITY,
            compute_dir="tmpier_tst",
        )

        self.assertFalse(os.path.isfile("completed/Si2.bib"))
        self.assertTrue(os.path.isfile("completed/Si2.check"))

        self.assertTrue(os.path.isfile("completed/Si2-out.cell_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2-out.cell_dos"))
        self.assertTrue(os.path.isfile("completed/Si2-out.cell_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.adaptive.agr"))
        self.assertTrue(os.path.isfile("completed/Si2.adaptive.dat"))

        self.assertTrue(os.path.isfile("completed/Si2.bands"))
        self.assertTrue(os.path.isfile("completed/Si2.bands_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.bands_dos"))

        self.assertTrue(os.path.isfile("completed/Si2.castep"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.cell"))
        self.assertTrue(os.path.isfile("completed/Si2.cell_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.cell_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.cell_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.cst_esp"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_pdos"))

        self.assertTrue(os.path.isfile("completed/Si2.odi"))
        self.assertTrue(os.path.isfile("completed/Si2.odi_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.odi_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.odi_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.odo"))
        self.assertTrue(os.path.isfile("completed/Si2.odo_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.odo_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.odo_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.param"))
        self.assertTrue(os.path.isfile("completed/Si2.param_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.param_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.param_scf"))
        self.assertTrue(os.path.isfile("completed/Si2.pdis.dat"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos.dat"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.res"))

    def test_dos_only_spectral(self):
        for _f in glob.glob(REAL_PATH + "data/spectral_workflow/*"):
            shutil.copy(_f, ".")

        cell_dict, _ = cell2dict("Si.cell", db=False)
        del cell_dict["spectral_kpoints_path_spacing"]
        param_dict, _ = param2dict("Si.param", db=False)
        _ = ComputeTask(
            res="Si2.res",
            ncores=NCORES,
            nnodes=None,
            node=None,
            cell_dict=cell_dict,
            param_dict=param_dict,
            verbosity=VERBOSITY,
            compute_dir="tmpier_tst",
        )

        self.assertFalse(os.path.isfile("completed/Si2.bib"))
        self.assertTrue(os.path.isfile("completed/Si2.check"))

        self.assertFalse(os.path.isfile("completed/Si2-out.cell_dispersion"))
        self.assertFalse(os.path.isfile("completed/Si2.bands_dispersion"))
        self.assertFalse(os.path.isfile("completed/Si2.castep_dispersion"))
        self.assertFalse(os.path.isfile("completed/Si2.cell_dispersion"))
        self.assertFalse(os.path.isfile("completed/Si2.dome_bin_dispersion"))
        self.assertFalse(os.path.isfile("completed/Si2.param_dispersion"))
        self.assertFalse(os.path.isfile("completed/Si2.pdos_bin_dispersion"))

        self.assertTrue(os.path.isfile("completed/Si2-out.cell_dos"))
        self.assertTrue(os.path.isfile("completed/Si2-out.cell_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.adaptive.agr"))
        self.assertTrue(os.path.isfile("completed/Si2.adaptive.dat"))

        self.assertTrue(os.path.isfile("completed/Si2.bands"))
        self.assertTrue(os.path.isfile("completed/Si2.bands_dos"))

        self.assertTrue(os.path.isfile("completed/Si2.castep"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.cell"))
        self.assertTrue(os.path.isfile("completed/Si2.cell_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.cell_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.cst_esp"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_pdos"))

        self.assertTrue(os.path.isfile("completed/Si2.odi"))
        self.assertTrue(os.path.isfile("completed/Si2.odi_broadening"))
        self.assertFalse(os.path.isfile("completed/Si2.odi_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.odi_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.odo"))
        self.assertTrue(os.path.isfile("completed/Si2.odo_broadening"))
        self.assertFalse(os.path.isfile("completed/Si2.odo_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.odo_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.param"))
        self.assertTrue(os.path.isfile("completed/Si2.param_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.param_scf"))
        self.assertFalse(os.path.isfile("completed/Si2.pdis.dat"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos.dat"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_dos"))
        self.assertFalse(os.path.isfile("completed/Si2.pdos_bin_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.res"))

    def test_full_spectral(self):
        for _f in glob.glob(REAL_PATH + "data/spectral_workflow/*"):
            shutil.copy(_f, ".")

        cell_dict, _ = cell2dict("Si.cell")
        param_dict, _ = param2dict("Si.param", db=False)
        _ = ComputeTask(
            res="Si2.res",
            ncores=NCORES,
            nnodes=None,
            node=None,
            cell_dict=cell_dict,
            param_dict=param_dict,
            verbosity=VERBOSITY,
            compute_dir=None,
        )

        self.assertFalse(os.path.isfile("completed/Si2.bib"))
        self.assertTrue(os.path.isfile("completed/Si2.check"))

        self.assertTrue(os.path.isfile("completed/Si2-out.cell_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2-out.cell_dos"))
        self.assertTrue(os.path.isfile("completed/Si2-out.cell_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.adaptive.agr"))
        self.assertTrue(os.path.isfile("completed/Si2.adaptive.dat"))

        self.assertTrue(os.path.isfile("completed/Si2.bands"))
        self.assertTrue(os.path.isfile("completed/Si2.bands_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.bands_dos"))

        self.assertTrue(os.path.isfile("completed/Si2.castep"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.castep_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.cell"))
        self.assertTrue(os.path.isfile("completed/Si2.cell_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.cell_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.cell_scf"))

        self.assertTrue(os.path.isfile("completed/Si2.cst_esp"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.dome_bin_pdos"))

        self.assertTrue(os.path.isfile("completed/Si2.odi"))
        self.assertTrue(os.path.isfile("completed/Si2.odi_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.odi_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.odi_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.odo"))
        self.assertTrue(os.path.isfile("completed/Si2.odo_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.odo_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.odo_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.param"))
        self.assertTrue(os.path.isfile("completed/Si2.param_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.param_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.param_scf"))
        self.assertTrue(os.path.isfile("completed/Si2.pdis.dat"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos.dat"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_broadening"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_dispersion"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_dos"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_pdis"))
        self.assertTrue(os.path.isfile("completed/Si2.pdos_bin_pdos"))
        self.assertTrue(os.path.isfile("completed/Si2.res"))


if __name__ == "__main__":
    unittest.main()
