#!/usr/bin/env python
""" Test file scraping and writing functionality. """

import json
import os
import glob
import itertools
import numpy as np

from matador.scrapers import castep2dict, res2dict, cell2dict
from matador.scrapers import (
    cif2dict,
    param2dict,
    phonon2dict,
    optados2dict,
    phonon_dos2dict,
)
from matador.scrapers import arbitrary2dict, bands2dict, pwout2dict, magres2dict
from matador.scrapers.castep_scrapers import usp2dict, get_seed_metadata
from matador.export import doc2res, doc2param, doc2cell, query2files
from matador.orm.spectral import (
    ElectronicDispersion,
    ElectronicDOS,
    VibrationalDispersion,
    VibrationalDOS,
)
from matador.utils.chem_utils import INVERSE_CM_TO_EV
from .utils import REAL_PATH, MatadorUnitTest

VERBOSITY = 10


class CellScraperTests(MatadorUnitTest):
    """ Test cell scraper functions. """

    def test_standard_cell_scraper(self):
        cell_fname = REAL_PATH + "data/LiP2Zn-0bm995-a_9-out.cell"
        self.assertTrue(
            os.path.isfile(cell_fname),
            msg="Failed to open test case {} - please check installation.".format(
                cell_fname
            ),
        )
        test_dict, s = cell2dict(
            cell_fname, db=False, lattice=True, verbosity=VERBOSITY
        )
        self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
        self.assertEqual(
            test_dict["lattice_cart"][0][0],
            9.83262140721165,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][1][1],
            5.96357780025648,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][2][2],
            4.39895761828278,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][1][0],
            -0.115688800302997,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["symmetry_tol"], 0.001, msg="Failed to read symmetry tolerance."
        )
        self.assertEqual(
            test_dict["kpoints_mp_grid"],
            [2, 3, 4],
            msg="Failed to read kpoint grid {}".format(test_dict["kpoints_mp_grid"]),
        )
        self.assertEqual(
            test_dict["species_pot"]["Li"], "Li_00PBE.usp", msg="Failed to read pspots."
        )
        self.assertEqual(
            test_dict["species_pot"]["P"], "P_00PBE.usp", msg="Failed to read pspots."
        )
        self.assertEqual(
            test_dict["species_pot"]["Zn"], "Zn_00PBE.usp", msg="Failed to read pspots."
        )
        # test that lattice_vec only read when outcell is true
        test_dict, s = cell2dict(
            cell_fname, db=False, lattice=False, verbosity=VERBOSITY
        )
        self.assertTrue(test_dict.get("lattice_cart") is None)

    def test_cell_outcell(self):
        cell_fname = REAL_PATH + "data/Li2C2-out.cell"
        self.assertTrue(
            os.path.isfile(cell_fname),
            msg="Failed to open test case {} - please check installation.".format(
                cell_fname
            ),
        )
        test_dict, s = cell2dict(
            cell_fname, db=False, lattice=True, verbosity=VERBOSITY
        )
        self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
        self.assertEqual(test_dict["cell_constraints"], [[1, 1, 3], [4, 4, 6]])
        self.assertEqual(
            test_dict["external_pressure"],
            [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
        )
        tmp_name = "tmp.cell"
        doc2cell(test_dict, tmp_name)
        new_test_dict, s = cell2dict(
            tmp_name, db=False, lattice=True, verbosity=VERBOSITY
        )
        new_test_dict["source"] = test_dict["source"]
        self.assertEqual(
            test_dict["external_pressure"], new_test_dict["external_pressure"]
        )

    def test_cell_phonon(self):
        cell_fname = REAL_PATH + "data/K5P4-phonon.cell"
        self.assertTrue(
            os.path.isfile(cell_fname),
            msg="Failed to open test case {} - please check installation.".format(
                cell_fname
            ),
        )
        test_dict, s = cell2dict(
            cell_fname, db=False, lattice=True, verbosity=VERBOSITY
        )
        self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
        self.assertEqual(
            test_dict["lattice_cart"][0][0],
            11.4518745146637,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][1][1],
            5.09448137301246,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][2][2],
            9.18378851243459,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][1][0], 0.0, msg="Failed to read lattice vectors."
        )
        self.assertEqual(
            test_dict["symmetry_tol"], 0.0001, msg="Failed to read symmetry tolerance."
        )
        self.assertEqual(
            test_dict["kpoints_mp_spacing"],
            0.03,
            msg="Failed to read kpoint grid {}".format(test_dict["kpoints_mp_spacing"]),
        )
        self.assertEqual(
            test_dict["phonon_kpoint_mp_grid"],
            [2, 2, 2],
            msg="Failed to read kpoint grid {}".format(
                test_dict["phonon_kpoint_mp_grid"]
            ),
        )
        self.assertEqual(
            test_dict["phonon_kpoint_mp_offset"],
            [0.25, 0.25, 0.25],
            msg="Failed to read kpoint grid {}".format(
                test_dict["phonon_kpoint_mp_offset"]
            ),
        )
        self.assertEqual(
            test_dict["phonon_fine_kpoint_mp_spacing"],
            0.02,
            msg="Failed to read kpoint {}".format(
                test_dict["phonon_fine_kpoint_mp_spacing"]
            ),
        )
        self.assertEqual(
            test_dict["phonon_fine_kpoint_path_spacing"],
            0.01,
            msg="Failed to read kpoint {}".format(
                test_dict["phonon_fine_kpoint_path_spacing"]
            ),
        )
        self.assertEqual(
            test_dict["species_pot"]["K"],
            "2|1.5|9|10|11|30U:40:31(qc=6)",
            msg="Failed to read pspots.",
        )
        self.assertEqual(
            test_dict["species_pot"]["P"],
            "3|1.8|4|4|5|30:31:32",
            msg="Failed to read pspots.",
        )
        self.assertEqual(
            test_dict["hubbard_u"]["K"]["s"], 2, msg="Failed to read Hubbard U block."
        )
        self.assertEqual(
            test_dict["hubbard_u"]["P"]["p"], 3, msg="Failed to read Hubbard U block."
        )
        self.assertEqual(
            test_dict["hubbard_u"]["U"]["d"],
            10.101,
            msg="Failed to read Hubbard U block.",
        )
        self.assertTrue(test_dict["snap_to_symmetry"])
        self.assertTrue(test_dict["symmetry_generate"])
        self.assertEqual(test_dict["phonon_supercell_matrix"][0], [3, 0, 1])
        self.assertEqual(test_dict["phonon_supercell_matrix"][1], [0, 3, 0])
        self.assertEqual(test_dict["phonon_supercell_matrix"][2], [0, 0, 9])

    def test_cell_failure(self):
        cell_fname = REAL_PATH + "data/K5P4-phonon_bodged.cell"
        self.assertTrue(
            os.path.isfile(cell_fname),
            msg="Failed to open test case {} - please check installation.".format(
                cell_fname
            ),
        )
        test_dict, s = cell2dict(cell_fname, db=True, lattice=True, verbosity=VERBOSITY)
        self.assertFalse(s, msg=test_dict)

    def test_cell_spin(self):
        cell_fname = REAL_PATH + "data/cell_files/spin_test.cell"
        self.assertTrue(
            os.path.isfile(cell_fname),
            msg="Failed to open test case {} - please check installation.".format(
                cell_fname
            ),
        )
        test_dict, s = cell2dict(
            cell_fname, db=False, lattice=True, positions=True, verbosity=VERBOSITY
        )
        self.assertEqual(test_dict["species_pot"]["library"], "QC5")
        self.assertEqual(test_dict["atom_types"], ["H", "C", "H", "C", "H", "H"])
        self.assertEqual(test_dict["lattice_cart"][0], [10, 0, 0])
        self.assertEqual(test_dict["lattice_cart"][1], [0, 10, 0])
        self.assertEqual(test_dict["lattice_cart"][2], [0, 0, 10])
        self.assertEqual(test_dict["lattice_abc"][0], [10, 10, 10])
        self.assertEqual(test_dict["lattice_abc"][1], [90, 90, 90])
        self.assertEqual(
            test_dict["atomic_init_spins"], [None, 0.32675521, None, -0.1234, None, 1.0]
        )

    def test_cell_kpoint_path(self):
        cell_name = REAL_PATH + "data/cell_files/kpoint_path.cell"
        cell, s = cell2dict(cell_name, db=False)

        self.assertTrue(s)
        self.assertEqual(
            cell["spectral_kpoints_path"],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0], [0.5, 0.0, 0.0]],
        )
        self.assertEqual(
            cell["spectral_kpoints_path_labels"], ["$\\Gamma$", "Z", "$Y$", "X"]
        )
        self.assertEqual(cell["spectral_kpoints_path_spacing"], 0.02)
        self.assertEqual(
            cell["phonon_fine_kpoint_path"],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0], [0.5, 0.0, 0.0]],
        )
        self.assertEqual(
            cell["phonon_fine_kpoint_path_labels"], ["$\\Gamma$", "Z", "$Y$", "X"]
        )
        self.assertEqual(cell["phonon_fine_kpoint_path_spacing"], 0.01)

    def test_cell_positions_abs(self):
        cell_name = REAL_PATH + "data/cell_files/npm.cell"
        cell, s = cell2dict(cell_name, db=False)

        self.assertTrue(s, msg="Failed entirely: {}".format(cell))
        self.assertEqual(cell["lattice_cart"][0], [21.84, 0, 0])
        self.assertEqual(cell["lattice_cart"][1], [0, 16.38, 0])
        self.assertEqual(cell["lattice_cart"][2], [0, 0, 40.46])

        self.assertEqual(cell["positions_abs"][0], [0, 0, 0])
        self.assertEqual(cell["positions_abs"][95], [17.745, 15.015, -4.095])
        self.assertEqual(cell["positions_abs"][-1], [13.65, 13.65, 10.92])

        np.testing.assert_array_almost_equal(cell["positions_frac"][0], [0, 0, 0])
        np.testing.assert_array_almost_equal(
            cell["positions_frac"][95], [0.8125, 0.9166666, 0.8987889]
        )
        np.testing.assert_array_almost_equal(
            cell["positions_frac"][-1], [0.625, 0.8333333, 0.26989619]
        )

    def test_cell_ionic_cell_constraints(self):
        cell_name = REAL_PATH + "data/cell_files/ionic_constraints.cell"
        cell, s = cell2dict(cell_name, db=False)

        self.assertTrue(s)
        self.assertEqual(cell["ionic_constraints"][0], "1 C 1 0 0 1")
        self.assertEqual(cell["ionic_constraints"][1], "2 C 1 0 1 0")
        self.assertEqual(cell["ionic_constraints"][2], "3 C 1 1 0 0")
        self.assertEqual(cell["cell_constraints"][0], [1, 1, 3])
        self.assertEqual(cell["cell_constraints"][1], [0, 0, 0])


class CastepScraperTests(MatadorUnitTest):
    """ Test CASTEP scrapers. """

    def test_castep16(self):
        castep_fname = REAL_PATH + "data/Na3Zn4-swap-ReOs-OQMD_759599.castep"
        failed_open = False
        try:
            f = open(castep_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    castep_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
            self.assertEqual(
                test_dict["pressure"], 0.0763, msg="Failed to read pressure!"
            )
            self.assertEqual(
                test_dict["enthalpy"], -2.15036930e4, msg="Failed to read enthalpy!"
            )
            self.assertEqual(test_dict["num_atoms"], 14, msg="Wrong number of atoms!")
            self.assertTrue(
                ["Na", 3] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertTrue(
                ["Zn", 4] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertEqual(
                test_dict["cell_volume"], 288.041941, msg="Wrong cell volume!"
            )
            self.assertEqual(test_dict["space_group"], "Pm", msg="Wrong space group!")
            self.assertEqual(
                test_dict["lattice_abc"][0][0], 9.039776, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][1], 9.045651, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][2], 4.068682, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][0], 90, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][1], 90, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][2],
                59.971185,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["geom_force_tol"], 0.05, msg="Wrong geom force tol"
            )
            self.assertEqual(test_dict["castep_version"], "16.11")
            self.assertEqual(test_dict["_castep_commit"], "203e84763863+")
            self.assertAlmostEqual(test_dict["total_time_secs"], 1291.14, places=2)
            self.assertEqual(test_dict["geom_iter"], 8)
            self.assertEqual(
                test_dict["external_pressure"],
                [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            )
            self.assertEqual(test_dict["estimated_mem_per_process_MB"], 345.1)
            self.assertEqual(test_dict["peak_mem_MB"], int(675372 / 1024))
            self.assertEqual(test_dict["num_mpi_processes"], 15)

    def test_castep17(self):
        castep_fname = REAL_PATH + "data/KP-castep17.castep"
        failed_open = False
        try:
            f = open(castep_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    castep_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
            self.assertEqual(
                test_dict["pressure"], 0.0180, msg="Failed to read pressure!"
            )
            self.assertEqual(
                test_dict["enthalpy"], -5.98055077e3, msg="Failed to read enthalpy!"
            )
            self.assertEqual(test_dict["num_atoms"], 9, msg="Wrong number of atoms!")
            self.assertTrue(
                ["P", 2] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertTrue(
                ["K", 7] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertEqual(
                test_dict["cell_volume"], 522.226927, msg="Wrong cell volume!"
            )
            self.assertEqual(test_dict["space_group"], "Pm", msg="Wrong space group!")
            self.assertEqual(
                test_dict["lattice_abc"][0][0],
                10.231976,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][1], 5.024837, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][2],
                10.186949,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][0],
                90.000000,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][1],
                94.373377,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][2],
                90.000000,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["geom_force_tol"], 0.01, msg="Wrong geom force tol"
            )
            self.assertEqual(test_dict["castep_version"], "17.21")
            self.assertEqual(
                test_dict["_compiler_architecture"], "linux_x86_64_ifort17"
            )
            self.assertEqual(test_dict["_castep_commit"], "056e886bd5a1+")
            self.assertEqual(test_dict["optimised"], True)
            self.assertEqual(test_dict["estimated_mem_per_process_MB"], 300.1)
            self.assertEqual(
                test_dict["species_pot"]["K"],
                "2|1.5|9|10|11|30U:40:31(qc=6)",
                msg="Failed to scrape K_OTF.usp file",
            )
            self.assertEqual(
                test_dict["species_pot"]["P"],
                "3|1.8|4|4|5|30:31:32",
                msg="Failed to scrape P_OTF.usp file",
            )

    def test_castep_single_atom_edgecase(self):
        castep_fname = REAL_PATH + "data/castep_files/Na-edgecase-CollCode10101.castep"
        failed_open = False
        try:
            f = open(castep_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    castep_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
            self.assertEqual(
                test_dict["pressure"], -0.0966, msg="Failed to read pressure!"
            )
            self.assertEqual(
                test_dict["enthalpy"], -1.30423371e3, msg="Failed to read enthalpy!"
            )
            self.assertEqual(test_dict["positions_frac"], [[0, 0, 0]])
            self.assertEqual(test_dict["forces"], [[0, 0, 0]])
            self.assertEqual(
                test_dict["enthalpy"], -1.30423371e3, msg="Failed to read enthalpy!"
            )
            self.assertEqual(
                test_dict["total_energy"],
                -1304.223019263,
                msg="Failed to read total energy!",
            )
            self.assertEqual(
                test_dict["total_energy_per_atom"],
                -1304.223019263,
                msg="Failed to read total energy!",
            )
            self.assertEqual(
                test_dict["smeared_free_energy"],
                -1304.233706274,
                msg="Failed to read free energy!",
            )
            self.assertEqual(test_dict["num_atoms"], 1, msg="Wrong number of atoms!")
            self.assertTrue(
                ["Na", 1] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertEqual(
                len(test_dict["stoichiometry"]), 1, msg="Wrong stoichiometry!"
            )
            self.assertEqual(
                test_dict["cell_volume"], 36.761902, msg="Wrong cell volume!"
            )
            self.assertEqual(
                test_dict["space_group"], "Im-3m", msg="Wrong space group!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][0], 3.628050, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][1], 3.628050, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][2], 3.628050, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][0],
                109.471221,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][1],
                109.471221,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][2],
                109.471221,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["geom_force_tol"], 0.05, msg="Wrong geom force tol"
            )
            self.assertEqual(test_dict["castep_version"], "16.1")
            self.assertEqual(test_dict["species_pot"]["Na"], "Na_00PBE.usp")
            self.assertEqual(test_dict["icsd"], 10101)
            self.assertEqual(
                test_dict["_compiler_architecture"], "linux_x86_64_ifort14"
            )
            self.assertEqual(test_dict["_castep_commit"], "2756eb6097bf+")

            int_dict, s = castep2dict(
                castep_fname, db=False, intermediates=True, verbosity=VERBOSITY
            )
            for key in test_dict:
                self.assertEqual(test_dict[key], int_dict[key])

            self.assertEqual(len(int_dict["intermediates"]), 45)
            for i in range(45):
                self.assertEqual(int_dict["intermediates"][i]["forces"], [[0, 0, 0]])
                self.assertEqual(
                    int_dict["intermediates"][i]["positions_frac"], [[0, 0, 0]]
                )
                self.assertEqual(int_dict["intermediates"][i]["atom_types"], ["Na"])
            self.assertEqual(
                int_dict["intermediates"][-1]["total_energy"], -1304.223019263
            )
            self.assertEqual(
                int_dict["intermediates"][-1]["smeared_free_energy"], -1304.233706274
            )
            self.assertEqual(
                int_dict["intermediates"][-7]["total_energy"], -1304.222982442
            )
            self.assertEqual(
                int_dict["intermediates"][-7]["smeared_free_energy"], -1304.233677344
            )
            self.assertEqual(
                int_dict["intermediates"][-1]["total_energy_per_atom"], -1304.223019263
            )
            self.assertEqual(
                int_dict["intermediates"][-1]["smeared_free_energy_per_atom"],
                -1304.233706274,
            )
            self.assertEqual(
                int_dict["intermediates"][-7]["total_energy_per_atom"], -1304.222982442
            )
            self.assertEqual(
                int_dict["intermediates"][-7]["smeared_free_energy_per_atom"],
                -1304.233677344,
            )
            self.assertEqual(int_dict["geom_iter"], 44)

    def test_castep_unoptimised(self):
        castep_fname = REAL_PATH + "data/castep_files/TiO2_unconverged-MP-10101.castep"
        failed_open = False
        try:
            f = open(castep_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    castep_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertFalse(s, msg="Should have failed with db=True, but didn't!")
            self.assertTrue(
                isinstance(test_dict, Exception),
                msg="Should have returned error message!",
            )
            test_dict, s = castep2dict(castep_fname, db=False, verbosity=VERBOSITY)
            self.assertTrue(s, msg="Should have succeeded with db=False, but didn't!")
            self.assertTrue(
                isinstance(test_dict, dict), msg="Should have returned dict!"
            )
            self.assertEqual(test_dict["total_energy"], -12479.86611705)
            self.assertEqual(test_dict["num_atoms"], 12)
            self.assertEqual(
                test_dict["pressure"], 0.9455, msg="Failed to read pressure!"
            )
            self.assertTrue(
                ["Ti", 1] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertTrue(
                ["O", 2] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertEqual(
                test_dict["cell_volume"], 127.269750, msg="Wrong cell volume!"
            )
            self.assertEqual(test_dict["space_group"], "Pmmm")
            self.assertEqual(
                test_dict["lattice_abc"][0][0], 4.026041, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][1], 7.906524, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][0][2], 3.998172, msg="Wrong lattice constants!"
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][0],
                90.000000,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][1],
                90.000000,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(
                test_dict["lattice_abc"][1][2],
                90.000000,
                msg="Wrong lattice constants!",
            )
            self.assertEqual(test_dict["optimised"], False)
            self.assertEqual(test_dict["geom_force_tol"], 0.05)
            self.assertEqual(test_dict["castep_version"], "18.1")
            self.assertEqual(
                test_dict["species_pot"]["Ti"], "3|1.9|8|9|10|30U:40:31:32(qc=5)"
            )
            self.assertEqual(
                test_dict["species_pot"]["O"], "2|1.5|12|13|15|20:21(qc=5)"
            )
            self.assertEqual(test_dict["mp_id"], 10101)

    def test_file_not_found(self):
        """ Ensure that FileNotFound errors fail gracefully. """
        error = False
        try:
            res, s = res2dict("___not_a_file")
        except FileNotFoundError:
            error = True
        self.assertTrue(error)

        castep_fname = []
        castep_fname += [REAL_PATH + "data/castep_files/NaP_intermediates.castep"]
        castep_fname += [REAL_PATH + "data/___not_a_file"]
        castep_fname += [REAL_PATH + "data/KP-castep17.castep"]
        castep_fname += [REAL_PATH + "data/Na3Zn4-swap-ReOs-OQMD_759599.castep"]

        error = False
        try:
            cursor, failures = castep2dict(castep_fname, db=True)
        except FileNotFoundError:
            error = True

    def test_history(self):
        castep_fname = (
            REAL_PATH + "data/castep_files/Na3Zn4-swap-ReOs-OQMD_759599.history"
        )
        test_dict, s = castep2dict(castep_fname, db=True)
        self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
        self.assertEqual(test_dict["pressure"], 0.0763, msg="Failed to read pressure!")
        self.assertEqual(
            test_dict["enthalpy"], -2.15036930e4, msg="Failed to read enthalpy!"
        )
        self.assertEqual(test_dict["num_atoms"], 14, msg="Wrong number of atoms!")
        self.assertTrue(
            ["Na", 3] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
        )
        self.assertTrue(
            ["Zn", 4] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
        )
        self.assertEqual(test_dict["cell_volume"], 288.041941, msg="Wrong cell volume!")
        self.assertEqual(test_dict["space_group"], "Pm", msg="Wrong space group!")
        self.assertEqual(
            test_dict["lattice_abc"][0][0], 9.039776, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][0][1], 9.045651, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][0][2], 4.068682, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][1][0], 90, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][1][1], 90, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][1][2], 59.971185, msg="Wrong lattice constants!"
        )
        self.assertEqual(test_dict["geom_force_tol"], 0.05, msg="Wrong geom force tol")
        self.assertEqual(test_dict["castep_version"], "16.11")
        self.assertEqual(test_dict["estimated_mem_per_process_MB"], 345.1)

    def test_history_gz(self):
        castep_fname = (
            REAL_PATH + "data/castep_files/Na3Zn4-swap-ReOs-OQMD_759599.history.gz"
        )
        test_dict, s = castep2dict(castep_fname, db=True)
        self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
        self.assertEqual(test_dict["pressure"], 0.0763, msg="Failed to read pressure!")
        self.assertEqual(
            test_dict["enthalpy"], -2.15036930e4, msg="Failed to read enthalpy!"
        )
        self.assertEqual(test_dict["num_atoms"], 14, msg="Wrong number of atoms!")
        self.assertTrue(
            ["Na", 3] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
        )
        self.assertTrue(
            ["Zn", 4] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
        )
        self.assertEqual(test_dict["cell_volume"], 288.041941, msg="Wrong cell volume!")
        self.assertEqual(test_dict["space_group"], "Pm", msg="Wrong space group!")
        self.assertEqual(
            test_dict["lattice_abc"][0][0], 9.039776, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][0][1], 9.045651, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][0][2], 4.068682, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][1][0], 90, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][1][1], 90, msg="Wrong lattice constants!"
        )
        self.assertEqual(
            test_dict["lattice_abc"][1][2], 59.971185, msg="Wrong lattice constants!"
        )
        self.assertEqual(test_dict["geom_force_tol"], 0.05, msg="Wrong geom force tol")
        self.assertEqual(test_dict["castep_version"], "16.11")

    def test_castep_intermediates(self):
        castep_fname = REAL_PATH + "data/castep_files/NaP_intermediates.castep"
        self.assertTrue(
            os.path.isfile(castep_fname),
            msg="Failed to open test case {} - please check installation.".format(
                castep_fname
            ),
        )
        test_dict, s = castep2dict(
            castep_fname, db=False, intermediates=True, verbosity=VERBOSITY
        )
        self.assertTrue(s, msg="Should have succeeded with db=False, but didn't!")
        final_dict, s = castep2dict(
            castep_fname, db=True, intermediates=False, verbosity=VERBOSITY
        )
        self.assertTrue(s)
        for key in final_dict:
            self.assertEqual(
                final_dict[key], test_dict[key], msg="{} didn't match".format(key)
            )
        self.assertEqual(test_dict["intermediates"][0]["total_energy"], -8537.190779552)
        self.assertEqual(test_dict["intermediates"][1]["total_energy"], -8538.161269966)
        self.assertEqual(
            test_dict["intermediates"][-1]["total_energy"], -8546.922111847
        )
        self.assertEqual(
            test_dict["intermediates"][0]["smeared_free_energy"], -8537.247551883
        )
        self.assertEqual(
            test_dict["intermediates"][1]["smeared_free_energy"], -8538.215032441
        )
        self.assertEqual(
            test_dict["intermediates"][-1]["smeared_free_energy"], -8546.922614706
        )
        self.assertEqual(test_dict["geom_iter"], 70)
        self.assertEqual(len(test_dict["intermediates"]), 141)
        self.assertEqual(test_dict["smeared_free_energy"], -8546.922614706)
        self.assertEqual(final_dict["smeared_free_energy"], -8546.922614706)

    def test_castep_parameter_change(self):
        castep_fname = REAL_PATH + "data/castep_files/input-mzs7x1.castep"
        self.assertTrue(
            os.path.isfile(castep_fname),
            msg="Failed to open test case {} - please check installation.".format(
                castep_fname
            ),
        )
        test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
        self.assertTrue(s)
        self.assertTrue(test_dict["optimised"])
        self.assertEqual(test_dict["enthalpy"], -6.16805339e003)
        self.assertEqual(test_dict["total_energy"], -6168.053386094)

    def test_castep_mulliken_scraper(self):
        castep_fname = REAL_PATH + "data/castep_files/Fe-spin.castep"
        self.assertTrue(
            os.path.isfile(castep_fname),
            msg="Failed to open test case {} - please check installation.".format(
                castep_fname
            ),
        )
        test_dict, s = castep2dict(castep_fname, db=False, verbosity=VERBOSITY)
        self.assertTrue(s)
        self.assertEqual(test_dict["task"], "singlepointenergy")
        self.assertEqual(test_dict["atom_types"], ["Fe", "Fe"])
        self.assertAlmostEqual(test_dict["integrated_spin_density"], 4.27207)
        self.assertAlmostEqual(test_dict["integrated_mod_spin_density"], 4.44521)
        self.assertEqual(test_dict["mulliken_spins"], [2.14, 2.14])
        self.assertEqual(test_dict["mulliken_net_spin"], 4.28)
        self.assertEqual(test_dict["mulliken_abs_spin"], 4.28)

    def test_castep_beef_scraper(self):
        from matador.utils.chem_utils import HARTREE_TO_EV

        castep_fname = REAL_PATH + "data/beef_files/K3P_BEEF.castep"
        self.assertTrue(
            os.path.isfile(castep_fname),
            msg="Failed to open test case {} - please check installation".format(
                castep_fname
            ),
        )
        test_dict, s = castep2dict(castep_fname, db=False, verbosity=VERBOSITY)
        self.assertTrue(s)
        self.assertEqual(test_dict["task"], "singlepointenergy")
        self.assertEqual(
            test_dict["atom_types"], ["P", "P", "K", "K", "K", "K", "K", "K"]
        )
        self.assertEqual(len(test_dict["_beef"]["thetas"]), 5000)
        self.assertEqual(len(test_dict["_beef"]["total_energy_per_atom"]), 5000)
        self.assertAlmostEqual(
            test_dict["_beef"]["total_energy"][-1], -1.9029640520e02 * HARTREE_TO_EV
        )
        self.assertAlmostEqual(
            test_dict["_beef"]["total_energy_per_atom"][-1],
            -1.9029640520e02 * HARTREE_TO_EV / 8,
        )
        self.assertAlmostEqual(
            test_dict["_beef"]["mean_total_energy"], -190.6571830577 * HARTREE_TO_EV
        )
        self.assertAlmostEqual(
            test_dict["_beef"]["std_dev_total_energy"], 2.2674151843 * HARTREE_TO_EV
        )
        self.assertAlmostEqual(
            test_dict["_beef"]["mean_total_energy_per_atom"],
            -190.6571830577 * HARTREE_TO_EV / 8,
        )
        self.assertAlmostEqual(
            test_dict["_beef"]["std_dev_total_energy_per_atom"],
            2.2674151843 * HARTREE_TO_EV / 8,
        )

    def test_castep_encap_scraper(self):
        castep_fname = REAL_PATH + "data/encap_files/Se.castep"
        self.assertTrue(
            os.path.isfile(castep_fname),
            msg="Failed to find test case {}, please check installation".format(
                castep_fname
            ),
        )
        test_dict, s = castep2dict(castep_fname, db=False, verbosity=VERBOSITY)
        self.assertTrue(s)
        self.assertEqual(test_dict["task"], "geometryoptimization")
        self.assertEqual(test_dict["atom_types"], 12 * ["Se"])
        self.assertEqual(test_dict["encapsulated"], True)
        self.assertEqual(test_dict["cnt_radius"], 4.69825)
        self.assertEqual(len(test_dict["devel_code"]), 195)

    def test_castep_fixed_cell_scraper(self):
        castep_fname = REAL_PATH + "data/fix_cell_test/TiNb2O7-JVAa6LNI-0K-prim.castep"
        self.assertTrue(
            os.path.isfile(castep_fname),
            msg="Failed to find test case {}, please check installation".format(
                castep_fname
            ),
        )
        test_dict, s = castep2dict(castep_fname, db=False, verbosity=VERBOSITY)
        self.assertTrue(s)
        self.assertEqual(test_dict["positions_frac"][0], [0.501184, 0.501184, 0.997063])
        self.assertEqual(
            test_dict["positions_frac"][-1], [0.182069, 0.182069, 0.989013]
        )
        self.assertEqual(test_dict["lattice_abc"][0], [10.360830, 10.360830, 11.883000])
        self.assertEqual(
            test_dict["lattice_abc"][1], [119.631200, 119.631200, 21.144990]
        )
        self.assertTrue(test_dict["fix_all_cell"])
        self.assertTrue("cell_constraints" not in test_dict)


class ResScraperTests(MatadorUnitTest):
    def test_res(self):
        failed_open = False
        res_fname = REAL_PATH + "data/LiPZn-r57des.res"
        try:
            f = open(res_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    res_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = res2dict(res_fname)
            self.assertTrue(s, "Failed entirely, oh dear!")
            self.assertEqual(
                test_dict["pressure"], 0.0106, msg="Failed to read pressure!"
            )
            self.assertEqual(
                test_dict["enthalpy"], -7600.06148, msg="Failed to read enthalpy!"
            )
            self.assertEqual(test_dict["num_atoms"], 8, msg="Wrong number of atoms!")
            self.assertTrue(
                ["Li", 1] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertTrue(
                ["Zn", 1] in test_dict["stoichiometry"], msg="Wrong stoichiometry!"
            )
            self.assertTrue(
                sorted(test_dict["stoichiometry"]) == test_dict["stoichiometry"],
                msg="Wrong stoichiometry!",
            )
            self.assertEqual(
                test_dict["cell_volume"], 105.918342, msg="Wrong cell volume!"
            )
            self.assertEqual(
                test_dict["space_group"], "Pmc2_1", msg="Wrong space group!"
            )
            self.assertEqual(
                test_dict["lattice_abc"],
                [[5.057429, 4.93404, 4.244619], [90.0, 90.0, 90.0]],
                msg="Wrong lattice constants!",
            )

        res_fname = (
            REAL_PATH
            + "data/hull-NaFeP-afh41_new_Na+Fe+P/FeP2-OQMD_2958-CollCode15027-nospin.res"
        )
        failed_open = False
        try:
            f = open(res_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    res_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = res2dict(res_fname)
            self.assertTrue(s)
            self.assertEqual(test_dict["icsd"], 15027)

        res_fname = REAL_PATH + "data/LiPZn-r57des_bodged.res"
        failed_open = False
        try:
            f = open(res_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    res_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = res2dict(res_fname)
            self.assertFalse(s, "This wasn't meant to succeed!")

    def test_c2x_shelx_res(self):
        res_fname = REAL_PATH + "data/structures/npm.res"
        res, s = res2dict(res_fname, db=False)
        self.assertTrue(s, msg="Failed entirely: {}".format(res))


class ParamScraperTests(MatadorUnitTest):
    """ Test CASTEP param scrapers. """

    def test_param(self):
        param_fname = REAL_PATH + "data/KX.param"
        failed_open = False
        try:
            f = open(param_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    param_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = param2dict(param_fname, db=True)
            self.assertTrue(s, "Failed entirely, oh dear!")
            self.assertEqual(
                test_dict["source"][0].split("/")[-1], "KX.param", msg="Wrong source!"
            )
            self.assertEqual(
                test_dict["task"], "geometryoptimization", msg="Failed to read task!"
            )
            self.assertEqual(
                test_dict["xc_functional"], "PBE", msg="Failed to read xc!"
            )
            self.assertEqual(
                test_dict["perc_extra_bands"], 40.0, msg="Failed to read extra bands!"
            )
            self.assertEqual(
                test_dict["cut_off_energy"], 500, msg="Failed to read cut_off_energy"
            )

            test_dict, s = param2dict(param_fname, db=False)
            self.assertTrue(s, "Failed db=False test entirely, oh dear!")
            self.assertEqual(
                test_dict["source"][0].split("/")[-1],
                "KX.param",
                msg="Wrong db=False source!",
            )
            self.assertEqual(
                test_dict["task"],
                "geometryoptimization",
                msg="Failed to read db=False task!",
            )
            self.assertEqual(
                test_dict["xc_functional"], "PBE", msg="Failed to read db=False xc!"
            )
            self.assertEqual(
                test_dict["fix_occupancy"],
                False,
                msg="Failed to read db=False occupancy!",
            )
            self.assertEqual(
                test_dict["perc_extra_bands"],
                40.0,
                msg="Failed to read db=False extra bands!",
            )
            self.assertEqual(
                test_dict["geom_max_iter"], 200, msg="Wrong db=False geom_max_iter!"
            )
            self.assertEqual(
                test_dict["fixed_npw"], False, msg="Wrong db=False fixed_npw!"
            )
            self.assertEqual(
                test_dict["write_checkpoint"],
                "none",
                msg="Wrong db=False checkpointing!",
            )
            self.assertEqual(
                test_dict["write_cell_structure"],
                True,
                msg="Wrong db=False cell_structure!",
            )

    def test_tricky_param(self):
        param_fname = REAL_PATH + "data/tricky_param.param"
        failed_open = False
        try:
            f = open(param_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    param_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = param2dict(param_fname, db=False, debug=True, verbosity=4)
            self.assertTrue(s, "Failed entirely, oh dear!")
            self.assertEqual(
                test_dict["source"][0].split("/")[-1],
                "tricky_param.param",
                msg="Wrong source!",
            )
            self.assertEqual(
                test_dict["task"],
                "spectral",
                msg="Failed to read non colon delimited field task!",
            )
            self.assertEqual(
                test_dict["perc_extra_bands"], 40.0, msg="Failed to read extra bands!"
            )
            self.assertEqual(
                test_dict["fix_occupancy"],
                True,
                msg="Failed to read lowercase bool fix_occupancy",
            )
            self.assertEqual(
                test_dict["spin_polarized"],
                True,
                msg="Failed to read Anglicised spelling of polarised",
            )
            self.assertFalse("spin_polarised" in test_dict)
            self.assertEqual(
                test_dict["write_cell_structure"],
                True,
                msg="Failed to read = delimited field write_cell_structure",
            )
            self.assertEqual(
                test_dict["cut_off_energy"],
                "50.0 ry",
                msg="Failed to non-eV cut_off_energy.",
            )
            self.assertEqual(
                test_dict["devel_code"],
                "xc_bee: true\nxc_bee_rand_seed: 2\n# including comment\nxc_bee_num_trials: 100\n",
                msg="Failed to read devel code",
            )
            self.assertEqual(len(test_dict), 14)


class ScraperMiscTest(MatadorUnitTest):
    """ Test miscellaneous other scrapers. """

    def test_batch_loading(self):
        """ Test passing a list of files to scraper function, which
        should be handled by decorator.

        """
        castep_fname = []
        castep_fname += [REAL_PATH + "data/castep_files/NaP_intermediates.castep"]
        castep_fname += [
            REAL_PATH + "data/castep_files/Na-edgecase-CollCode10101.castep"
        ]
        castep_fname += [REAL_PATH + "data/castep_files/KP-castep17.castep"]
        castep_fname += [
            REAL_PATH + "data/castep_files/Na3Zn4-swap-ReOs-OQMD_759599.castep"
        ]
        castep_fname += [
            REAL_PATH + "data/castep_files/TiO2_unconverged-MP-10101.castep"
        ]

        cursor, failures = castep2dict(castep_fname, db=True)
        self.assertEqual(len(cursor), 4)
        self.assertEqual(len(failures), 1)

        cursor, failures = castep2dict(
            REAL_PATH + "data/castep_files/*.castep", db=True
        )
        self.assertEqual(len(cursor), 5)
        self.assertEqual(len(failures), 3)

        res_fname = []
        res_fname += [REAL_PATH + "data/LiPZn-r57des.res"]
        res_fname += [REAL_PATH + "data/LiPZn-r57des_bodged.res"]
        cursor, failures = res2dict(res_fname, db=True)
        self.assertEqual(len(cursor), 1)
        self.assertEqual(len(failures), 1)

        res_fname = []
        res_fname += [REAL_PATH + "data/LiPZn-r57des.res"]
        res_fname += [REAL_PATH + "data/LiPZn-r57des_bodged.res"]
        with self.assertRaises(Exception):
            cursor, failures = res2dict(res_fname, db=True, fail_fast=True)

    def test_phonon_scraper(self):
        phonon_fname = REAL_PATH + "data/phonon_dispersion/K3P.phonon"
        self.assertTrue(
            os.path.isfile(phonon_fname),
            msg="Failed to open test case {} - please check installation.".format(
                phonon_fname
            ),
        )
        ph_dict, s = phonon2dict(phonon_fname, verbosity=VERBOSITY)
        self.assertTrue(s, msg="Failed to read phonon file")
        self.assertEqual(ph_dict["num_atoms"], 8)
        self.assertEqual(ph_dict["num_branches"], 24)
        self.assertEqual(ph_dict["num_modes"], 24)
        self.assertEqual(ph_dict["num_kpoints"], 250)
        self.assertEqual(ph_dict["freq_unit"], "cm-1")
        self.assertEqual(ph_dict["lattice_cart"][0], [4.961529, 2.864318, -0.00000])
        self.assertEqual(ph_dict["lattice_cart"][1], [-4.961529, 2.864318, 0.00000])
        self.assertEqual(ph_dict["lattice_cart"][2], [0.000000, 0.000000, 10.127257])
        self.assertEqual(ph_dict["positions_frac"][0], [0.666699, 0.333301, 0.750129])
        self.assertEqual(ph_dict["atom_types"][0], "P")
        self.assertEqual(ph_dict["atom_types"][2], "K")
        self.assertEqual(ph_dict["atom_masses"][0], 30.97376)
        self.assertEqual(ph_dict["atom_masses"][2], 39.0983)
        self.assertEqual(ph_dict["softest_mode_freq"], -23.654487 * INVERSE_CM_TO_EV)

        disp, s = phonon2dict(phonon_fname, verbosity=VERBOSITY, as_model=True)
        self.assertTrue(isinstance(disp, VibrationalDispersion))
        ph_dict["kpoint_branches"] = disp.kpoint_branches
        ph_dict["kpoint_path_spacing"] = disp.kpoint_path_spacing
        self.assertAlmostEqual(ph_dict["kpoint_path_spacing"], 0.021, places=2)
        self.assertEqual(ph_dict["kpoint_branches"][0][0], 0)
        self.assertEqual(ph_dict["kpoint_branches"][0][-1], 35)
        self.assertEqual(ph_dict["kpoint_branches"][1][0], 36)
        self.assertEqual(ph_dict["kpoint_branches"][1][-1], 134)
        self.assertEqual(ph_dict["kpoint_branches"][-2][0], 135)
        self.assertEqual(ph_dict["kpoint_branches"][-1][0], 185)
        self.assertEqual(ph_dict["kpoint_branches"][-1][-1], 249)

    def test_phonon_scraper_ir(self):
        phonon_fname = REAL_PATH + "data/phonon_ir/h-BN_IRR.phonon"
        self.assertTrue(
            os.path.isfile(phonon_fname),
            msg="Failed to open test case {} - please check installation.".format(
                phonon_fname
            ),
        )
        data, s = phonon2dict(phonon_fname, VERBOSITY=VERBOSITY)
        self.assertTrue(s)
        self.assertTrue("infrared_intensity" in data)
        self.assertTrue("raman_intensity" in data)
        self.assertEqual(
            np.shape(data["eigenvalues_q"]), np.shape(data["infrared_intensity"])
        )
        self.assertEqual(
            np.shape(data["eigenvalues_q"]), np.shape(data["raman_intensity"])
        )

    def test_phonon_dos_scraper(self):
        phonon_fname = REAL_PATH + "data/phonon_dispersion/K3P.phonon_dos"
        self.assertTrue(
            os.path.isfile(phonon_fname),
            msg="Failed to open test case {} - please check installation.".format(
                phonon_fname
            ),
        )
        dos_data, s = phonon_dos2dict(phonon_fname)
        self.assertTrue(s)
        self.assertEqual(dos_data["source"], [phonon_fname])
        self.assertEqual(len(dos_data["dos"]), 10001)
        self.assertEqual(len(dos_data["energies"]), 10001)
        self.assertEqual(len(dos_data["pdos"]), 2)
        self.assertEqual(len(dos_data["pdos"][("K", None)]), 10001)
        self.assertEqual(len(dos_data["pdos"][("P", None)]), 10001)

        dos, s = phonon_dos2dict(phonon_fname, as_model=True)
        self.assertTrue(isinstance(dos, VibrationalDOS))

    def test_optados_dos_scraper(self):
        odo_fname = REAL_PATH + "data/optados_files/K3P.adaptive.dat"
        self.assertTrue(
            os.path.isfile(odo_fname),
            msg="Failed to open test case {} - please check installation.".format(
                odo_fname
            ),
        )
        od_dict, s = optados2dict(odo_fname)
        self.assertTrue(s)
        self.assertEqual(len(od_dict["dos"]), 529)
        self.assertEqual(len(od_dict["energies"]), 529)
        self.assertEqual(od_dict["dos_unit_label"], "DOS (electrons per eV/A^3)")
        od, s = optados2dict(odo_fname, as_model=True)
        self.assertTrue(isinstance(od, ElectronicDOS))

    def test_optados_pdos_scraper(self):
        odo_fname = REAL_PATH + "data/optados_files/KP.pdos.adaptive.dat"
        failed_open = False
        try:
            f = open(odo_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    odo_fname
                ),
            )
        if not failed_open:
            f.close()
            od_dict, s = optados2dict(odo_fname)
            self.assertTrue(s)
            self.assertEqual(len(od_dict["sum_pdos"]), 53684)
            self.assertEqual(len(od_dict["energies"]), 53684)
            self.assertEqual(od_dict["num_projectors"], 4)
            self.assertEqual(len(od_dict["pdos"][("K", "s", None)]), 53684)
            self.assertEqual(len(od_dict["pdos"][("K", "p", None)]), 53684)
            self.assertEqual(len(od_dict["pdos"][("P", "s", None)]), 53684)
            self.assertEqual(len(od_dict["pdos"][("P", "p", None)]), 53684)

    def test_optados_spin_pdos_scraper(self):
        odo_fname = REAL_PATH + "data/optados_files/EDASOS-Cr.pdos.dat"
        failed_open = False
        try:
            f = open(odo_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    odo_fname
                ),
            )
        if not failed_open:
            f.close()
            od_dict, s = optados2dict(odo_fname)
            self.assertTrue(s)
            self.assertEqual(len(od_dict["sum_pdos"]), 17366)
            self.assertEqual(len(od_dict["energies"]), 17366)
            self.assertEqual(od_dict["num_projectors"], 10)
            self.assertEqual(len(od_dict["pdos"][("H", None, "up")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("C", None, "up")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("N", None, "up")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("Cl", None, "up")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("Cr", None, "up")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("H", None, "down")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("C", None, "down")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("N", None, "down")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("Cl", None, "down")]), 17366)
            self.assertEqual(len(od_dict["pdos"][("Cr", None, "down")]), 17366)

    def test_optados_pdis_scraper(self):
        odo_fname = REAL_PATH + "data/optados_files/Si2.pdis.dat"
        failed_open = False
        try:
            f = open(odo_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    odo_fname
                ),
            )
        if not failed_open:
            f.close()
            od_dict, s = optados2dict(odo_fname)
            self.assertTrue(s)
            self.assertEqual(len(od_dict["kpoints"]), 166)
            self.assertEqual(od_dict["num_kpoints"], 166)
            self.assertEqual(od_dict["num_bands"], 23)
            self.assertEqual(od_dict["num_projectors"], 4)
            self.assertEqual(np.shape(od_dict["pdis"]), (166, 23, 4))
            self.assertEqual(np.shape(od_dict["eigenvalues"]), (166, 23))
            self.assertEqual(od_dict["projectors"][0], ("Si", "s", None))
            self.assertEqual(od_dict["projectors"][1], ("Si", "p", None))
            self.assertEqual(od_dict["projectors"][2], ("Si", "d", None))
            self.assertEqual(od_dict["projectors"][3], ("Si", "f", None))
            self.assertEqual(od_dict["pdis"][0][0][0], 0.99654675)
            self.assertEqual(od_dict["eigenvalues"][0][0], -12.110537)
            self.assertEqual(od_dict["eigenvalues"][0][-1], 24.862777)
            self.assertEqual(od_dict["eigenvalues"][-1][-1], 24.771165)
            self.assertEqual(od_dict["pdis"][0][0][-1], 0)
            self.assertEqual(od_dict["pdis"][0][-1][1], 0.028667372)
            self.assertEqual(od_dict["pdis"][-1][2][1], 0.99444594)

        odo_fname = REAL_PATH + "data/optados_files/graphite.pdis.dat"
        failed_open = False
        try:
            f = open(odo_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    odo_fname
                ),
            )
        if not failed_open:
            f.close()
            od_dict, s = optados2dict(odo_fname)
            self.assertTrue(s)
            self.assertEqual(len(od_dict["kpoints"]), 942)
            self.assertEqual(od_dict["num_kpoints"], 942)
            self.assertEqual(od_dict["num_bands"], 30)
            self.assertEqual(od_dict["num_projectors"], 4)
            self.assertEqual(np.shape(od_dict["pdis"]), (942, 30, 4))
            self.assertEqual(np.shape(od_dict["eigenvalues"]), (942, 30))
            self.assertEqual(od_dict["projectors"][0], ("C", "s", None))
            self.assertEqual(od_dict["projectors"][1], ("C", "p", None))
            self.assertEqual(od_dict["projectors"][2], ("C", "d", None))
            self.assertEqual(od_dict["projectors"][3], ("C", "f", None))
            self.assertEqual(od_dict["pdis"][29][3][1], 0.85401752)
            self.assertEqual(od_dict["pdis"][30][3][1], 0.84705066)
            self.assertEqual(od_dict["pdis"][31][3][1], 0.84004878)
            self.assertEqual(od_dict["pdis"][32][3][1], 0.83310338)
            self.assertEqual(od_dict["pdis"][33][3][1], 0.82617687)
            self.assertEqual(od_dict["pdis"][34][3][1], 0.81927189)
            self.assertEqual(od_dict["pdis"][35][3][1], 0.81239121)
            self.assertEqual(od_dict["pdis"][36][3][1], 0.80304369)
            self.assertEqual(od_dict["pdis"][37][3][1], 0.79613539)

    def test_arbitrary_scraper(self):
        odi_fname = REAL_PATH + "data/optados_files/testcase.odi"
        od_dict, s = arbitrary2dict(odi_fname)
        self.assertEqual(od_dict["pdispersion"], "species")
        self.assertEqual(od_dict["adaptive_smearing"], "1")
        self.assertEqual(od_dict["set_efermi_zero"], "True")
        self.assertEqual(od_dict["dos_per_volume"], "True")
        self.assertEqual(od_dict["broadening"], "adaptive")
        self.assertEqual(od_dict["dos_spacing"], "0.01")
        self.assertEqual(od_dict["task"], "pdispersion")
        self.assertTrue(od_dict["source"][0].endswith("testcase.odi"))
        self.assertEqual(len(od_dict["source"]), 1)
        self.assertEqual(len(od_dict), 8)

    def test_bands(self):
        from matador.utils.chem_utils import HARTREE_TO_EV

        bands_fname = REAL_PATH + "data/bands_files/KPSn.bands"
        bs_dict, s = bands2dict(bands_fname)
        self.assertTrue(s, msg=bs_dict)
        self.assertEqual(len(bs_dict["kpoint_path"]), 518)
        self.assertEqual(np.shape(bs_dict["eigs_s_k"]), (1, 71, 518))
        self.assertEqual(bs_dict["num_kpoints"], 518)
        self.assertEqual(bs_dict["num_bands"], 71)
        self.assertAlmostEqual(bs_dict["fermi_energy"], 4.0781, places=4)
        self.assertAlmostEqual(bs_dict["spin_fermi_energy"][0], 4.0781, places=4)

        dispersion = ElectronicDispersion(bs_dict)
        self.assertLessEqual(dispersion.kpoint_path_spacing, 0.01)
        self.assertGreaterEqual(dispersion.kpoint_path_spacing, 0.009)
        self.assertAlmostEqual(dispersion.band_gap, 0.760001, places=4)
        self.assertAlmostEqual(dispersion.spin_band_gap[0], 0.760001, places=4)
        self.assertEqual(dispersion.band_gap_path_inds, [246, 235])

        bands_fname = REAL_PATH + "data/bands_files/KPSn_2.bands"
        bs_dict, s = bands2dict(bands_fname)
        self.assertTrue(s)
        self.assertEqual(len(bs_dict["kpoint_path"]), 28)
        self.assertEqual(np.shape(bs_dict["eigs_s_k"]), (1, 71, 28))
        self.assertEqual(bs_dict["num_kpoints"], 28)
        self.assertEqual(bs_dict["num_bands"], 71)

        dispersion = ElectronicDispersion(bs_dict)
        self.assertAlmostEqual(dispersion.fermi_energy, 4.0781, places=4)
        self.assertLessEqual(dispersion.kpoint_path_spacing, 0.3)
        self.assertGreaterEqual(dispersion.kpoint_path_spacing, 0.29)
        self.assertEqual(len(dispersion.kpoint_branches), 2)
        self.assertEqual(
            dispersion["band_gap_path_inds"], dispersion["direct_gap_path_inds"]
        )
        self.assertAlmostEqual(
            bs_dict["eigs_s_k"][0][0][0], -0.99624287 * HARTREE_TO_EV, places=4
        )
        self.assertAlmostEqual(
            bs_dict["eigs_s_k"][-1][-1][-1], 0.74794320 * HARTREE_TO_EV, places=4
        )

        bands_fname = REAL_PATH + "data/bands_files/spin_polarised.bands"
        bs_dict, s = bands2dict(bands_fname)
        self.assertTrue(s)
        self.assertEqual(len(bs_dict["kpoint_path"]), 51)
        self.assertEqual(np.shape(bs_dict["eigs_s_k"]), (2, 462, 51))
        self.assertEqual(bs_dict["num_kpoints"], 51)
        self.assertEqual(bs_dict["num_bands"], 462)
        dispersion = ElectronicDispersion(bs_dict)
        self.assertAlmostEqual(bs_dict["fermi_energy"], 6.7507, places=4)
        self.assertAlmostEqual(dispersion.spin_fermi_energy[0], 6.7507, places=4)
        self.assertAlmostEqual(dispersion.spin_fermi_energy[1], 6.7507, places=4)
        self.assertLessEqual(dispersion.kpoint_path_spacing, 0.03)
        self.assertGreaterEqual(dispersion.kpoint_path_spacing, 0.01)
        self.assertEqual(len(dispersion["kpoint_branches"]), 1)
        self.assertAlmostEqual(
            bs_dict["eigs_s_k"][0][0][0], -1.84888124 * HARTREE_TO_EV, places=4
        )
        self.assertAlmostEqual(
            bs_dict["eigs_s_k"][1][0][0], -1.84666287 * HARTREE_TO_EV, places=4
        )
        self.assertAlmostEqual(
            bs_dict["eigs_s_k"][-1][-1][-1], 0.64283955 * HARTREE_TO_EV, places=4
        )
        self.assertAlmostEqual(
            bs_dict["eigs_s_k"][0][-1][-1], 0.63571135 * HARTREE_TO_EV, places=4
        )
        bs, s = bands2dict(bands_fname, as_model=True)
        self.assertTrue(isinstance(bs, ElectronicDispersion))

    def test_qe_magres(self):
        magres_fname = REAL_PATH + "data/NaP_QE6.magres"
        magres_dict, s = magres2dict(magres_fname)
        self.assertTrue(s)
        self.assertEqual(len(magres_dict["atom_types"]), 4)
        self.assertTrue(
            magres_dict["lattice_cart"][0] == [-2.503686, 2.503686, 3.540961]
        )
        self.assertTrue(
            magres_dict["lattice_cart"][1] == [2.503686, -2.503686, 3.540961]
        )
        self.assertTrue(
            magres_dict["lattice_cart"][2] == [2.503686, 2.503686, -3.540961]
        )

        np.testing.assert_almost_equal(
            magres_dict["susceptibility_tensor"],
            [
                [-2.3100, 0.0000, -0.0000],
                [-0.0000, -2.3100, -0.0000],
                [0.0000, -0.0000, 1.4354],
            ],
        )
        np.testing.assert_almost_equal(
            magres_dict["chemical_shifts"], [518.15, 467.61, 467.61, 275.34], decimal=2
        )

        self.assertEqual(magres_dict["calculator"], "QE-GIPAW")

    def test_castep_magres(self):
        magres_fname = REAL_PATH + "data/LiP_CASTEP18.magres"
        magres_dict, s = magres2dict(magres_fname)
        self.assertTrue(s)
        self.assertEqual(len(magres_dict["atom_types"]), 20)
        self.assertTrue(
            magres_dict["lattice_cart"][0]
            == [4.1332870000000002, 0.0000000000000000, 0.0000000000000000]
        )
        self.assertTrue(
            magres_dict["lattice_cart"][1]
            == [-8.9905292805212659e-4, 6.0637949333506347, 0.0000000000000000]
        )
        self.assertTrue(
            magres_dict["lattice_cart"][2]
            == [2.0677013018922552, 3.3924745014331725e-1, 12.368724395669441]
        )

        np.testing.assert_almost_equal(
            magres_dict["chemical_shifts"],
            [
                83.7,
                84.3,
                83.4,
                86.6,
                83.3,
                85.1,
                84.4,
                83.8,
                82.8,
                83.6,
                84.9,
                84.9,
                83.6,
                82.7,
                85.1,
                350.0,
                500.3,
                353.3,
                530.9,
                531.2,
            ],
            decimal=1,
        )
        np.testing.assert_almost_equal(
            magres_dict["chemical_shift_anisos"],
            [
                9.4,
                4.4,
                8.1,
                2.9,
                8.1,
                3.4,
                4.7,
                9.1,
                10.1,
                -9.5,
                8.7,
                8.8,
                -9.6,
                10.4,
                3.4,
                -393.0,
                162.7,
                -391.2,
                223.9,
                224.0,
            ],
            decimal=1,
        )
        np.testing.assert_almost_equal(
            magres_dict["chemical_shift_asymmetries"],
            [
                0.33,
                0.76,
                0.19,
                0.46,
                0.21,
                0.84,
                0.65,
                0.32,
                0.11,
                0.92,
                0.85,
                0.86,
                0.91,
                0.11,
                0.92,
                0.48,
                0.95,
                0.47,
                0.59,
                0.61,
            ],
            decimal=2,
        )

        self.assertEqual(magres_dict["calculator"], "CASTEP")
        self.assertEqual(magres_dict["calculator_version"], "18.1")

    def test_pwscfout(self):
        pwout_fname = REAL_PATH + "data/NaP.out"
        pwout_dict, s = pwout2dict(pwout_fname)
        self.assertTrue(s)
        self.assertEqual(len(pwout_dict["atom_types"]), 14)
        self.assertEqual(pwout_dict["num_atoms"], 14)
        self.assertTrue(
            pwout_dict["lattice_cart"][0] == [5.887513122, 0.011925355, 0.011971927]
        )
        self.assertTrue(
            pwout_dict["lattice_cart"][1] == [0.605472370, 5.817169640, -0.011329548]
        )
        self.assertTrue(
            pwout_dict["lattice_cart"][2] == [-4.543028478, 0.450282751, 10.044268095]
        )
        self.assertTrue(pwout_dict["source"][0].endswith("NaP.out"))

        self.assertEqual(pwout_dict["pressure"], 0)
        from matador.utils.chem_utils import RY_TO_EV

        np.testing.assert_equal(pwout_dict["enthalpy"], -RY_TO_EV * 97.6314378617)
        np.testing.assert_array_almost_equal(
            pwout_dict["positions_frac"][5], [0.779038368, 0.580790316, 0.631222097]
        )

    def test_usp(self):
        self.assertEqual(
            usp2dict(REAL_PATH + "data/K_OTF.usp")["K"],
            "2|1.5|9|10|11|30U:40:31(qc=6)",
            msg="Failed to scrape K_OTF.usp file",
        )
        self.assertEqual(
            usp2dict(REAL_PATH + "data/P_OTF.usp")["P"],
            "3|1.8|4|4|5|30:31:32",
            msg="Failed to scrape P_OTF.usp file",
        )
        self.assertEqual(
            usp2dict(REAL_PATH + "data/Sn_OTF.usp")["Sn"],
            "2|2|2|1.6|9.6|10.8|11.7|50U=-0.395U=+0.25:51U=-0.14U=+0.25",
            msg="Failed to scrape Sn_OTF.usp file",
        )

    def test_seed_metadata_scrape(self):
        doc = {}
        seed = "blah/blah/blah4/AgBiI4-spinel-Config5-DOI-10.17638__datacat.liverpool.ac.uk__240"
        get_seed_metadata(doc, seed)
        self.assertEqual(doc["doi"], "10.17638/datacat.liverpool.ac.uk/240")
        doc = {}
        seed = "blah/blah/blah4/AgBiI4-spinel-Config5-CollCode123456-from_polish_swaps_garbage"
        get_seed_metadata(doc, seed)
        self.assertEqual(doc["icsd"], 123456)
        doc = {}
        seed = "blah/blah/blah4/AgBiI4-spinel-Config5-CollCode-123456-from_polish_swaps_garbage"
        get_seed_metadata(doc, seed)
        self.assertEqual(doc["icsd"], 123456)
        doc = {}
        seed = "blah/blah/blah4/AgBiI4-spinel-Config5-ICSD-123456-from_polish_swaps_garbage"
        get_seed_metadata(doc, seed)
        self.assertEqual(doc["icsd"], 123456)
        doc = {}
        seed = "blah/blah/blah4/AgBiI4-spinel-Config5-MP-123456-blah-SnPQ"
        get_seed_metadata(doc, seed)
        self.assertEqual(doc["mp_id"], 123456)

    def test_thermo_castep(self):
        castep_fname = REAL_PATH + "data/CuP-thermo-test.castep"
        test_dict, s = castep2dict(castep_fname, db=False, verbosity=VERBOSITY)
        self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
        self.assertEqual(
            test_dict["task"].lower(),
            "thermodynamicscalculation",
            msg="This is not a Thermodynamics calculation...",
        )
        self.assertEqual(
            test_dict["thermo_temp_final"], 1000.0, msg="Wrong final temp!"
        )
        self.assertEqual(test_dict["thermo_temp_init"], 50.0, msg="Wrong initial temp!")
        self.assertEqual(
            test_dict["thermo_temp_spacing"], 100.0, msg="Wrong temp spacing!"
        )
        self.assertEqual(
            test_dict["thermo_num_temp_vals"], 11, msg="Wrong number of temps!"
        )
        self.assertEqual(
            test_dict["thermo_zero_point_energy"],
            0.093412,
            msg="Wrong zero point energy!",
        )

        thermo_db_compare = {
            "thermo_temps": [
                50.0,
                145.0,
                240.0,
                335.0,
                430.0,
                525.0,
                620.0,
                715.0,
                810.0,
                905.0,
                1000.0,
            ],
            "thermo_enthalpy": [
                0.098557,
                0.142535,
                0.204959,
                0.273022,
                0.343308,
                0.414672,
                0.486634,
                0.558962,
                0.63153,
                0.704262,
                0.777113,
            ],
            "thermo_free_energy": [
                0.089968,
                0.050865,
                -0.025747,
                -0.128941,
                -0.252035,
                -0.390909,
                -0.542824,
                -0.705838,
                -0.878507,
                -1.059717,
                -1.248581,
            ],
            "thermo_entropy": [
                16.573,
                60.998,
                92.749,
                115.772,
                133.586,
                148.051,
                160.206,
                170.678,
                179.872,
                188.064,
                195.45,
            ],
            "thermo_heat_cap": [
                24.686,
                57.799,
                67.215,
                70.549,
                72.047,
                72.836,
                73.301,
                73.596,
                73.795,
                73.936,
                74.039,
            ],
        }

        for num, i in enumerate(test_dict["thermo_temps"]):
            self.assertEqual(
                i,
                thermo_db_compare["thermo_temps"][num],
                msg="Wrong temperature %f" % test_dict["thermo_temps"][num],
            )
            self.assertEqual(
                test_dict["thermo_enthalpy"][i],
                thermo_db_compare["thermo_enthalpy"][num],
                msg="Wrong enthalpy %f" % test_dict["thermo_enthalpy"][i],
            )
            self.assertEqual(
                test_dict["thermo_free_energy"][i],
                thermo_db_compare["thermo_free_energy"][num],
                msg="Wrong free energy %f" % test_dict["thermo_free_energy"][i],
            )
            self.assertEqual(
                test_dict["thermo_entropy"][i],
                thermo_db_compare["thermo_entropy"][num],
                msg="Wrong entropy %f" % test_dict["thermo_entropy"][i],
            )
            self.assertEqual(
                test_dict["thermo_heat_cap"][i],
                thermo_db_compare["thermo_heat_cap"][num],
                msg="Wrong heat capacity %f" % test_dict["thermo_heat_cap"][i],
            )

        self.assertEqual(len(test_dict["phonon_fine_kpoint_list"]), 310)
        self.assertEqual(np.shape(test_dict["eigs_q"]), (1, 9, 310))

    def test_fortran_e100_bug(self):
        """ Test whether the scraper handles improperly formatted floats
        by Fortran when e.g. exponent < -99.

        """
        optados_fname = REAL_PATH + "data/fortran_e100_bug/fortran_e100_bug.pdis.dat"
        pdis, s = optados2dict(optados_fname, verbosity=VERBOSITY)
        self.assertTrue(s)

        cell_fname = REAL_PATH + "data/fortran_e100_bug/fortran_e100_bug.cell"
        pdis, s = cell2dict(cell_fname, db=False, lattice=True, verbosity=VERBOSITY)
        self.assertTrue(s)


class CifTests(MatadorUnitTest):
    """ These tests check the cif scraper for correctness. """

    def test_cif_partial_occ(self):
        cif_fname = REAL_PATH + "data/cif_files/AgBiI.cif"

        failed_open = False
        try:
            f = open(cif_fname, "r")
        except Exception:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    cif_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = cif2dict(cif_fname, verbosity=VERBOSITY)
            self.assertTrue(s, "Failed entirely, oh dear! {}".format(test_dict))
            self.assertAlmostEqual(
                test_dict["num_atoms"],
                46.623999999999995,
                msg="Failed to read num_atoms!",
                places=5,
            )
            Bi_ratio = [
                elem[1] for elem in test_dict["stoichiometry"] if elem[0] == "Bi"
            ][0]
            I_ratio = [
                elem[1] for elem in test_dict["stoichiometry"] if elem[0] == "I"
            ][0]
            self.assertEqual(I_ratio / Bi_ratio, 4)
            self.assertAlmostEqual(
                test_dict["cell_volume"],
                1826.0028753,
                msg="Wrong cell volume!",
                places=3,
            )
            self.assertEqual(
                test_dict["space_group"], "Fd-3m", msg="Wrong space group!"
            )
            self.assertEqual(len(test_dict["atom_types"]), 64)
            self.assertEqual(len(test_dict["positions_frac"]), 64)
            self.assertEqual(len(test_dict["site_occupancy"]), 64)

    def test_malicious_cif(self):
        cif_fname = REAL_PATH + "data/cif_files/malicious.cif"
        failed_open = False
        try:
            f = open(cif_fname, "r")
        except FileNotFoundError:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    cif_fname
                ),
            )

        with self.assertRaises(RuntimeError):
            if not failed_open:
                f.close()
            test_dict, s = cif2dict(cif_fname, verbosity=VERBOSITY)
            raise test_dict

    def test_high_symmetry_cif(self):
        cif_fname = REAL_PATH + "data/cif_files/SiO_n001_CollCode1109.cif"
        failed_open = False
        try:
            f = open(cif_fname, "r")
        except FileNotFoundError:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    cif_fname
                ),
            )
        if not failed_open:
            f.close()
            test_dict, s = cif2dict(cif_fname, verbosity=VERBOSITY)
            self.assertTrue(s, "Failed entirely, oh dear! {}".format(test_dict))
            self.assertAlmostEqual(
                test_dict["cell_volume"], 2110.2, msg="Wrong cell volume!", places=1
            )
            self.assertAlmostEqual(
                test_dict["lattice_abc"][0], [18.4940, 4.991, 23.758], places=3
            )
            self.assertAlmostEqual(
                test_dict["lattice_abc"][1], [90, 105.79, 90], places=3
            )
            self.assertEqual(len(test_dict["atom_types"]), 144)
            self.assertEqual(test_dict["num_atoms"], 144)
            self.assertEqual(len(test_dict["positions_frac"]), 144)
            self.assertEqual(len(test_dict["site_occupancy"]), 144)
            self.assertEqual(
                sum(test_dict["site_multiplicity"]), test_dict["num_atoms"]
            )

    def test_problematic_cif(self):
        cif_fname = REAL_PATH + "data/cif_files/SiO_n002_CollCode62404.cif"
        failed_open = False
        try:
            f = open(cif_fname, "r")
        except FileNotFoundError:
            failed_open = True
            self.assertFalse(
                failed_open,
                msg="Failed to open test case {} - please check installation.".format(
                    cif_fname
                ),
            )
        f.close()
        test_dict, s = cif2dict(cif_fname, verbosity=VERBOSITY)
        self.assertTrue(s, "Failed entirely, oh dear! {}".format(test_dict))
        self.assertEqual(sum(test_dict["site_multiplicity"]), test_dict["num_atoms"])
        self.assertEqual(sum(test_dict["site_occupancy"]), test_dict["num_atoms"])
        self.assertEqual(len(test_dict["positions_frac"]), test_dict["num_atoms"])

    def test_big_cif(self):
        cif_fname = REAL_PATH + "data/cif_files/1000001.cif"
        self.assertTrue(
            os.path.isfile(cif_fname),
            msg="Failed to open test case {} - please check installation.".format(
                cif_fname
            ),
        )

        cif, s = cif2dict(cif_fname)
        self.assertTrue(s)
        self.assertEqual(cif['stoichiometry'], [["C", 107], ["H", 142], ["N", 14], ["O", 26]])
        self.assertEqual(cif['space_group'], "P2_12_12_1")
        self.assertAlmostEqual(cif['cell_volume'], 11309.1, places=1)
        self.assertEqual(cif['num_atoms'], 1156)
        self.assertEqual(len(cif['positions_frac']), 1156)


class ExportTest(MatadorUnitTest):
    """ Test file export functions. """

    def test_doc2res(self):
        res_fname = REAL_PATH + "data/LiPZn-r57des.res"
        test_fname = "doc2res.res"
        doc, s = res2dict(res_fname)
        doc2res(doc, test_fname, hash_dupe=False, overwrite=True)
        doc_exported, s = res2dict(test_fname)
        self.assertTrue(s, msg="Failed entirely, oh dear!")
        self.compare_res_with_res(doc, doc_exported)

    def test_doc2param(self):
        param_fname = REAL_PATH + "data/param_test.param"
        test_fname = "dummy.param"
        doc, s = param2dict(param_fname, db=False, verbosity=VERBOSITY)
        self.assertTrue(s, msg="Failed entirely: {}".format(doc))
        doc2param(doc, test_fname, hash_dupe=False, overwrite=True)
        doc_exported, s = param2dict(test_fname, db=False)
        self.assertTrue(s, msg="Failed entirely: {}".format(doc_exported))
        self.assertEqual(len(doc_exported), len(doc))
        self.assertEqual(doc["devel_code"], doc_exported["devel_code"])

        param_fname = REAL_PATH + "data/nmr.param"
        test_fname = "dummy2.param"
        doc, s = param2dict(param_fname, db=False)
        doc2param(doc, test_fname, hash_dupe=False, overwrite=True)
        doc_exported, s = param2dict(test_fname, db=False)
        self.assertTrue(s, msg="Failed entirely, oh dear!")
        self.assertEqual(len(doc_exported), len(doc))

    def test_doc2cell(self):
        cell_fname = REAL_PATH + "data/K5P4-phonon.cell"
        test_fname = "dummy1.cell"

        doc, s = cell2dict(
            cell_fname, db=False, lattice=True, verbosity=VERBOSITY, positions=False
        )
        doc2cell(doc, test_fname)
        test_dict, s = cell2dict(test_fname, db=False, lattice=True, positions=False)
        self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(test_dict))
        self.assertEqual(
            test_dict["lattice_cart"][0][0],
            11.4518745146637,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][1][1],
            5.09448137301246,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][2][2],
            9.18378851243459,
            msg="Failed to read lattice vectors.",
        )
        self.assertEqual(
            test_dict["lattice_cart"][1][0], 0.0, msg="Failed to read lattice vectors."
        )
        self.assertEqual(
            test_dict["symmetry_tol"], 0.0001, msg="Failed to read symmetry tolerance."
        )
        self.assertEqual(
            test_dict["kpoints_mp_spacing"],
            0.03,
            msg="Failed to read kpoint grid {}".format(test_dict["kpoints_mp_spacing"]),
        )
        self.assertEqual(
            test_dict["phonon_kpoint_mp_grid"],
            [2, 2, 2],
            msg="Failed to read kpoint grid {}".format(
                test_dict["phonon_kpoint_mp_grid"]
            ),
        )
        self.assertEqual(
            test_dict["phonon_kpoint_mp_offset"],
            [0.25, 0.25, 0.25],
            msg="Failed to read kpoint grid {}".format(
                test_dict["phonon_kpoint_mp_offset"]
            ),
        )
        self.assertEqual(round(test_dict["phonon_fine_kpoint_mp_spacing"], 2), 0.02)
        self.assertEqual(
            test_dict["species_pot"]["K"],
            "2|1.5|9|10|11|30U:40:31(qc=6)",
            msg="Failed to read pspots.",
        )
        self.assertEqual(
            test_dict["species_pot"]["P"],
            "3|1.8|4|4|5|30:31:32",
            msg="Failed to read pspots.",
        )
        self.assertTrue(test_dict["snap_to_symmetry"])
        self.assertTrue(test_dict["symmetry_generate"])
        self.assertEqual(test_dict["phonon_supercell_matrix"][0], [3, 0, 1])
        self.assertEqual(test_dict["phonon_supercell_matrix"][1], [0, 3, 0])
        self.assertEqual(test_dict["phonon_supercell_matrix"][2], [0, 0, 9])
        self.assertEqual(test_dict["cell_constraints"], [[1, 2, 3], [4, 4, 4]])
        # test that overwrite overwrites
        doc["phonon_supercell_matrix"][2] = [0, 0, 140]
        doc2cell(doc, test_fname, overwrite=True)
        test_dict, s = cell2dict(test_fname, db=False, lattice=True, positions=False)
        self.assertEqual(test_dict["phonon_supercell_matrix"][2], [0, 0, 140])

        # test that hash dupe doesn't
        doc["phonon_supercell_matrix"][2] = [0, 0, 333]
        doc2cell(doc, test_fname, overwrite=False, hash_dupe=True)
        test_dict, s = cell2dict(test_fname, db=False, lattice=True, positions=False)
        self.assertEqual(test_dict["phonon_supercell_matrix"][2], [0, 0, 140])

        dummy_files = glob.glob("dummy*.cell")
        self.assertEqual(len(dummy_files), 2)

    def test_doc2cell_partial_occ_fail(self):
        cell_name = REAL_PATH + "data/cell_files/kpoint_path.cell"
        cell, s = cell2dict(cell_name, db=False)

        cell["site_occupancy"] = []
        for site in cell["positions_frac"]:
            cell["site_occupancy"].append(0.1)

        self.assertTrue(s)
        with self.assertRaises(RuntimeError):
            doc2cell(cell, "dummy")

    def test_doc2cell_kpoint_path(self):
        cell_name = REAL_PATH + "data/cell_files/kpoint_path.cell"
        dummy_name = "dummy2.cell"
        cell, s = cell2dict(cell_name, db=False)

        doc2cell(cell, dummy_name)
        cell, s = cell2dict(dummy_name, db=False)

        self.assertTrue(s)
        self.assertEqual(
            cell["spectral_kpoints_path"],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0], [0.5, 0.0, 0.0]],
        )
        self.assertEqual(
            cell["spectral_kpoints_path_labels"], ["$\\Gamma$", "Z", "$Y$", "X"]
        )
        self.assertEqual(cell["spectral_kpoints_path_spacing"], 0.02)

    def test_doc2cell_and_param_with_spin(self):

        cell_fname = REAL_PATH + "data/K5P4-phonon.cell"
        test_fname = "dummy3"
        test_param = {
            "xc_functional": "PBE",
            "task": "geometryoptimisation",
            "spin_polarized": False,
        }

        doc, s = cell2dict(
            cell_fname, db=False, lattice=True, verbosity=VERBOSITY, positions=True
        )
        doc2cell(doc, test_fname + ".cell", spin=10)

        doc2param(test_param, test_fname + ".param", spin=10)

        param_doc, s = param2dict(test_fname + ".param")
        cell_doc, s = cell2dict(
            test_fname + ".cell", db=False, lattice=True, positions=True
        )

        self.assertTrue(param_doc["spin_polarized"])
        self.assertEqual(param_doc["spin"], 10)

        self.assertEqual(cell_doc["atomic_init_spins"][0], 10)

    def test_doc2res_from_json(self):
        json_fname = REAL_PATH + "data/doc2res.json"
        test_fname = "doc2res.res"
        self.compare_json_with_res(json_fname, test_fname)

    def test_doc2res_from_json_encap(self):
        json_fname = REAL_PATH + "data/doc2res_encap.json"
        test_fname = "doc2res_encap.res"
        self.compare_json_with_res(json_fname, test_fname)

    def test_query2files(self):
        """ Test that MP/ICSD/OQMD structures get written correctly. """
        json_files = glob.glob(REAL_PATH + "data/json_query_files/*.json")
        cursor = []
        for f in json_files:
            with open(f, "r") as _f:
                cursor.append(json.load(_f))
        query2files(cursor, res=True, cell=True)
        self.assertTrue(os.path.isdir("query"))
        self.assertEqual(len(glob.glob("query/*.res")), 3)
        self.assertEqual(len(glob.glob("query/*.cell")), 3)
        fnames = [
            "CuSr-MP_1025402-CollCode629305",
            "H-MP_632250",
            "BaS3Te-OQMD_1606-CollCode8",
        ]
        exts = ["cell", "res"]
        for name, ext in itertools.product(fnames, exts):
            self.assertTrue(
                os.path.isfile("query/{}.{}".format(name, ext)),
                msg="Missing {}.{}".format(name, ext),
            )

    def test_large_writes(self):
        """ Fake some large queries and make sure they are not written. """
        fake_cursor = 100*[{'dummy': 'data'}]
        with self.assertRaises(RuntimeError):
            query2files(fake_cursor, res=True, max_files=99)
        with self.assertRaises(RuntimeError):
            query2files(fake_cursor, res=True, cell=True, max_files=199)
        with self.assertRaises(RuntimeError):
            query2files(fake_cursor, res=True, cell=True, pdb=True, max_files=299)

    def compare_json_with_res(self, json_fname, test_fname):
        with open(json_fname, "r") as f:
            doc = json.load(f)
            doc2res(doc, test_fname, hash_dupe=False, overwrite=True)
            doc_exported, s = res2dict(test_fname)
            self.assertTrue(s, msg="Failed entirely, oh dear!\n{}".format(doc_exported))
            self.compare_res_with_res(doc, doc_exported)

    def compare_res_with_res(self, doc, doc_exported):
        for key in doc_exported:
            if key not in [
                "source",
                "atom_types",
                "positions_frac",
                "stoichiometry",
                "user",
                "lattice_abc",
                "lattice_cart",
                "site_occupancy",
            ]:
                self.assertEqual(
                    doc_exported[key],
                    doc[key],
                    msg="Input and output of {} do not match after scraping.".format(
                        key
                    ),
                )
            elif key == "positions_frac":
                for ind, atom_pos in enumerate(doc_exported["positions_frac"]):
                    self.assertIn(
                        atom_pos,
                        doc["positions_frac"],
                        msg="Atom with this position is missing.",
                    )
                    self.assertEqual(
                        doc_exported["atom_types"][ind],
                        doc["atom_types"][doc["positions_frac"].index(atom_pos)],
                        msg="Atom has wrong type!",
                    )
            elif key == "stoichiometry":
                self.assertEqual(
                    sorted(doc["stoichiometry"]),
                    sorted(doc_exported["stoichiometry"]),
                    msg="Stoichs do not match!",
                )
            elif key == "atom_types":
                self.assertEqual(
                    sorted(doc["atom_types"]),
                    sorted(doc_exported["atom_types"]),
                    msg="Atom types do not match!",
                )
            elif key == "lattice_abc":
                np.testing.assert_almost_equal(
                    doc["lattice_abc"], doc_exported["lattice_abc"]
                )
            elif key == "lattice_cart":
                np.testing.assert_almost_equal(
                    doc["lattice_cart"], doc_exported["lattice_cart"]
                )
