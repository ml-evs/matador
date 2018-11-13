#!/usr/bin/env python
""" Test file scraping and writing functionality. """
import unittest
import json
from os.path import realpath
from os import system, remove

import numpy as np
# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'
VERBOSITY = 0


class ScrapeTest(unittest.TestCase):
    """ Test scraper functions. """
    def testCellScraper(self):
        from matador.scrapers.castep_scrapers import cell2dict
        cell_fname = REAL_PATH + 'data/LiP2Zn-0bm995-a_9-out.cell'
        failed_open = False
        try:
            f = open(cell_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(cell_fname))

        if not failed_open:
            f.close()
            test_dict, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
            self.assertEqual(test_dict['lattice_cart'][0][0], 9.83262140721165, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][1][1], 5.96357780025648, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][2][2], 4.39895761828278, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][1][0], -0.115688800302997, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['symmetry_tol'], 0.001, msg='Failed to read symmetry tolerance.')
            self.assertEqual(test_dict['kpoints_mp_grid'], [2, 3, 4], msg='Failed to read kpoint grid {}'.format(test_dict['kpoints_mp_grid']))
            self.assertEqual(test_dict['species_pot']['Li'], 'Li_00PBE.usp', msg='Failed to read pspots.')
            self.assertEqual(test_dict['species_pot']['P'], 'P_00PBE.usp', msg='Failed to read pspots.')
            self.assertEqual(test_dict['species_pot']['Zn'], 'Zn_00PBE.usp', msg='Failed to read pspots.')
            # test that lattice_vec only read when outcell is true
            test_dict, s = cell2dict(cell_fname, db=False, outcell=False, verbosity=VERBOSITY)
            self.assertTrue(test_dict.get('lattice_cart') is None)

        cell_fname = REAL_PATH + 'data/Li2C2-out.cell'
        failed_open = False
        try:
            f = open(cell_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(cell_fname))

        if not failed_open:
            f.close()
            test_dict, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(test_dict))

        cell_fname = REAL_PATH + 'data/K5P4-phonon.cell'
        failed_open = False
        try:
            f = open(cell_fname, 'r')
        except Exception:
            failed_open = True
        if not failed_open:
            f.close()
            test_dict, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
            self.assertEqual(test_dict['lattice_cart'][0][0], 11.4518745146637, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][1][1], 5.09448137301246, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][2][2], 9.18378851243459, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][1][0], 0.0, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['symmetry_tol'], 0.0001, msg='Failed to read symmetry tolerance.')
            self.assertEqual(test_dict['kpoints_mp_spacing'], 0.03, msg='Failed to read kpoint grid {}'.format(test_dict['kpoints_mp_spacing']))
            self.assertEqual(test_dict['phonon_kpoint_mp_grid'], [2, 2, 2], msg='Failed to read kpoint grid {}'.format(test_dict['phonon_kpoint_mp_grid']))
            self.assertEqual(test_dict['phonon_kpoint_mp_offset'], [0.25, 0.25, 0.25], msg='Failed to read kpoint grid {}'.format(test_dict['phonon_kpoint_mp_offset']))
            self.assertEqual(test_dict['phonon_fine_kpoint_mp_spacing'], 0.02, msg='Failed to read kpoint {}'.format(test_dict['phonon_fine_kpoint_mp_spacing']))
            self.assertEqual(test_dict['phonon_fine_kpoint_path_spacing'], 0.01, msg='Failed to read kpoint {}'.format(test_dict['phonon_fine_kpoint_path_spacing']))
            self.assertEqual(test_dict['species_pot']['K'], '2|1.5|9|10|11|30U:40:31(qc=6)', msg='Failed to read pspots.')
            self.assertEqual(test_dict['species_pot']['P'], '3|1.8|4|4|5|30:31:32', msg='Failed to read pspots.')
            self.assertEqual(test_dict['hubbard_u']['K']['s'], 2, msg='Failed to read Hubbard U block.')
            self.assertEqual(test_dict['hubbard_u']['P']['p'], 3, msg='Failed to read Hubbard U block.')
            self.assertEqual(test_dict['hubbard_u']['U']['d'], 10.101, msg='Failed to read Hubbard U block.')
            self.assertTrue(test_dict['snap_to_symmetry'])
            self.assertTrue(test_dict['symmetry_generate'])
            self.assertEqual(test_dict['phonon_supercell_matrix'][0], [3, 0, 1])
            self.assertEqual(test_dict['phonon_supercell_matrix'][1], [0, 3, 0])
            self.assertEqual(test_dict['phonon_supercell_matrix'][2], [0, 0, 9])

        cell_fname = REAL_PATH + 'data/K5P4-phonon_bodged.cell'
        failed_open = False
        try:
            f = open(cell_fname, 'r')
        except Exception:
            failed_open = True
        if not failed_open:
            f.close()
            test_dict, s = cell2dict(cell_fname, db=True, lattice=True, verbosity=VERBOSITY)
            self.assertFalse(s)

    def testCastep16(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
            self.assertEqual(test_dict['pressure'], 0.0763, msg='Failed to read pressure!')
            self.assertEqual(test_dict['enthalpy'], -2.15036930e4, msg='Failed to read enthalpy!')
            self.assertEqual(test_dict['num_atoms'], 14, msg='Wrong number of atoms!')
            self.assertTrue(['Na', 3] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertTrue(['Zn', 4] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertEqual(test_dict['cell_volume'], 288.041941, msg='Wrong cell volume!')
            self.assertEqual(test_dict['space_group'], 'Pm', msg='Wrong space group!')
            self.assertEqual(test_dict['lattice_abc'][0][0], 9.039776, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][0][1], 9.045651, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][0][2], 4.068682, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][0], 90, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][1], 90, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][2], 59.971185, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['geom_force_tol'], 0.05, msg='Wrong geom force tol')
            self.assertEqual(test_dict['castep_version'], '16.11')
            self.assertEqual(test_dict['estimated_mem_MB'], 345.1)

    def testCastep17(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/KP-castep17.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
            self.assertEqual(test_dict['pressure'], 0.0180, msg='Failed to read pressure!')
            self.assertEqual(test_dict['enthalpy'], -5.98055077e3, msg='Failed to read enthalpy!')
            self.assertEqual(test_dict['num_atoms'], 9, msg='Wrong number of atoms!')
            self.assertTrue(['P', 2] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertTrue(['K', 7] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertEqual(test_dict['cell_volume'], 522.226927, msg='Wrong cell volume!')
            self.assertEqual(test_dict['space_group'], 'Pm', msg='Wrong space group!')
            self.assertEqual(test_dict['lattice_abc'][0][0], 10.231976, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][0][1], 5.024837, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][0][2], 10.186949, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][0], 90.000000, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][1], 94.373377, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][2], 90.000000, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['geom_force_tol'], 0.01, msg='Wrong geom force tol')
            self.assertEqual(test_dict['castep_version'], '17.21')
            self.assertEqual(test_dict['optimised'], True)
            self.assertEqual(test_dict['estimated_mem_MB'], 300.1)
            self.assertEqual(test_dict['species_pot']['K'], '2|1.5|9|10|11|30U:40:31(qc=6)', msg='Failed to scrape K_OTF.usp file')
            self.assertEqual(test_dict['species_pot']['P'], '3|1.8|4|4|5|30:31:32', msg='Failed to scrape P_OTF.usp file')

    def testCastepSingleAtomEdgeCase(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/castep_files/Na-edgecase-CollCode10101.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
            self.assertEqual(test_dict['pressure'], -0.0966, msg='Failed to read pressure!')
            self.assertEqual(test_dict['enthalpy'], -1.30423371e3, msg='Failed to read enthalpy!')
            self.assertEqual(test_dict['positions_frac'], [[0, 0, 0]])
            self.assertEqual(test_dict['forces'], [[0, 0, 0]])
            self.assertEqual(test_dict['enthalpy'], -1.30423371e3, msg='Failed to read enthalpy!')
            self.assertEqual(test_dict['total_energy'], -1304.223019263, msg='Failed to read total energy!')
            self.assertEqual(test_dict['total_energy_per_atom'], -1304.223019263, msg='Failed to read total energy!')
            self.assertEqual(test_dict['free_energy'], -1304.233706274, msg='Failed to read free energy!')
            self.assertEqual(test_dict['num_atoms'], 1, msg='Wrong number of atoms!')
            self.assertTrue(['Na', 1] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertEqual(len(test_dict['stoichiometry']), 1, msg='Wrong stoichiometry!')
            self.assertEqual(test_dict['cell_volume'], 36.761902, msg='Wrong cell volume!')
            self.assertEqual(test_dict['space_group'], 'Im-3m', msg='Wrong space group!')
            self.assertEqual(test_dict['lattice_abc'][0][0], 3.628050, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][0][1], 3.628050, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][0][2], 3.628050, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][0], 109.471221, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][1], 109.471221, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][2], 109.471221, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['geom_force_tol'], 0.05, msg='Wrong geom force tol')
            self.assertEqual(test_dict['castep_version'], '16.1')
            self.assertEqual(test_dict['species_pot']['Na'], 'Na_00PBE.usp')
            self.assertEqual(test_dict['icsd'], 10101)

            int_dict, s = castep2dict(castep_fname, db=False, intermediates=True, verbosity=VERBOSITY)
            for key in test_dict:
                self.assertEqual(test_dict[key], int_dict[key])

            self.assertEqual(len(int_dict['intermediates']), 45)
            for i in range(45):
                self.assertEqual(int_dict['intermediates'][i]['forces'], [[0, 0, 0]])
                self.assertEqual(int_dict['intermediates'][i]['positions_frac'], [[0, 0, 0]])
                self.assertEqual(int_dict['intermediates'][i]['atom_types'], ['Na'])
            self.assertEqual(int_dict['intermediates'][-1]['total_energy'], -1304.223019263)
            self.assertEqual(int_dict['intermediates'][-1]['free_energy'], -1304.233706274)
            self.assertEqual(int_dict['intermediates'][-7]['total_energy'], -1304.222982442)
            self.assertEqual(int_dict['intermediates'][-7]['free_energy'], -1304.233677344)
            self.assertEqual(int_dict['intermediates'][-1]['total_energy_per_atom'], -1304.223019263)
            self.assertEqual(int_dict['intermediates'][-1]['free_energy_per_atom'], -1304.233706274)
            self.assertEqual(int_dict['intermediates'][-7]['total_energy_per_atom'], -1304.222982442)
            self.assertEqual(int_dict['intermediates'][-7]['free_energy_per_atom'], -1304.233677344)
            self.assertEqual(int_dict['geom_iter'], 44)

    def testCastepUnoptimised(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/castep_files/TiO2_unconverged-MP-10101.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertFalse(s, msg='Should have failed with db=True, but didn\'t!')
            self.assertTrue(isinstance(test_dict, str), msg='Should have returned error message!')
            test_dict, s = castep2dict(castep_fname, db=False, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Should have succeeded with db=False, but didn\'t!')
            self.assertTrue(isinstance(test_dict, dict), msg='Should have returned dict!')
            self.assertEqual(test_dict['total_energy'], -12479.86611705)
            self.assertEqual(test_dict['num_atoms'], 12)
            self.assertEqual(test_dict['pressure'], 0.9455, msg='Failed to read pressure!')
            self.assertTrue(['Ti', 1] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertTrue(['O', 2] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertEqual(test_dict['cell_volume'], 127.269750, msg='Wrong cell volume!')
            self.assertEqual(test_dict['space_group'], 'Pmmm')
            self.assertEqual(test_dict['lattice_abc'][0][0], 4.026041, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][0][1], 7.906524, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][0][2], 3.998172, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][0], 90.000000, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][1], 90.000000, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['lattice_abc'][1][2], 90.000000, msg='Wrong lattice constants!')
            self.assertEqual(test_dict['optimised'], False)
            self.assertEqual(test_dict['geom_force_tol'], 0.05)
            self.assertEqual(test_dict['castep_version'], '18.1')
            self.assertEqual(test_dict['species_pot']['Ti'], '3|1.9|8|9|10|30U:40:31:32(qc=5)')
            self.assertEqual(test_dict['species_pot']['O'], '2|1.5|12|13|15|20:21(qc=5)')
            self.assertEqual(test_dict['mp-id'], 10101)

    def testFileNotFound(self):
        """ Ensure that FileNotFound errors fail gracefully. """
        from matador.scrapers.castep_scrapers import res2dict, castep2dict
        error = False
        try:
            res, s = res2dict('___not_a_file')
        except FileNotFoundError:
            error = True
        self.assertTrue(error)

        castep_fname = []
        castep_fname += [REAL_PATH + 'data/castep_files/NaP_intermediates.castep']
        castep_fname += [REAL_PATH + 'data/___not_a_file']
        castep_fname += [REAL_PATH + 'data/KP-castep17.castep']
        castep_fname += [REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep']

        error = False
        try:
            cursor, failures = castep2dict(castep_fname, db=True)
        except FileNotFoundError:
            error = True

    def testBatchLoading(self):
        """ Test passing a list of files to scraper function, which
        should be handled by decorator.

        """
        from matador.scrapers.castep_scrapers import castep2dict, res2dict
        castep_fname = []
        castep_fname += [REAL_PATH + 'data/castep_files/NaP_intermediates.castep']
        castep_fname += [REAL_PATH + 'data/castep_files/Na-edgecase-CollCode10101.castep']
        castep_fname += [REAL_PATH + 'data/castep_files/KP-castep17.castep']
        castep_fname += [REAL_PATH + 'data/castep_files/Na3Zn4-OQMD_759599.castep']
        castep_fname += [REAL_PATH + 'data/castep_files/TiO2_unconverged-MP-10101.castep']

        cursor, failures = castep2dict(castep_fname, db=True)
        self.assertEqual(len(cursor), 4)
        self.assertEqual(len(failures), 1)

        res_fname = []
        res_fname += [REAL_PATH + 'data/LiPZn-r57des.res']
        res_fname += [REAL_PATH + 'data/LiPZn-r57des_bodged.res']
        cursor, failures = res2dict(res_fname, db=True)
        self.assertEqual(len(cursor), 1)
        self.assertEqual(len(failures), 1)

    def testCastepIntermediates(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/castep_files/NaP_intermediates.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=False, intermediates=True, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Should have succeeded with db=False, but didn\'t!')
            final_dict, s = castep2dict(castep_fname, db=True, intermediates=False, verbosity=VERBOSITY)
            self.assertTrue(s)
            for key in final_dict:
                self.assertEqual(final_dict[key], test_dict[key], msg='{} didn\'t match'.format(key))
            self.assertEqual(test_dict['intermediates'][0]['total_energy'], -8537.190779552)
            self.assertEqual(test_dict['intermediates'][1]['total_energy'], -8538.161269966)
            self.assertEqual(test_dict['intermediates'][-1]['total_energy'], -8546.922111847)
            self.assertEqual(test_dict['intermediates'][0]['free_energy'], -8537.247551883)
            self.assertEqual(test_dict['intermediates'][1]['free_energy'], -8538.215032441)
            self.assertEqual(test_dict['intermediates'][-1]['free_energy'], -8546.922614706)
            self.assertEqual(test_dict['geom_iter'], 70)
            self.assertEqual(len(test_dict['intermediates']), 141)
            self.assertEqual(test_dict['free_energy'], -8546.922614706)
            self.assertEqual(final_dict['free_energy'], -8546.922614706)

    def testCastepParameterChange(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/castep_files/input-mzs7x1.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            self.assertTrue(s)
            self.assertTrue(test_dict['optimised'])
            self.assertEqual(test_dict['enthalpy'], -6.16805339E+003)
            self.assertEqual(test_dict['total_energy'], -6168.053386094)

    def testRes(self):
        from matador.scrapers.castep_scrapers import res2dict
        failed_open = False
        res_fname = REAL_PATH + 'data/LiPZn-r57des.res'
        try:
            f = open(res_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(res_fname))
        if not failed_open:
            f.close()
            test_dict, s = res2dict(res_fname)
            self.assertTrue(s, 'Failed entirely, oh dear!')
            self.assertEqual(test_dict['pressure'], 0.0106, msg='Failed to read pressure!')
            self.assertEqual(test_dict['enthalpy'], -7600.06148, msg='Failed to read enthalpy!')
            self.assertEqual(test_dict['num_atoms'], 8, msg='Wrong number of atoms!')
            self.assertTrue(['Li', 1] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertTrue(['Zn', 1] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertTrue(sorted(test_dict['stoichiometry']) == test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertEqual(test_dict['cell_volume'], 105.918342, msg='Wrong cell volume!')
            self.assertEqual(test_dict['space_group'], 'Pmc2_1', msg='Wrong space group!')
            self.assertEqual(test_dict['lattice_abc'], [[5.057429, 4.93404, 4.244619], [90.0, 90.0, 90.0]], msg='Wrong lattice constants!')

        res_fname = REAL_PATH + 'data/hull-NaFeP-afh41_new_Na+Fe+P/FeP2-OQMD_2958-CollCode15027-nospin.res'
        failed_open = False
        try:
            f = open(res_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(res_fname))
        if not failed_open:
            f.close()
            test_dict, s = res2dict(res_fname)
            self.assertTrue(s)
            self.assertEqual(test_dict['icsd'], 15027)

        res_fname = REAL_PATH + 'data/LiPZn-r57des_bodged.res'
        failed_open = False
        try:
            f = open(res_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(res_fname))
        if not failed_open:
            f.close()
            test_dict, s = res2dict(res_fname)
            self.assertFalse(s, 'This wasn\'t meant to succeed!')

    @unittest.skipIf(True, 'CIF tests temporarily disabled...')
    def testCIF(self):
        from matador.scrapers import cif2dict
        cif_fname = REAL_PATH + 'data/cif_files/AgBiI.cif'
        failed_open = False
        try:
            f = open(cif_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(cif_fname))
        if not failed_open:
            f.close()
            test_dict, s = cif2dict(cif_fname)
            self.assertTrue(s, 'Failed entirely, oh dear! {}'.format(test_dict))
            self.assertAlmostEqual(test_dict['num_atoms'], 46.623999999999995, msg='Failed to read num_atoms!', places=5)
            self.assertTrue(['Bi', 1.0] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertTrue(['I', 4.0] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertTrue(sorted(test_dict['stoichiometry']) == test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertAlmostEqual(test_dict['cell_volume'], 1826.0028753, msg='Wrong cell volume!', places=3)
            self.assertEqual(test_dict['space_group'], 'Fd-3m', msg='Wrong space group!')
            self.assertEqual(len(test_dict['atom_types']), 64)
            self.assertEqual(len(test_dict['positions_frac']), 64)
            self.assertEqual(len(test_dict['site_occupancy']), 64)

        cif_fname = REAL_PATH + 'data/cif_files/malicious.cif'
        failed_open = False
        try:
            import os
            f = open(cif_fname, 'r')
        except FileNotFoundError:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(cif_fname))
        if not failed_open:
            f.close()
            errored = False
            test_dict, s = cif2dict(cif_fname)
            errored = isinstance(test_dict, str) and 'RuntimeError' in test_dict
            self.assertTrue(errored, 'WARNING: malicious attack is possible through symops')
            self.assertFalse(s, 'This should have failed entirely, oh dear!')

    def testParam(self):
        from matador.scrapers.castep_scrapers import param2dict
        param_fname = REAL_PATH + 'data/KX.param'
        failed_open = False
        try:
            f = open(param_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(param_fname))
        if not failed_open:
            f.close()
            test_dict, s = param2dict(param_fname, db=True)
            self.assertTrue(s, 'Failed entirely, oh dear!')
            self.assertEqual(test_dict['source'][0].split('/')[-1], 'KX.param', msg='Wrong source!')
            self.assertEqual(test_dict['task'], 'geometryoptimization', msg='Failed to read task!')
            self.assertEqual(test_dict['xc_functional'], 'PBE', msg='Failed to read xc!')
            self.assertEqual(test_dict['perc_extra_bands'], 40.0, msg='Failed to read extra bands!')

            test_dict, s = param2dict(param_fname, db=False)
            self.assertTrue(s, 'Failed db=False test entirely, oh dear!')
            self.assertEqual(test_dict['source'][0].split('/')[-1], 'KX.param', msg='Wrong db=False source!')
            self.assertEqual(test_dict['task'], 'geometryoptimization', msg='Failed to read db=False task!')
            self.assertEqual(test_dict['xc_functional'], 'PBE', msg='Failed to read db=False xc!')
            self.assertEqual(test_dict['fix_occupancy'], False, msg='Failed to read db=False occupancy!')
            self.assertEqual(test_dict['perc_extra_bands'], 40.0, msg='Failed to read db=False extra bands!')
            self.assertEqual(test_dict['geom_max_iter'], 200, msg='Wrong db=False geom_max_iter!')
            self.assertEqual(test_dict['fixed_npw'], False, msg='Wrong db=False fixed_npw!')
            self.assertEqual(test_dict['write_checkpoint'], 'none', msg='Wrong db=False checkpointing!')
            self.assertEqual(test_dict['write_cell_structure'], True, msg='Wrong db=False cell_structure!')

    def testPhononScraper(self):
        from matador.scrapers import phonon2dict
        phonon_fname = REAL_PATH + 'data/K8SnP4.phonon'
        failed_open = False
        try:
            f = open(phonon_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(phonon_fname))
        if not failed_open:
            f.close()
            ph_dict, s = phonon2dict(phonon_fname, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Failed to read phonon file')
            self.assertEqual(ph_dict['num_atoms'], 26)
            self.assertEqual(ph_dict['num_branches'], 78)
            self.assertEqual(ph_dict['num_qpoints'], 110)
            self.assertEqual(ph_dict['freq_unit'], 'cm-1')
            self.assertEqual(ph_dict['lattice_cart'][0], [7.621987, 7.621987, 0.00000])
            self.assertEqual(ph_dict['lattice_cart'][1], [-7.621987, 7.621987, 0.00000])
            self.assertEqual(ph_dict['lattice_cart'][2], [-7.621987, 0.000000, 7.621987])
            self.assertEqual(ph_dict['positions_frac'][0], [0.725087, 0.725075, 0.549843])
            self.assertEqual(ph_dict['atom_types'][0], 'P')
            self.assertEqual(ph_dict['atom_types'][14], 'K')
            self.assertEqual(ph_dict['atom_types'][-1], 'Sn')
            self.assertEqual(ph_dict['atom_masses'][0], 30.97376)
            self.assertEqual(ph_dict['atom_masses'][14], 39.0983)
            self.assertEqual(ph_dict['atom_masses'][-1], 118.710)
            self.assertEqual(ph_dict['softest_mode_freq'], -0.021599)
            self.assertAlmostEqual(ph_dict['qpoint_path_spacing'], 0.021, places=2)
            self.assertEqual(ph_dict['qpoint_branches'][0][0], 0)
            self.assertEqual(ph_dict['qpoint_branches'][0][-1], 28)
            self.assertEqual(ph_dict['qpoint_branches'][1][0], 29)
            self.assertEqual(ph_dict['qpoint_branches'][1][-1], 76)
            self.assertEqual(ph_dict['qpoint_branches'][-1][0], 77)
            self.assertEqual(ph_dict['qpoint_branches'][-1][-1], 109)

    def testOptadosDOSScraper(self):
        from matador.scrapers import optados2dict
        odo_fname = REAL_PATH + 'data/K3P.adaptive.dat'
        failed_open = False
        try:
            f = open(odo_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(odo_fname))
        if not failed_open:
            f.close()
            od_dict, s = optados2dict(odo_fname)
            self.assertTrue(s)
            self.assertEqual(len(od_dict['dos']), 529)
            self.assertEqual(len(od_dict['energies']), 529)

    def testOptadosPDOSScraper(self):
        from matador.scrapers import optados2dict
        odo_fname = REAL_PATH + 'data/KP.pdos.adaptive.dat'
        failed_open = False
        try:
            f = open(odo_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(odo_fname))
        if not failed_open:
            f.close()
            od_dict, s = optados2dict(odo_fname)
            self.assertTrue(s)
            self.assertEqual(len(od_dict['sum_pdos']), 53684)
            self.assertEqual(len(od_dict['energies']), 53684)
            self.assertEqual(od_dict['num_projectors'], 4)
            self.assertEqual(len(od_dict['pdos'][('K', 's')]), 53684)
            self.assertEqual(len(od_dict['pdos'][('K', 'p')]), 53684)
            self.assertEqual(len(od_dict['pdos'][('P', 's')]), 53684)
            self.assertEqual(len(od_dict['pdos'][('P', 'p')]), 53684)

    def testOptadosPDISScraper(self):
        from matador.scrapers import optados2dict
        odo_fname = REAL_PATH + 'data/Si2.pdis.dat'
        failed_open = False
        try:
            f = open(odo_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(odo_fname))
        if not failed_open:
            f.close()
            od_dict, s = optados2dict(odo_fname)
            self.assertTrue(s)
            self.assertEqual(len(od_dict['kpoints']), 166)
            self.assertEqual(od_dict['num_kpoints'], 166)
            self.assertEqual(od_dict['num_bands'], 23)
            self.assertEqual(od_dict['num_projectors'], 4)
            self.assertEqual(np.shape(od_dict['pdis']), (166, 23, 4))
            self.assertEqual(np.shape(od_dict['eigenvalues']), (166, 23))
            self.assertEqual(od_dict['projectors'][0], ('Si', 's'))
            self.assertEqual(od_dict['projectors'][1], ('Si', 'p'))
            self.assertEqual(od_dict['projectors'][2], ('Si', 'd'))
            self.assertEqual(od_dict['projectors'][3], ('Si', 'f'))
            self.assertEqual(od_dict['pdis'][0][0][0], 0.99654675)
            self.assertEqual(od_dict['eigenvalues'][0][0], -12.110537)
            self.assertEqual(od_dict['eigenvalues'][0][-1], 24.862777)
            self.assertEqual(od_dict['eigenvalues'][-1][-1], 24.771165)
            self.assertEqual(od_dict['pdis'][0][0][-1], 0)
            self.assertEqual(od_dict['pdis'][0][-1][1], 0.028667372)
            self.assertEqual(od_dict['pdis'][-1][2][1], 0.99444594)

        odo_fname = REAL_PATH + 'data/graphite.pdis.dat'
        failed_open = False
        try:
            f = open(odo_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(odo_fname))
        if not failed_open:
            f.close()
            od_dict, s = optados2dict(odo_fname)
            self.assertTrue(s)
            self.assertEqual(len(od_dict['kpoints']), 942)
            self.assertEqual(od_dict['num_kpoints'], 942)
            self.assertEqual(od_dict['num_bands'], 30)
            self.assertEqual(od_dict['num_projectors'], 4)
            self.assertEqual(np.shape(od_dict['pdis']), (942, 30, 4))
            self.assertEqual(np.shape(od_dict['eigenvalues']), (942, 30))
            self.assertEqual(od_dict['projectors'][0], ('C', 's'))
            self.assertEqual(od_dict['projectors'][1], ('C', 'p'))
            self.assertEqual(od_dict['projectors'][2], ('C', 'd'))
            self.assertEqual(od_dict['projectors'][3], ('C', 'f'))
            self.assertEqual(od_dict['pdis'][29][3][1], 0.85401752)
            self.assertEqual(od_dict['pdis'][30][3][1], 0.84705066)
            self.assertEqual(od_dict['pdis'][31][3][1], 0.84004878)
            self.assertEqual(od_dict['pdis'][32][3][1], 0.83310338)
            self.assertEqual(od_dict['pdis'][33][3][1], 0.82617687)
            self.assertEqual(od_dict['pdis'][34][3][1], 0.81927189)
            self.assertEqual(od_dict['pdis'][35][3][1], 0.81239121)
            self.assertEqual(od_dict['pdis'][36][3][1], 0.80304369)
            self.assertEqual(od_dict['pdis'][37][3][1], 0.79613539)

    def testBands(self):
        from matador.scrapers.castep_scrapers import bands2dict
        from matador.utils.chem_utils import HARTREE_TO_EV
        bands_fname = REAL_PATH + 'data/bands_files/KPSn.bands'
        failed_open = False
        try:
            f = open(bands_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(bands_fname))
        if not failed_open:
            f.close()
            bs_dict, s = bands2dict(bands_fname, gap=True)
            self.assertTrue(s, msg=bs_dict)
            self.assertEqual(len(bs_dict['kpoint_path']), 518)
            self.assertEqual(np.shape(bs_dict['eigenvalues_k_s']), (1, 71, 518))
            self.assertEqual(bs_dict['num_kpoints'], 518)
            self.assertEqual(bs_dict['num_bands'], 71)
            self.assertAlmostEqual(bs_dict['fermi_energy'], 4.0781, places=4)
            self.assertLessEqual(bs_dict['kpoint_path_spacing'], 0.01)
            self.assertGreaterEqual(bs_dict['kpoint_path_spacing'], 0.009)
            self.assertAlmostEqual(bs_dict['direct_gap'], 0.7807715152197994, places=4)
            self.assertEqual(bs_dict['direct_gap_path_inds'], [0, 0])
            self.assertAlmostEqual(bs_dict['band_gap'], 0.760001, places=4)
            self.assertEqual(bs_dict['band_gap_path_inds'], [246, 235])

        bands_fname = REAL_PATH + 'data/bands_files/KPSn_2.bands'
        failed_open = False
        try:
            f = open(bands_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(bands_fname))
        if not failed_open:
            f.close()
            bs_dict, s = bands2dict(bands_fname, gap=True)
            self.assertTrue(s)
            self.assertEqual(len(bs_dict['kpoint_path']), 28)
            self.assertEqual(np.shape(bs_dict['eigenvalues_k_s']), (1, 71, 28))
            self.assertEqual(bs_dict['num_kpoints'], 28)
            self.assertEqual(bs_dict['num_bands'], 71)
            self.assertAlmostEqual(bs_dict['fermi_energy'], 4.0781, places=4)
            self.assertLessEqual(bs_dict['kpoint_path_spacing'], 0.3)
            self.assertGreaterEqual(bs_dict['kpoint_path_spacing'], 0.29)
            self.assertEqual(len(bs_dict['kpoint_branches']), 2)
            self.assertAlmostEqual(bs_dict['direct_gap'], 0.7807715152197994, places=4)
            self.assertAlmostEqual(bs_dict['direct_gap'], bs_dict['band_gap'], places=4)
            self.assertEqual(bs_dict['direct_gap_path_inds'], [0, 0])
            self.assertEqual(bs_dict['band_gap_path_inds'], bs_dict['direct_gap_path_inds'])
            print((bs_dict['eigenvalues_k_s'][0][0][0]+bs_dict['fermi_energy']) / HARTREE_TO_EV)
            print((bs_dict['eigenvalues_k_s'][0][0][0]-bs_dict['fermi_energy']) / HARTREE_TO_EV)
            self.assertAlmostEqual(bs_dict['eigenvalues_k_s'][0][0][0], -0.99624287*HARTREE_TO_EV-bs_dict['fermi_energy'], places=4)
            self.assertAlmostEqual(bs_dict['eigenvalues_k_s'][-1][-1][-1], 0.74794320*HARTREE_TO_EV-bs_dict['fermi_energy'], places=4)

        bands_fname = REAL_PATH + 'data/bands_files/spin_polarised.bands'
        failed_open = False
        try:
            f = open(bands_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(bands_fname))
        if not failed_open:
            f.close()
            bs_dict, s = bands2dict(bands_fname, gap=True)
            self.assertTrue(s)
            self.assertEqual(len(bs_dict['kpoint_path']), 51)
            self.assertEqual(np.shape(bs_dict['eigenvalues_k_s']), (2, 462, 51))
            self.assertEqual(bs_dict['num_kpoints'], 51)
            self.assertEqual(bs_dict['num_bands'], 462)
            self.assertAlmostEqual(bs_dict['fermi_energy'], 6.7507, places=4)
            self.assertLessEqual(bs_dict['kpoint_path_spacing'], 0.03)
            self.assertGreaterEqual(bs_dict['kpoint_path_spacing'], 0.01)
            self.assertEqual(len(bs_dict['kpoint_branches']), 1)
            self.assertAlmostEqual(bs_dict['eigenvalues_k_s'][0][0][0], -1.84888124*HARTREE_TO_EV-bs_dict['fermi_energy'], places=4)
            self.assertAlmostEqual(bs_dict['eigenvalues_k_s'][1][0][0], -1.84666287*HARTREE_TO_EV-bs_dict['fermi_energy'], places=4)
            self.assertAlmostEqual(bs_dict['eigenvalues_k_s'][-1][-1][-1], 0.64283955*HARTREE_TO_EV-bs_dict['fermi_energy'], places=4)
            self.assertAlmostEqual(bs_dict['eigenvalues_k_s'][0][-1][-1], 0.63571135*HARTREE_TO_EV-bs_dict['fermi_energy'], places=4)

    def testQEMagres(self):
        from matador.scrapers.magres_scrapers import magres2dict
        magres_fname = REAL_PATH + 'data/NaP_QE6.magres'
        failed_open = False
        try:
            f = open(magres_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(magres_fname))
        if not failed_open:
            f.close()
            magres_dict, s = magres2dict(magres_fname)
            self.assertTrue(s)
            self.assertEqual(len(magres_dict['atom_types']), 4)
            self.assertTrue(magres_dict['lattice_cart'][0] == [-2.503686, 2.503686, 3.540961])
            self.assertTrue(magres_dict['lattice_cart'][1] == [2.503686, -2.503686, 3.540961])
            self.assertTrue(magres_dict['lattice_cart'][2] == [2.503686, 2.503686, -3.540961])

            np.testing.assert_almost_equal(magres_dict['susceptibility_tensor'], [[-2.3100, 0.0000, -0.0000], [-0.0000, -2.3100, -0.0000], [0.0000, -0.0000, 1.4354]])
            np.testing.assert_almost_equal(magres_dict['chemical_shifts'], [518.15, 467.61, 467.61, 275.34], decimal=2)

            self.assertEqual(magres_dict['calculator'], 'QE-GIPAW')

    def testCASTEPMagres(self):
        from matador.scrapers.magres_scrapers import magres2dict
        magres_fname = REAL_PATH + 'data/LiP_CASTEP18.magres'
        failed_open = False
        try:
            f = open(magres_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(magres_fname))
        if not failed_open:
            f.close()
            magres_dict, s = magres2dict(magres_fname)
            self.assertTrue(s)
            self.assertEqual(len(magres_dict['atom_types']), 20)
            self.assertTrue(magres_dict['lattice_cart'][0] == [4.1332870000000002, 0.0000000000000000, 0.0000000000000000])
            self.assertTrue(magres_dict['lattice_cart'][1] == [-8.9905292805212659e-4, 6.0637949333506347, 0.0000000000000000])
            self.assertTrue(magres_dict['lattice_cart'][2] == [2.0677013018922552, 3.3924745014331725e-1, 12.368724395669441])

            np.testing.assert_almost_equal(magres_dict['chemical_shifts'],
                                           [83.7, 84.3, 83.4, 86.6, 83.3, 85.1, 84.4, 83.8, 82.8, 83.6, 84.9, 84.9, 83.6, 82.7, 85.1, 350.0, 500.3, 353.3, 530.9, 531.2],
                                           decimal=1)
            np.testing.assert_almost_equal(magres_dict['chemical_shift_anisos'],
                                           [9.4, 4.4, 8.1, 2.9, 8.1, 3.4, 4.7, 9.1, 10.1, -9.5, 8.7, 8.8, -9.6, 10.4, 3.4, -393.0, 162.7, -391.2, 223.9, 224.0],
                                           decimal=1)
            np.testing.assert_almost_equal(magres_dict['chemical_shift_asymmetries'],
                                           [0.33, 0.76, 0.19, 0.46, 0.21, 0.84, 0.65, 0.32, 0.11, 0.92, 0.85, 0.86, 0.91, 0.11, 0.92, 0.48, 0.95, 0.47, 0.59, 0.61],
                                           decimal=2)

            self.assertEqual(magres_dict['calculator'], 'CASTEP')
            self.assertEqual(magres_dict['calculator_version'], '18.1')

    def testPWSCF(self):
        from matador.scrapers.qe_scrapers import pwout2dict
        pwout_fname = REAL_PATH + 'data/NaP.out'
        failed_open = False
        try:
            f = open(pwout_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(pwout_fname))
        if not failed_open:
            f.close()
            pwout_dict, s = pwout2dict(pwout_fname)
            self.assertTrue(s)
            self.assertEqual(len(pwout_dict['atom_types']), 14)
            self.assertEqual(pwout_dict['num_atoms'], 14)
            self.assertTrue(pwout_dict['lattice_cart'][0] == [5.887513122, 0.011925355, 0.011971927])
            self.assertTrue(pwout_dict['lattice_cart'][1] == [0.605472370, 5.817169640, -0.011329548])
            self.assertTrue(pwout_dict['lattice_cart'][2] == [-4.543028478, 0.450282751, 10.044268095])
            self.assertTrue(pwout_dict['source'][0].endswith('NaP.out'))

            self.assertEqual(pwout_dict['pressure'], 0)
            from matador.utils.chem_utils import RY_TO_EV
            np.testing.assert_equal(pwout_dict['enthalpy'], -RY_TO_EV*97.6314378617)
            np.testing.assert_array_almost_equal(pwout_dict['positions_frac'][5], [0.779038368, 0.580790316, 0.631222097])

    def testUSP(self):
        from matador.scrapers.castep_scrapers import usp2dict
        self.assertEqual(usp2dict(REAL_PATH + 'data/K_OTF.usp')['K'], '2|1.5|9|10|11|30U:40:31(qc=6)', msg='Failed to scrape K_OTF.usp file')
        self.assertEqual(usp2dict(REAL_PATH + 'data/P_OTF.usp')['P'], '3|1.8|4|4|5|30:31:32', msg='Failed to scrape P_OTF.usp file')
        self.assertEqual(usp2dict(REAL_PATH + 'data/Sn_OTF.usp')['Sn'], '2|2|2|1.6|9.6|10.8|11.7|50U=-0.395U=+0.25:51U=-0.14U=+0.25', msg='Failed to scrape Sn_OTF.usp file')

    def testSeedMetadataScrape(self):
        from matador.scrapers.castep_scrapers import get_seed_metadata
        doc = {}
        seed = 'blah/blah/blah4/AgBiI4-spinel-Config5-DOI-10.17638__datacat.liverpool.ac.uk__240'
        get_seed_metadata(doc, seed)
        self.assertEqual(doc['doi'], '10.17638/datacat.liverpool.ac.uk/240')
        doc = {}
        seed = 'blah/blah/blah4/AgBiI4-spinel-Config5-CollCode123456-from_polish_swaps_garbage'
        get_seed_metadata(doc, seed)
        self.assertEqual(doc['icsd'], 123456)
        doc = {}
        seed = 'blah/blah/blah4/AgBiI4-spinel-Config5-CollCode-123456-from_polish_swaps_garbage'
        get_seed_metadata(doc, seed)
        self.assertEqual(doc['icsd'], 123456)
        doc = {}
        seed = 'blah/blah/blah4/AgBiI4-spinel-Config5-ICSD-123456-from_polish_swaps_garbage'
        get_seed_metadata(doc, seed)
        self.assertEqual(doc['icsd'], 123456)
        doc = {}
        seed = 'blah/blah/blah4/AgBiI4-spinel-Config5-MP-123456-blah-SnPQ'
        get_seed_metadata(doc, seed)
        self.assertEqual(doc['mp-id'], 123456)


class ExportTest(unittest.TestCase):
    """ Test file export functions. """
    def testDoc2Res(self):
        from matador.scrapers.castep_scrapers import res2dict
        from matador.export import doc2res
        res_fname = REAL_PATH + 'data/LiPZn-r57des.res'
        test_fname = REAL_PATH + 'data/doc2res.res'
        failed_open = False
        try:
            f = open(res_fname, 'r')
        except Exception:
            failed_open = True
            print('Failed to open test case', res_fname, '- please check installation.')
        if not failed_open:
            f.close()
            doc, s = res2dict(res_fname)
            doc2res(doc, test_fname, hash_dupe=False, overwrite=True)
            doc_exported, s = res2dict(test_fname)
            self.assertTrue(s, msg='Failed entirely, oh dear!')
            self.compareResDocwithResDoc(doc, doc_exported)
        system('rm {}'.format(test_fname))

    def testDoc2Param(self):
        from matador.scrapers.castep_scrapers import param2dict
        from matador.export import doc2param
        param_fname = REAL_PATH + 'data/param_test.param'
        test_fname = REAL_PATH + 'data/dummy.param'
        failed_open = False
        try:
            f = open(param_fname, 'r')
        except Exception:
            failed_open = True
            print('Failed to open test case', param_fname, '- please check installation.')
        if not failed_open:
            f.close()
            doc, s = param2dict(param_fname, db=False)
            doc2param(doc, test_fname, hash_dupe=False, overwrite=True)
            doc_exported, s = param2dict(test_fname, db=False)
            remove(test_fname)
            self.assertTrue(s, msg='Failed entirely, oh dear!')
            self.assertEqual(len(doc_exported), len(doc))

        param_fname = REAL_PATH + 'data/nmr.param'
        test_fname = REAL_PATH + 'data/dummy.param'
        failed_open = False
        try:
            f = open(param_fname, 'r')
        except Exception:
            failed_open = True
            print('Failed to open test case', param_fname, '- please check installation.')
        if not failed_open:
            f.close()
            doc, s = param2dict(param_fname, db=False)
            doc2param(doc, test_fname, hash_dupe=False, overwrite=True)
            doc_exported, s = param2dict(test_fname, db=False)
            remove(test_fname)
            self.assertTrue(s, msg='Failed entirely, oh dear!')
            self.assertEqual(len(doc_exported), len(doc))

    def testDoc2Cell(self):
        from matador.scrapers.castep_scrapers import cell2dict
        from matador.export import doc2cell
        cell_fname = REAL_PATH + 'data/K5P4-phonon.cell'
        test_fname = REAL_PATH + 'data/dummy.cell'
        failed_open = False
        try:
            f = open(cell_fname, 'r')
        except Exception:
            failed_open = True
        if not failed_open:
            f.close()
            doc, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=VERBOSITY, positions=False)
            doc2cell(doc, test_fname)
            test_dict, s = cell2dict(test_fname, db=False, outcell=True, positions=False)
            remove(test_fname)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(test_dict))
            self.assertEqual(test_dict['lattice_cart'][0][0], 11.4518745146637, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][1][1], 5.09448137301246, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][2][2], 9.18378851243459, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['lattice_cart'][1][0], 0.0, msg='Failed to read lattice vectors.')
            self.assertEqual(test_dict['symmetry_tol'], 0.0001, msg='Failed to read symmetry tolerance.')
            self.assertEqual(test_dict['kpoints_mp_spacing'], 0.03, msg='Failed to read kpoint grid {}'.format(test_dict['kpoints_mp_spacing']))
            self.assertEqual(test_dict['phonon_kpoint_mp_grid'], [2, 2, 2], msg='Failed to read kpoint grid {}'.format(test_dict['phonon_kpoint_mp_grid']))
            self.assertEqual(test_dict['phonon_kpoint_mp_offset'], [0.25, 0.25, 0.25], msg='Failed to read kpoint grid {}'.format(test_dict['phonon_kpoint_mp_offset']))
            self.assertEqual(round(test_dict['phonon_fine_kpoint_mp_spacing'], 2), 0.02)
            self.assertEqual(test_dict['species_pot']['K'], '2|1.5|9|10|11|30U:40:31(qc=6)', msg='Failed to read pspots.')
            self.assertEqual(test_dict['species_pot']['P'], '3|1.8|4|4|5|30:31:32', msg='Failed to read pspots.')
            self.assertTrue(test_dict['snap_to_symmetry'])
            self.assertTrue(test_dict['symmetry_generate'])
            self.assertEqual(test_dict['phonon_supercell_matrix'][0], [3, 0, 1])
            self.assertEqual(test_dict['phonon_supercell_matrix'][1], [0, 3, 0])
            self.assertEqual(test_dict['phonon_supercell_matrix'][2], [0, 0, 9])

    def testDoc2ResFromJson(self):
        json_fname = REAL_PATH + 'data/doc2res.json'
        test_fname = REAL_PATH + 'data/doc2res.res'
        self.compareJsonWithRes(json_fname, test_fname)

    def testDoc2ResFromJsonWithEncapsulatedStructure(self):
        json_fname = REAL_PATH + 'data/doc2res_encap.json'
        test_fname = REAL_PATH + 'data/doc2res_encap.res'
        self.compareJsonWithRes(json_fname, test_fname)

    def compareJsonWithRes(self, json_fname, test_fname):
        from matador.scrapers.castep_scrapers import res2dict
        from matador.export import doc2res
        failed_open = False
        try:
            f = open(json_fname, 'r')
        except Exception:
            failed_open = True
            print('Failed to open test case', json_fname, '- please check installation.')
            raise AssertionError
        if not failed_open:
            doc = json.load(f)
            f.close()
            doc2res(doc, test_fname, hash_dupe=False, overwrite=True)
            doc_exported, s = res2dict(test_fname)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
            self.compareResDocwithResDoc(doc, doc_exported)
        system('rm {}'.format(test_fname))

    def compareResDocwithResDoc(self, doc, doc_exported):
        for key in doc_exported:
            if key not in ['source', 'atom_types', 'positions_frac', 'stoichiometry', 'user', 'lattice_abc', 'lattice_cart', 'site_occupancy']:
                self.assertEqual(doc_exported[key], doc[key],
                                 msg='Input and output of {} do not match after scraping.'.format(key))
            elif key == 'positions_frac':
                for ind, atom_pos in enumerate(doc_exported['positions_frac']):
                    self.assertIn(atom_pos, doc['positions_frac'], msg='Atom with this position is missing.')
                    self.assertEqual(doc_exported['atom_types'][ind], doc['atom_types'][doc['positions_frac'].index(atom_pos)], msg='Atom has wrong type!')
            elif key == 'stoichiometry':
                self.assertEqual(sorted(doc['stoichiometry']), sorted(doc_exported['stoichiometry']), msg='Stoichs do not match!')
            elif key == 'atom_types':
                self.assertEqual(sorted(doc['atom_types']), sorted(doc_exported['atom_types']), msg='Atom types do not match!')
            elif key == 'lattice_abc':
                np.testing.assert_almost_equal(doc['lattice_abc'], doc_exported['lattice_abc'])
            elif key == 'lattice_cart':
                np.testing.assert_almost_equal(doc['lattice_cart'], doc_exported['lattice_cart'])

    def testThermoCastep(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/CuP-thermo-test.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except Exception:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=False, verbosity=VERBOSITY)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
            self.assertEqual(test_dict['task'].lower(), 'thermodynamicscalculation', msg='This is not a Thermodynamics calculation...')
            self.assertEqual(test_dict['temp_final'], 1000.0, msg='Wrong final temp!')
            self.assertEqual(test_dict['temp_init'], 50.0, msg='Wrong initial temp!')
            self.assertEqual(test_dict['temp_spacing'], 100.0, msg='Wrong temp spacing!')
            self.assertEqual(test_dict['num_temp_vals'], 11, msg='Wrong number of temps!')
            self.assertEqual(test_dict['zero_point_E'], 0.093412, msg='Wrong zero point energy!')

            thermo_db_compare = {'thermo_temps': [50.0, 145.0, 240.0, 335.0, 430.0, 525.0, 620.0, 715.0, 810.0, 905.0, 1000.0],
                                 'thermo_enthalpy': [0.098557, 0.142535, 0.204959, 0.273022, 0.343308, 0.414672, 0.486634, 0.558962, 0.63153, 0.704262, 0.777113],
                                 'thermo_free_energy': [0.089968, 0.050865, -0.025747, -0.128941, -0.252035, -0.390909, -0.542824, -0.705838, -0.878507, -1.059717, -1.248581],
                                 'thermo_entropy': [16.573, 60.998, 92.749, 115.772, 133.586, 148.051, 160.206, 170.678, 179.872, 188.064, 195.45],
                                 'thermo_heat_cap': [24.686, 57.799, 67.215, 70.549, 72.047, 72.836, 73.301, 73.596, 73.795, 73.936, 74.039]}

            for num, i in enumerate(test_dict['thermo_temps']):
                self.assertEqual(i, thermo_db_compare['thermo_temps'][num],
                                 msg='Wrong temperature %f' % test_dict['thermo_temps'][num])
                self.assertEqual(test_dict['thermo_enthalpy'][i], thermo_db_compare['thermo_enthalpy'][num],
                                 msg='Wrong enthalpy %f' % test_dict['thermo_enthalpy'][i])
                self.assertEqual(test_dict['thermo_free_energy'][i], thermo_db_compare['thermo_free_energy'][num],
                                 msg='Wrong free energy %f' % test_dict['thermo_free_energy'][i])
                self.assertEqual(test_dict['thermo_entropy'][i], thermo_db_compare['thermo_entropy'][num],
                                 msg='Wrong entropy %f' % test_dict['thermo_entropy'][i])
                self.assertEqual(test_dict['thermo_heat_cap'][i], thermo_db_compare['thermo_heat_cap'][num],
                                 msg='Wrong heat capacity %f' % test_dict['thermo_heat_cap'][i])


if __name__ == '__main__':
    unittest.main(buffer=True, verbosity=2)
