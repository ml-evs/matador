#!/usr/bin/env python
import unittest
import json
import numpy as np
from os.path import realpath
from os import system, remove

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class ScrapeTest(unittest.TestCase):
    """ Test scraper functions. """
    def testCellScraper(self):
        from matador.scrapers.castep_scrapers import cell2dict
        cell_fname = REAL_PATH + 'data/LiP2Zn-0bm995-a_9-out.cell'
        failed_open = False
        try:
            f = open(cell_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(cell_fname))

        if not failed_open:
            f.close()
            test_dict, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=0)
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
            test_dict, s = cell2dict(cell_fname, db=False, outcell=False, verbosity=0)
            self.assertTrue(test_dict.get('lattice_cart') is None)

        cell_fname = REAL_PATH + 'data/K5P4-phonon.cell'
        failed_open = False
        try:
            f = open(cell_fname, 'r')
        except:
            failed_open = True
        if not failed_open:
            f.close()
            test_dict, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=0)
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
            self.assertEqual(test_dict['species_pot']['K'], '2|1.5|9|10|11|30U:40:31(qc=6)', msg='Failed to read pspots.')
            self.assertEqual(test_dict['species_pot']['P'], '3|1.8|4|4|5|30:31:32', msg='Failed to read pspots.')
            self.assertTrue(test_dict['snap_to_symmetry'])
            self.assertTrue(test_dict['symmetry_generate'])

    def testCastep16(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=0)
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
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=0)
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
            self.assertEqual(test_dict['estimated_mem_MB'], 300.1)
            self.assertEqual(test_dict['species_pot']['K'], '2|1.5|9|10|11|30U:40:31(qc=6)', msg='Failed to scrape K_OTF.usp file')
            self.assertEqual(test_dict['species_pot']['P'], '3|1.8|4|4|5|30:31:32', msg='Failed to scrape P_OTF.usp file')

    def testCastepSingleAtomEdgeCase(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/Na-edgecase.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=0)
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

            int_dict, s = castep2dict(castep_fname, db=False, intermediates=True, verbosity=10)
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

    def testCastepUnoptimised(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/TiO2_unconverged.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=0)
            self.assertFalse(s, msg='Should have failed with db=True, but didn\'t!')
            self.assertTrue(isinstance(test_dict, str), msg='Should have returned error message!')
            test_dict, s = castep2dict(castep_fname, db=False, verbosity=0)
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

    def testCastepIntermediates(self):
        from matador.scrapers.castep_scrapers import castep2dict
        castep_fname = REAL_PATH + 'data/NaP_intermediates.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(castep_fname))
        if not failed_open:
            f.close()
            fallover = False
            try:
                test_dict, s = castep2dict(castep_fname, db=True, intermediates=True, verbosity=0)
            except:
                fallover = True
            self.assertTrue(fallover, msg='Didn\'t fallover when using db and intermediates')

            test_dict, s = castep2dict(castep_fname, db=False, intermediates=True, verbosity=10)
            self.assertTrue(s, msg='Should have succeeded with db=False, but didn\'t!')
            final_dict, s = castep2dict(castep_fname, db=True, intermediates=False, verbosity=10)
            self.assertTrue(s)
            for key in final_dict:
                self.assertEqual(final_dict[key], test_dict[key], msg='{} didn\'t match'.format(key))
            self.assertEqual(test_dict['intermediates'][0]['total_energy'], -8537.190779552)
            self.assertEqual(test_dict['intermediates'][1]['total_energy'], -8538.161269966)
            self.assertEqual(test_dict['intermediates'][-1]['total_energy'], -8546.922111847)
            self.assertEqual(test_dict['intermediates'][0]['free_energy'], -8537.247551883)
            self.assertEqual(test_dict['intermediates'][1]['free_energy'], -8538.215032441)
            self.assertEqual(test_dict['intermediates'][-1]['free_energy'], -8546.922614706)
            self.assertEqual(len(test_dict['intermediates']), 141)
            self.assertEqual(test_dict['free_energy'], -8546.922614706)
            self.assertEqual(final_dict['free_energy'], -8546.922614706)

    def testRes(self):
        from matador.scrapers.castep_scrapers import res2dict
        res_fname = REAL_PATH + 'data/LiPZn-r57des.res'
        failed_open = False
        try:
            f = open(res_fname, 'r')
        except:
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

    def testParam(self):
        from matador.scrapers.castep_scrapers import param2dict
        param_fname = REAL_PATH + 'data/KX.param'
        failed_open = False
        try:
            f = open(param_fname, 'r')
        except:
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
            self.assertEqual(test_dict['fix_occupancy'], 'false', msg='Failed to read db=False occupancy!')
            self.assertEqual(test_dict['perc_extra_bands'], 40.0, msg='Failed to read db=False extra bands!')
            self.assertEqual(test_dict['geom_max_iter'], '200', msg='Wrong db=False geom_max_iter!')
            self.assertEqual(test_dict['fixed_npw'], 'false', msg='Wrong db=False fixed_npw!')
            self.assertEqual(test_dict['write_checkpoint'], 'none', msg='Wrong db=False checkpointing!')
            self.assertEqual(test_dict['write_cell_structure'], True, msg='Wrong db=False cell_structure!')

    def testBands(self):
        from matador.scrapers.castep_scrapers import bands2dict
        bands_fname = REAL_PATH + 'data/KPSn.bands'
        failed_open = False
        try:
            f = open(bands_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(bands_fname))
        if not failed_open:
            f.close()
            bs_dict, s = bands2dict(bands_fname)
            self.assertEqual(len(bs_dict['kpoint_path']), 518)
            self.assertEqual(np.shape(bs_dict['eigenvalues_k_s']), (1, 71, 518))
            self.assertEqual(bs_dict['num_kpoints'], 518)
            self.assertEqual(bs_dict['num_bands'], 71)
            self.assertAlmostEqual(bs_dict['fermi_energy'], 4.0781, places=4)
            self.assertLessEqual(bs_dict['kpoint_path_spacing'], 0.01)
            self.assertGreaterEqual(bs_dict['kpoint_path_spacing'], 0.009)

        bands_fname = REAL_PATH + 'data/KPSn_2.bands'
        failed_open = False
        try:
            f = open(bands_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(bands_fname))
        if not failed_open:
            f.close()
            bs_dict, s = bands2dict(bands_fname)
            self.assertEqual(len(bs_dict['kpoint_path']), 28)
            self.assertEqual(np.shape(bs_dict['eigenvalues_k_s']), (1, 71, 28))
            self.assertEqual(bs_dict['num_kpoints'], 28)
            self.assertEqual(bs_dict['num_bands'], 71)
            self.assertAlmostEqual(bs_dict['fermi_energy'], 4.0781, places=4)
            self.assertLessEqual(bs_dict['kpoint_path_spacing'], 0.3)
            self.assertGreaterEqual(bs_dict['kpoint_path_spacing'], 0.29)
            self.assertEqual(len(bs_dict['kpoint_branches']), 1)

    def testMagres(self):
        from matador.scrapers.magres_scrapers import magres2dict
        magres_fname = REAL_PATH + 'data/NaP.magres'
        failed_open = False
        try:
            f = open(magres_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(magres_fname))
        if not failed_open:
            f.close()
            magres_dict, s = magres2dict(magres_fname)
            self.assertEqual(len(magres_dict['atom_types']), 4)
            self.assertTrue(magres_dict['lattice_cart'][0] == [-2.503686, 2.503686, 3.540961])
            self.assertTrue(magres_dict['lattice_cart'][1] == [2.503686, -2.503686, 3.540961])
            self.assertTrue(magres_dict['lattice_cart'][2] == [2.503686, 2.503686, -3.540961])

            np.testing.assert_almost_equal(magres_dict['susceptibility_tensor'], [[-2.3100, 0.0000, -0.0000], [-0.0000, -2.3100, -0.0000], [0.0000, -0.0000, 1.4354]])
            np.testing.assert_almost_equal(magres_dict['chemical_shifts'], [518.15, 467.61, 467.61, 275.34], decimal=2)

            self.assertEqual(magres_dict['calculator'], 'QE-GIPAW')

    def testPWSCF(self):
        from matador.scrapers.qe_scrapers import pwout2dict
        pwout_fname = REAL_PATH + 'data/NaP.out'
        failed_open = False
        try:
            f = open(pwout_fname, 'r')
        except:
            failed_open = True
            self.assertFalse(failed_open, msg='Failed to open test case {} - please check installation.'.format(pwout_fname))
        if not failed_open:
            f.close()
            pwout_dict, s = pwout2dict(pwout_fname)
            self.assertEqual(len(pwout_dict['atom_types']), 14)
            self.assertEqual(pwout_dict['num_atoms'], 14)
            self.assertTrue(pwout_dict['lattice_cart'][0] == [5.887513122, 0.011925355, 0.011971927])
            self.assertTrue(pwout_dict['lattice_cart'][1] == [0.605472370, 5.817169640, -0.011329548])
            self.assertTrue(pwout_dict['lattice_cart'][2] == [-4.543028478, 0.450282751, 10.044268095])

            self.assertEqual(pwout_dict['pressure'], 0)
            from matador.utils.chem_utils import RY_TO_EV
            np.testing.assert_equal(pwout_dict['enthalpy'], -RY_TO_EV*97.6314378617)
            np.testing.assert_array_almost_equal(pwout_dict['positions_frac'][5], [0.779038368, 0.580790316, 0.631222097])

    def testUSP(self):
        from matador.scrapers.castep_scrapers import usp2dict
        self.assertEqual(usp2dict(REAL_PATH + 'data/K_OTF.usp')['K'], '2|1.5|9|10|11|30U:40:31(qc=6)', msg='Failed to scrape K_OTF.usp file')
        self.assertEqual(usp2dict(REAL_PATH + 'data/P_OTF.usp')['P'], '3|1.8|4|4|5|30:31:32', msg='Failed to scrape P_OTF.usp file')
        self.assertEqual(usp2dict(REAL_PATH + 'data/Sn_OTF.usp')['Sn'], '2|2|2|1.6|9.6|10.8|11.7|50U=-0.395U=+0.25:51U=-0.14U=+0.25', msg='Failed to scrape Sn_OTF.usp file')


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
        except:
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
        except:
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
        except:
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
        except:
            failed_open = True
        if not failed_open:
            f.close()
            doc, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=0, positions=False)
            doc2cell(doc, test_fname)
            test_dict, s = cell2dict(test_fname, db=False, outcell=True, positions=False)
            remove(test_fname)
            self.assertTrue(s)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
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
        except:
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
            if key not in ['source', 'atom_types', 'positions_frac', 'stoichiometry', 'user', 'lattice_abc', 'lattice_cart']:
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


if __name__ == '__main__':
    unittest.main()
