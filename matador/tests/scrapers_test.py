#!/usr/bin/env python
import unittest
from matador.scrapers.castep_scrapers import castep2dict, cell2dict
from matador.utils.print_utils import print_warning
from traceback import print_exc
from os.path import realpath

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class ScrapeTest(unittest.TestCase):
    """ Test scraper functions. """
    def testCellScraper(self):
        cell_fname = REAL_PATH + 'data/LiP2Zn-0bm995-a_9-out.cell'
        failed_open = False
        failed = False
        try:
            f = open(cell_fname, 'r')
        except:
            failed_open = True
            print_warning('Failed to open test case ' + cell_fname + ' - please check installation.')

        if not failed_open:
            f.close()
            test_dict, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=5)
            try:
                self.assertTrue(s, msg='Failed entirely, oh dear!')
            except Exception as oops:
                print_exc()
                failed = True
                raise AssertionError
            try:
                self.assertEqual(test_dict['lattice_cart'][0][0], 9.83262140721165, msg='Failed to read lattice vectors.')
                self.assertEqual(test_dict['lattice_cart'][1][1], 5.96357780025648, msg='Failed to read lattice vectors.')
                self.assertEqual(test_dict['lattice_cart'][2][2], 4.39895761828278, msg='Failed to read lattice vectors.')
                self.assertEqual(test_dict['lattice_cart'][1][0], -0.115688800302997, msg='Failed to read lattice vectors.')
            except Exception as oops:
                print_exc()
                failed = True
                pass
            try:
                self.assertEqual(test_dict['symmetry_tol'], 0.001, msg='Failed to read symmetry tolerance.')
            except Exception as oops:
                print_exc()
                failed = True
                pass
            try:
                self.assertEqual(test_dict['kpoints_mp_grid'], [2, 3, 4], msg='Failed to read kpoint grid {}'.format(test_dict['kpoints_mp_grid']))
            except Exception as oops:
                print_exc()
                failed = True
                pass
            try:
                self.assertEqual(test_dict['species_pot']['Li'], 'Li_00PBE.usp', msg='Failed to read pspots.')
                self.assertEqual(test_dict['species_pot']['P'], 'P_00PBE.usp', msg='Failed to read pspots.')
                self.assertEqual(test_dict['species_pot']['Zn'], 'Zn_00PBE.usp', msg='Failed to read pspots.')
            except Exception as oops:
                print_exc()
                failed = True
                pass
            try:
                # test that lattice_vec only read when outcell is true
                test_dict, s = cell2dict(cell_fname, db=False, outcell=False, verbosity=5)
                self.assertTrue(test_dict.get('lattice_cart') is None)
            except:
                print_exc()
                failed = True
            if failed:
                raise(AssertionError, 'Cell test failed!')

    def testCastep(self):
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        failed = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', castep_fname, '- please check installation.')
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=5)
            try:
                self.assertTrue(s, 'Failed entirely, oh dear!')
            except Exception as oops:
                print_exc()
                failed = True
                raise AssertionError
            try:
                self.assertEqual(test_dict['pressure'], 0.0763, msg='Failed to read pressure!')
            except Exception as oops:
                print_exc()
                failed = True
            try:
                self.assertEqual(test_dict['enthalpy'], -2.15036930e4, msg='Failed to read enthalpy!')
            except Exception as oops:
                print_exc()
                failed = True
            try:
                self.assertEqual(test_dict['num_atoms'], 14, msg='Wrong number of atoms!')
            except Exception as oops:
                print_exc()
                failed = True
            try:
                self.assertTrue(['Na', 3] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
                self.assertTrue(['Zn', 4] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            except Exception as oops:
                print_exc()
                failed = True
            try:
                self.assertEqual(test_dict['cell_volume'], 288.041941, msg='Wrong cell volume!')
            except Exception as oops:
                print_exc()
                failed = True
            try:
                self.assertEqual(test_dict['space_group'], 'Pm', msg='Wrong space group!')
            except Exception as oops:
                print_exc()
                failed = True
            try:
                self.assertEqual(test_dict['lattice_abc'][0][0], 9.039776, msg='Wrong lattice constants!')
                self.assertEqual(test_dict['lattice_abc'][0][1], 9.045651, msg='Wrong lattice constants!')
                self.assertEqual(test_dict['lattice_abc'][0][2], 4.068682, msg='Wrong lattice constants!')
                self.assertEqual(test_dict['lattice_abc'][1][0], 90, msg='Wrong lattice constants!')
                self.assertEqual(test_dict['lattice_abc'][1][1], 90, msg='Wrong lattice constants!')
                self.assertEqual(test_dict['lattice_abc'][1][2], 59.971185, msg='Wrong lattice constants!')
            except Exception as oops:
                print_exc()
                failed = True

        if failed:
            raise(AssertionError, 'Castep test failed!')

if __name__ == '__main__':
    unittest.main()
