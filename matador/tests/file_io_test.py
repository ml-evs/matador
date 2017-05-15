#!/usr/bin/env python
import unittest
import json
from matador.scrapers.castep_scrapers import castep2dict, cell2dict, res2dict, param2dict
from matador.export import doc2res
from matador.utils.print_utils import print_warning
from os.path import realpath
from os import system

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class ScrapeTest(unittest.TestCase):
    """ Test scraper functions. """
    def testCellScraper(self):
        cell_fname = REAL_PATH + 'data/LiP2Zn-0bm995-a_9-out.cell'
        failed_open = False
        try:
            f = open(cell_fname, 'r')
        except:
            failed_open = True
            print_warning('Failed to open test case ' + cell_fname + ' - please check installation.')

        if not failed_open:
            f.close()
            test_dict, s = cell2dict(cell_fname, db=False, outcell=True, verbosity=5)
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
            test_dict, s = cell2dict(cell_fname, db=False, outcell=False, verbosity=5)
            self.assertTrue(test_dict.get('lattice_cart') is None)

    def testCastep(self):
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', castep_fname, '- please check installation.')
        if not failed_open:
            f.close()
            test_dict, s = castep2dict(castep_fname, db=True, verbosity=5)
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

    def testRes(self):
        res_fname = REAL_PATH + 'data/LiPZn-r57des.res'
        failed_open = False
        try:
            f = open(res_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', res_fname, '- please check installation.')
        if not failed_open:
            f.close()
            test_dict, s = res2dict(res_fname)
            self.assertTrue(s, 'Failed entirely, oh dear!')
            self.assertEqual(test_dict['pressure'], 0.0106, msg='Failed to read pressure!')
            self.assertEqual(test_dict['enthalpy'], -7600.06148, msg='Failed to read enthalpy!')
            self.assertEqual(test_dict['num_atoms'], 8, msg='Wrong number of atoms!')
            self.assertTrue(['Li', 1] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertTrue(['Zn', 1] in test_dict['stoichiometry'], msg='Wrong stoichiometry!')
            self.assertEqual(test_dict['cell_volume'], 105.918342, msg='Wrong cell volume!')
            self.assertEqual(test_dict['space_group'], 'Pmc2_1', msg='Wrong space group!')
            self.assertEqual(test_dict['lattice_abc'], [[5.057429, 4.93404, 4.244619], [90.0, 90.0, 90.0]], msg='Wrong lattice constants!')

    def testParam(self):
        param_fname = REAL_PATH + 'data/KX.param'
        failed_open = False
        try:
            f = open(param_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', param_fname, '- please check installation.')
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


class ExportTest(unittest.TestCase):
    """ Test file export functions. """
    def testDoc2Res(self):
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
            for key in doc:
                if key not in ['source', 'positions_frac']:
                    self.assertEqual(doc_exported[key], doc[key], msg='Input and output of {} do not match after scraping.'.format(key))
                if key == 'positions_frac':
                    for ind, atom in enumerate(doc['positions_frac']):
                        self.assertIn(atom, doc_exported['positions_frac'], msg='Atom with this position is missing.')
                        self.assertAlmostEqual(doc['atom_types'][ind], doc_exported['atom_types'][doc_exported['positions_frac'].index(atom)], msg='Atom has wrong type!')
            system('rm {}'.format(test_fname))

    def testDoc2ResFromJson(self):
        json_fname = REAL_PATH + 'data/doc2res.json'
        test_fname = REAL_PATH + 'data/doc2res.res'
        failed_open = False
        try:
            f = open(json_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', json_fname, '- please check installation.')
        if not failed_open:
            doc = json.load(f)
            f.close()
            doc2res(doc, test_fname, hash_dupe=False, overwrite=True)
            with open(test_fname, 'r') as f:
                print(f.readlines())
            doc_exported, s = res2dict(test_fname)
            self.assertTrue(s, msg='Failed entirely, oh dear!\n{}'.format(s))
            for key in doc:
                if key == 'positions_frac':
                    for ind, atom in enumerate(doc['positions_frac']):
                        self.assertIn(atom, doc_exported['positions_frac'], msg='Atom with this position is missing.')
                        self.assertEqual(doc['atom_types'][ind], doc_exported['atom_types'][doc_exported['positions_frac'].index(atom)], msg='Atom has wrong type!')
            system('rm {}'.format(test_fname))


if __name__ == '__main__':
    unittest.main()
