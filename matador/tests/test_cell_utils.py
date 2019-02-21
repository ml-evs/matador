#!/usr/bin/env python
import unittest
import numpy as np
from os.path import realpath
from matador.utils.cell_utils import abc2cart, cart2abc, cart2volume, create_simple_supercell
from matador.utils.cell_utils import cart2frac, frac2cart
from matador.utils.cell_utils import doc2spg, cart2abcstar, real2recip
from matador.scrapers.castep_scrapers import castep2dict, res2dict, cell2dict, bands2dict
from matador.similarity.pdf_similarity import PDF, PDFOverlap
from matador.export import doc2cell

VERBOSITY = 0

try:
    from matador.utils.cell_utils import get_seekpath_kpoint_path
    from seekpath import get_path
    IMPORTED_SEEKPATH = True
except:
    IMPORTED_SEEKPATH = False

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class CellUtilTest(unittest.TestCase):
    """ Tests cell util functions. """

    def testCart2AbcConversions(self):
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', castep_fname, '- please check installation.')
        if not failed_open:
            f.close()
            test_doc, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            try:
                self.assertTrue(np.allclose(test_doc['lattice_abc'], cart2abc(test_doc['lattice_cart'])),
                                msg='Conversion cart2abc failed.')
                self.assertTrue(np.allclose(cart2abc(test_doc['lattice_cart']), cart2abc(abc2cart(test_doc['lattice_abc']))),
                                msg='Conversion abc2cart failed.')
                self.assertAlmostEqual(test_doc['cell_volume'], cart2volume(test_doc['lattice_cart']),
                                       msg='Failed to calculate volume from lattice vectors.', places=5)
                self.assertIsInstance(test_doc['lattice_abc'], list, msg='Failed abc numpy cast to list')
                self.assertIsInstance(test_doc['lattice_cart'], list, msg='Failed cartesian numpy cast to list')
                cart_pos = frac2cart(test_doc['lattice_cart'], test_doc['positions_frac'])
                back2frac = cart2frac(test_doc['lattice_cart'], cart_pos)
                np.testing.assert_array_almost_equal(back2frac, test_doc['positions_frac'])
            except AssertionError:
                print('cart:', test_doc['lattice_cart'], abc2cart(test_doc['lattice_abc']))
                print('abc:', test_doc['lattice_abc'], cart2abc(test_doc['lattice_cart']))
                print('volume:', test_doc['cell_volume'], cart2volume(test_doc['lattice_cart']))
                raise AssertionError

    def testCart2AbcStar(self):
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        try:
            with open(castep_fname, 'r'):
                pass
        except FileNotFoundError:
            failed_open = True
        self.assertFalse(failed_open,
                         msg='Failed to open test case {} - please check installation'
                         .format(castep_fname))
        test_doc, success = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
        self.assertTrue(success)
        self.assertTrue(np.allclose(real2recip(test_doc['lattice_cart']), 2*np.pi*np.asarray(cart2abcstar(test_doc['lattice_cart']))),
                        msg='Conversion cart2abc failed.')

    def testFrac2Cart(self):
        lattice_cart = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        positions_frac = [0.5, 0.5, 0.5]
        self.assertEqual(frac2cart(lattice_cart, positions_frac), [1, 1, 1])

        lattice_cart = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        positions_frac = [[0.5, 0.5, 0.5]]
        self.assertEqual(frac2cart(lattice_cart, positions_frac), [[1, 1, 1]])

        lattice_cart = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        positions_frac = [[1, 1, 1], [0.5, 0.5, 0.5]]
        self.assertEqual(frac2cart(lattice_cart, positions_frac), [[2, 2, 2], [1, 1, 1]])

    def testCart2Frac(self):
        lattice_cart = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        positions_abs = [1, 1, 1]
        self.assertEqual(cart2frac(lattice_cart, positions_abs), [0.5, 0.5, 0.5])

        lattice_cart = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        positions_abs = [[1, 1, 1]]
        self.assertEqual(cart2frac(lattice_cart, positions_abs), [[0.5, 0.5, 0.5]])

        lattice_cart = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        positions_abs = [[2, 2, 2], [1, 1, 1]]
        self.assertEqual(cart2frac(lattice_cart, positions_abs), [[1, 1, 1], [0.5, 0.5, 0.5]])

    def testConversionTransitivity(self):
        """ Test that cart2frac(frac2cart(A)) == A. """
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        test_doc, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
        lattice_cart = test_doc['lattice_cart']
        positions_frac = test_doc['positions_frac']
        np.testing.assert_almost_equal(cart2frac(lattice_cart, frac2cart(lattice_cart, positions_frac)), positions_frac, decimal=10)

    def testRandomNoise(self):
        """ Test that zero amplitude random noise doesn't do anything, but any other
        amplitude does.
        """
        from matador.utils.cell_utils import add_noise
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        test_doc, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)

        init_positions_frac = test_doc['positions_frac']
        test_doc = add_noise(test_doc, amplitude=0)
        np.testing.assert_almost_equal(init_positions_frac, test_doc['positions_frac'], decimal=10)

        test_doc = add_noise(test_doc, amplitude=0.1)
        np.testing.assert_almost_equal(init_positions_frac, test_doc['positions_frac'], decimal=1)
        diff = np.sqrt(np.sum((np.asarray(init_positions_frac) - np.asarray(test_doc['positions_frac']))**2, axis=-1))
        self.assertTrue(np.max(diff) <= 0.1)
        self.assertTrue(np.max(diff) >= 0)

    def testRecipToReal(self):
        real_lattice = [[5.5902240, 0, 0], [3.7563195, 4.1401290, 0], [-2.9800295, -1.3200288, 8.5321695]]
        recip_lattice = real2recip(real_lattice)
        np.testing.assert_array_almost_equal(np.asarray(recip_lattice),
                                             np.asarray([[1.1239595, -1.0197632, 0.2347956], [0.0, 1.5176303, 0.2347956], [0, 0, 0.7364112]]))

        real_lattice = [[6.0235150, 0, 0], [0.0, 5.6096010, 0], [-5.0202472, 0, 10.0218337]]
        recip_lattice = real2recip(real_lattice)
        np.testing.assert_array_almost_equal(np.asarray(recip_lattice),
                                             np.asarray([[1.0431094, 0, 0.5225256], [0, 1.1200770, 0], [0, 0, 0.6269494]]))

    def testCalcMPSpacing(self):
        from matador.utils.cell_utils import calc_mp_spacing
        real_lattice = [[6.0235150, 0, 0], [0.0, 5.6096010, 0], [-5.0202472, 0, 10.0218337]]
        mp_grid = [4, 4, 2]
        spacing = calc_mp_spacing(real_lattice, mp_grid, prec=3)
        self.assertEqual(spacing, 0.05)
        spacing = calc_mp_spacing(real_lattice, mp_grid, prec=2)
        self.assertAlmostEqual(spacing, 0.05)
        spacing = calc_mp_spacing(real_lattice, mp_grid, prec=5)
        self.assertAlmostEqual(spacing, 0.05, places=3)

    def testCalcMPGrid(self):
        from matador.utils.cell_utils import calc_mp_grid
        real_lattice = [[6.0235150, 0, 0], [0.0, 5.6096010, 0], [-5.0202472, 0, 10.0218337]]
        spacing = 0.05
        self.assertEqual(calc_mp_grid(real_lattice, spacing), [4, 4, 2])

    def testShiftGrid(self):
        from matador.utils.cell_utils import calc_mp_grid, shift_to_include_gamma, abc2cart
        lattice_abc = [[5.57068, 10.222092, 10.222039], [90.0, 90.0, 90.0]]
        lattice_cart = abc2cart(lattice_abc)
        spacing = 0.1
        mp_grid = calc_mp_grid(lattice_cart, spacing)
        self.assertEqual(mp_grid, [2, 1, 1])
        self.assertEqual(shift_to_include_gamma(mp_grid), [0.25, 0, 0])

        lattice_abc = [[5.57068, 10.222092, 10.222039], [90.0, 90.0, 90.0]]
        lattice_cart = abc2cart(lattice_abc)
        spacing = 0.05
        mp_grid = calc_mp_grid(lattice_cart, spacing)
        self.assertEqual(mp_grid, [4, 2, 2])
        self.assertEqual(shift_to_include_gamma(mp_grid), [0.125, 0.25, 0.25])


class SymmetriesAndSupercellsTest(unittest.TestCase):
    """ Tests cell util functions. """

    def testSupercellCreator(self):
        castep_fname = REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep'
        failed_open = False
        num_tests = 3
        try:
            f = open(castep_fname, 'r')
        except:
            failed_open = True
            print('Failed to open test case', castep_fname, '- please check installation.')
        if not failed_open:
            f.close()
            test_doc, s = castep2dict(castep_fname, db=True, verbosity=VERBOSITY)
            _iter = 0
            while _iter < num_tests:
                extension = np.random.randint(low=1, high=5, size=(3)).tolist()
                if extension == [1, 1, 1]:
                    extension[np.random.randint(low=0, high=2)] += 1
                num_images = np.prod(extension)

                standardize = bool(_iter % 2)
                symmetric = bool(_iter % 2)

                supercell = create_simple_supercell(test_doc, tuple(extension),
                                                    standardize=standardize,
                                                    symmetric=symmetric)
                self.assertEqual(supercell['num_atoms'], num_images * test_doc['num_atoms'])
                self.assertAlmostEqual(supercell['cell_volume'], num_images * test_doc['cell_volume'], places=3)
                self.assertEqual(len(supercell['positions_frac']), num_images * len(test_doc['positions_frac']))
                for i in range(3):
                    if not standardize:
                        np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][i]), extension[i]*np.asarray(test_doc['lattice_cart'][i]))
                self.assertLess(pdf_sim_dist(test_doc, supercell), 1e-3)
                _iter += 1

            # test error for 1x1x1
            try:
                supercell = create_simple_supercell(test_doc, (1, 1, 1))
                error = False
            except:
                error = True
            self.assertTrue(error)

        res_fname = REAL_PATH + 'data/parent2.res'
        failed_open = False
        try:
            f = open(res_fname, 'r')
        except:
            failed_open = True
        if not failed_open:
            f.close()
            test_doc, s = res2dict(res_fname, db=False, verbosity=VERBOSITY)
            while _iter < num_tests:
                extension = np.random.randint(low=1, high=5, size=(3, 1)).tolist()
                num_images = np.prod(extension)

                supercell = create_simple_supercell(test_doc, tuple(extension))
                self.assertEqual(supercell['num_atoms'], num_images * test_doc['num_atoms'])
                self.assertAlmostEqual(supercell['cell_volume'], num_images * test_doc['cell_volume'], places=3)
                self.assertEqual(len(supercell['positions_frac']), num_images * len(test_doc['positions_frac']))
                for i in range(3):
                    np.testing.assert_array_equal(np.asarray(supercell['lattice_cart'][i]), extension[i]*np.asarray(test_doc['lattice_cart'][i]))
                self.assertLess(pdf_sim_dist(test_doc, supercell), 1e-3)
                _iter += 1

            # test error for 1x1x1
            try:
                supercell = create_simple_supercell(test_doc, (1, 1, 1))
                error = False
            except:
                error = True
            self.assertTrue(error)

    @unittest.skipIf(not IMPORTED_SEEKPATH, 'Seekpath package not found in this distribution')
    def testKPointPath(self):

        cell, s = castep2dict(REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep')
        std_cell, path, seekpath_results = get_seekpath_kpoint_path(cell, spacing=0.01, debug=False)
        self.assertEqual(539, len(path))

        self.assertLess(pdf_sim_dist(cell, std_cell), 0.05)

        import glob
        from os import remove
        from matador.utils.cell_utils import frac2cart
        fnames = glob.glob(REAL_PATH + 'data/bs_test/*.res')
        spacing = 0.01
        for fname in fnames:
            doc, s = res2dict(fname, db=False)
            doc['cell_volume'] = cart2volume(doc['lattice_cart'])

            std_doc, path, seekpath_results = get_seekpath_kpoint_path(doc, spacing=spacing, debug=False)
            seekpath_results_path = get_path(doc2spg(doc))

            rel_path = seekpath_results['explicit_kpoints_rel']
            abs_path = seekpath_results['explicit_kpoints_abs']

            cart_kpts = np.asarray(frac2cart(real2recip(std_doc['lattice_cart']), path))
            diffs = np.zeros((len(cart_kpts[:-1])))
            np.testing.assert_array_almost_equal(cart_kpts, abs_path)
            np.testing.assert_array_almost_equal(path, rel_path)
            # for ind, kpt in enumerate(cart_kpts[:-1]):
            for ind, kpt in enumerate(cart_kpts[:-1]):
                diffs[ind] = np.sqrt(np.sum((kpt - cart_kpts[ind + 1])**2))
            self.assertLess(len(np.where(diffs > 1.1 * spacing)[0]), len(seekpath_results['explicit_segments']))

            if 'flrys4-1x109' in fname:
                bs, s = bands2dict(fname.replace('.res', '.bands'))
                np.testing.assert_array_almost_equal(bs['kpoint_path'], rel_path)
                np.testing.assert_array_almost_equal(bs['lattice_cart'], std_doc['lattice_cart'])
            self.assertLess(len(np.where(diffs > 1.1 * spacing)[0]), len(seekpath_results['explicit_segments']))

            cell_path = fname.replace('.res', '.cell')
            doc2cell(std_doc, cell_path)
            new_doc, s = cell2dict(cell_path, lattice=True, positions=True, db=False)
            assert 'positions_frac' in new_doc
            remove(cell_path)
            seekpath_new_results = get_path(doc2spg(new_doc))
            self.assertEqual(seekpath_new_results['bravais_lattice_extended'], seekpath_results_path['bravais_lattice_extended'])

            dist = pdf_sim_dist(doc, std_doc)
            self.assertLess(dist, 0.01)
            dist = pdf_sim_dist(doc, new_doc)
            self.assertLess(dist, 0.01)

    def testSpgStandardize(self):
        from matador.utils.cell_utils import standardize_doc_cell
        import glob
        doc, s = castep2dict(REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep')
        std_doc = standardize_doc_cell(doc)
        dist = pdf_sim_dist(doc, std_doc)
        self.assertLess(dist, 0.01)

        fnames = glob.glob(REAL_PATH + 'data/bs_test/*.res')
        for fname in fnames:
            doc, s = res2dict(fname, db=False)
            doc['cell_volume'] = cart2volume(doc['lattice_cart'])
            std_doc = standardize_doc_cell(doc)
            dist = pdf_sim_dist(doc, std_doc)
            self.assertLess(dist, 0.01)


def pdf_sim_dist(doc_test, doc_supercell):
    doc_test['text_id'] = ['test', 'cell']
    doc_supercell['text_id'] = ['super', 'cell']
    pdf_test = PDF(doc_test, low_mem=True)
    pdf_supercell = PDF(doc_supercell, low_mem=True)
    overlap = PDFOverlap(pdf_test, pdf_supercell)
    return overlap.similarity_distance


if __name__ == '__main__':
    unittest.main()
