#!/usr/bin/env python
import unittest
from scrapers.castep_scrapers import castep2dict, res2dict, param2dict, cell2dict


class ScrapeTest(unittest.TestCase):
    """ Test scraper functions. """
    def test(self):
        print('Testing scraping functions...')
        print(70*'-')
        print('Testing cell scraper...')
        cell_fname = 'data/LiP2Zn-0bm995-a_9-out.cell'
        test_cell, s = cell2dict(cell_fname, db=False)
        assert s, 'Cell scraping failed entirely, oh dear!'
        assert test_cell['lattice_cart'][0][0] == 9.83262140721165, 'Cell scraping failed to read lattice vectors.'
        assert test_cell['lattice_cart'][1][1] == 5.96357780025648, 'Cell scraping failed to read lattice vectors.'
        assert test_cell['lattice_cart'][2][2] == 4.39895761828278, 'Cell scraping failed to read lattice vectors.'
        assert test_cell['lattice_cart'][1][0] == -0.115688800302997, 'Cell scraping failed to read lattice vectors.'
        assert test_cell['symmetry_tol'] == 0.001, 'Cell scraping failed to read symmetry tolerance.'
        assert test_cell['kpoints_mp_grid'] == [2, 3, 4], 'Cell scraping failed to read kpoint grid {}'.format(test_cell['kpoints_mp_grid'])
        assert test_cell['species_pot']['Li'] == 'Li_00PBE.usp', 'Cell scraping failed to read pspots.'
        assert test_cell['species_pot']['P'] == 'P_00PBE.usp', 'Cell scraping failed to read pspots.'
        assert test_cell['species_pot']['Zn'] == 'Zn_00PBE.usp', 'Cell scraping failed to read pspots.'

if __name__ == '__main__':
    unittest.main()
