#!/usr/bin/env python
# standard library
import unittest
from os.path import realpath
# matador modules
from matador.crystal import Crystal
from matador.scrapers.castep_scrapers import castep2dict
from matador.scrapers.magres_scrapers import magres2dict
from matador.voronoi_interface import get_voronoi_substructure

# grab abs path for accessing test data
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class CrystalTest(unittest.TestCase):
    def testSites(self):
        doc, s = castep2dict(REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep')
        del doc['lattice_cart']
        crystal = Crystal(doc)
        assert [atom for atom in crystal] == [atom[1] for atom in enumerate(crystal)]

    def testSpg(self):
        doc, s = castep2dict(REAL_PATH + 'data/Na3Zn4-OQMD_759599.castep')
        crystal = Crystal(doc)
        print(crystal.space_group(symprec=0.001))
        print(crystal.space_group(symprec=0.001))

    def testFromMagres(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP.magres')
        crystal = Crystal(doc)
        for atom in crystal:
            print(atom, atom.chemical_shift)

    def testCoordination(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP.magres')
        crystal = Crystal(doc)
        for atom in crystal:
            print(atom, atom.coordination)
        print(crystal.coordination_lists)
        print(crystal.coordination_stats)

    def testVoronoi(self):
        doc, s = magres2dict(REAL_PATH + 'data/NaP.magres')
        crystal = Crystal(doc)
        print(crystal.unique_sites)



if __name__ == '__main__':
    unittest.main(buffer=False, verbosity=2)
