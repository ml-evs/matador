# coding: utf-8
""" This file implements a thin wrapper to the Voronoi code
written by Can Kocer, for use with matador queries and docs.

https://bitbucket.org/can_kocer/ajm_group_voronoi_code
"""

def doc2voronoi(doc):
    """ Run Can's Voronoi analysis on a matador doc. """
    from Structure import Structure
    # fake init of Structure that normally requires a file
    voro = Structure.__new__(Structure)
    voro.noatoms = doc['num_atoms']
    voro.unitcell = dict()
    voro.unitcell['va'] = doc['lattice_cart'][0]
    voro.unitcell['vb'] = doc['lattice_cart'][1]
    voro.unitcell['vc'] = doc['lattice_cart'][2]
    voro.atomfracpos = doc['positions_frac']
    voro.atomID = dict()
    for ind, atom in enumerate(doc['atom_types']):
        voro.atomID[ind] = atom
    voro.FtCM, voro.CtFM = voro.initConvMat()
    voro.stdizeCell()
    voro.writeFile(fform='cell', fname='temp')

if __name__ == '__main__':
    from matador import DBQuery
    # test with LiAs
    doc = DBQuery(composition=['LiAs'], top=1).cursor[0]
    doc2voronoi(doc)
