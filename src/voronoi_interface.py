# coding: utf-8
""" This file implements a thin wrapper to the Voronoi code
written by Can Kocer, for use with matador queries and docs.

https://bitbucket.org/can_kocer/ajm_group_voronoi_code
"""
from export import generate_hash


def get_voronoi_substructure(doc):
    """ Run Can's Voronoi analysis on a matador doc. """
    from Structure import Structure
    from Voronizer import Voronizer
    from Vornetclass import VoronoiNetwork
    fname = generate_hash()
    # fake init of Structure that normally requires a file
    voro_struc = Structure.__new__(Structure)
    voro_struc.noatoms = doc['num_atoms']
    voro_struc.unitcell = dict()
    voro_struc.unitcell['va'] = doc['lattice_cart'][0]
    voro_struc.unitcell['vb'] = doc['lattice_cart'][1]
    voro_struc.unitcell['vc'] = doc['lattice_cart'][2]
    voro_struc.atomfracpos = doc['positions_frac']
    voro_struc.atomID = dict()
    for ind, atom in enumerate(doc['atom_types']):
        voro_struc.atomID[ind] = atom
    voro_struc.FtCM, voro_struc.CtFM = voro_struc.initConvMat()
    # voro_struc.stdizeCell()
    voronizer = Voronizer.__new__(Voronizer)
    voronizer.struc = voro_struc
    voronizer.Vornet = VoronoiNetwork(filename=fname, struc=voronizer.struc, fromfile=False)
    voronizer.printInitInfo()
    voronizer.Vornet.computeVorNet()
    doc['substruc'] = [vc.getSubStruc(use_area=True) for vc in voronizer.Vornet.VoronoiCells]
    return doc['substruc']


if __name__ == '__main__':
    from matador import DBQuery, QueryConvexHull
    # test with LiAs
    query = DBQuery(composition=['LiAs'], subcmd='hull')
    hull = QueryConvexHull(query, no_plot=True, hull_cutoff=0)
    hull_cursor = hull.hull_cursor
    most_lithiated = hull_cursor[-2]
    most_lithiated_substruc = get_voronoi_substructure(most_lithiated)
    other_substruc = []
    for doc in hull_cursor:
        other_substruc.append(get_voronoi_substructure(doc))
    print(other_substruc)
    print(most_lithiated_substruc)
