# coding: utf-8
""" This file implements a thin wrapper to the Voronoi code
written by Can Kocer, for use with matador queries and docs.

https://bitbucket.org/can_kocer/ajm_group_voronoi_code

Requires Can's code to be on your PYTHONPATH.
"""

from matador.export import doc2res
from os import remove
from Vornetclass import VoronoiNetwork


def get_voronoi_substructure(doc):
    """ Run Can's Voronoi analysis on a matador doc. """
    doc2res(doc, 'Vropple.res', hash_dupe=False)
    vornet = VoronoiNetwork(filename='Vropple.res')
    vornet.computeVorNet()
    doc['substruc'] = [vc.getSubStruc(use_area=False) for vc in vornet.VoronoiCells]
    remove('Vropple.res')
    return doc['substruc']


if __name__ == '__main__':
    from matador.query import DBQuery
    from matador.hull import QueryConvexHull
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
