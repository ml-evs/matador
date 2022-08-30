# coding: utf-8
""" This file implements a thin wrapper to the Voronoi code
written by Can Kocer, for use with matador queries and docs.

https://bitbucket.org/can_kocer/ajm_group_voronoi_code

Requires Can's code to be on your PYTHONPATH.
"""

from matador.export import doc2res
from os import remove
from os.path import isfile


def get_voronoi_substructure(doc):
    """Run Can's Voronoi analysis on a matador doc."""
    try:
        from Vornetclass import VoronoiNetwork

        doc2res(
            doc,
            "Vropple.res",
            hash_dupe=False,
            overwrite=True,
            info=False,
            sort_atoms=False,
        )
        vornet = VoronoiNetwork(filename="Vropple.res")
        vornet.computeSubStrucs()
        doc["voronoi_substruc"] = [
            vc.getSubStruc(use_area=False) for vc in vornet.VoronoiCells
        ]
        if isfile("Vropple.res"):
            remove("Vropple.res")
        return doc["voronoi_substruc"]
    except:
        if isfile("Vropple.res"):
            remove("Vropple.res")
        return False


def get_voronoi_points(doc, debug=False):
    """Run Can's Voronoi analysis on a matador doc
    and return nodes, face midpoints and edge midpoints.
    """
    try:
        from Vornetclass import VoronoiNetwork

        doc2res(
            doc,
            "Vropple.res",
            hash_dupe=False,
            overwrite=True,
            info=False,
            sort_atoms=False,
        )
        vornet = VoronoiNetwork(filename="Vropple.res")
        if debug:
            print(vornet.struc)
        vornet.computeVorNet()
        doc["voronoi_nodes"] = vornet.getNodeFracPos()
        doc["voronoi_face_midpoints"] = vornet.getFracFaceMidpoints()
        doc["voronoi_edge_midpoints"] = vornet.getFracEdgeMidpoints()
        if isfile("Vropple.res"):
            remove("Vropple.res")
        return doc["voronoi_nodes"]
    except:
        if isfile("Vropple.res"):
            remove("Vropple.res")
        return False


if __name__ == "__main__":
    from matador.query import DBQuery
    from matador.hull import QueryConvexHull

    # test with LiAs
    query = DBQuery(composition=["LiAs"], subcmd="hull")
    hull = QueryConvexHull(query, no_plot=True, hull_cutoff=0)
    hull_cursor = hull.hull_cursor
    most_lithiated = hull_cursor[-2]
    most_lithiated_substruc = get_voronoi_substructure(most_lithiated)
    other_substruc = []
    for doc in hull_cursor:
        other_substruc.append(get_voronoi_substructure(doc))
    print(other_substruc)
    print(most_lithiated_substruc)
