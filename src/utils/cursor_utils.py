# coding: utf-8
""" This file defines some useful generic cursor methods. """
import numpy as np
from traceback import print_exc


def set_cursor_from_array(cursor, array, key):
    """ Updates the key-value pair for documents in
    internal cursor from a numpy array.
    """
    assert(len(array) == len(cursor) or len(array) - 1 == len(cursor))
    for ind, doc in enumerate(cursor):
        cursor[ind][key] = array[ind]
    return


def get_array_from_cursor(cursor, key):
    """ Returns a numpy array of the values of a key
    in a cursor.
    """
    array = []
    try:
        for doc in cursor:
            array.append(doc[key])
    except:
        print_exc()
    array = np.asarray(array)
    assert(len(array) == len(cursor))
    return array


def get_spg_uniq(cursor, symprec=1e-2, latvecprec=1e-3, posprec=1e-3):
    """ Use spglib to find duplicate structures in a cursor.
    Returns uniq_list and same_list.

    * cursor     : list of matador structure docs.
    * symprec    : spglib symmetry precision for cell standardisation.
    * latvecprec : tolerance on lattice vectors.
    * posprec    : tolerance on fractional atomic positions.
    * uniq_list  : list of indices of the unique structures in cursor.
    * same_list  : list of pairs indices of duplicate structures.
    """

    from utils.cell_utils import doc2spg
    import spglib as spg

    spg_cursor = list()
    for doc in cursor:
        spg_cursor.append(doc2spg(doc))

    refined_list = []
    for crystal in spg_cursor:
        refined_list.append(spg.standardize_cell(crystal, to_primitive=False, no_idealize=False, symprec=symprec))
    for i in range(len(refined_list)):
        for j in range(len(refined_list[i][1])):
            for k in range(len(refined_list[i][1][j])):
                if refined_list[i][1][j][k] > 1-1e-10:
                    refined_list[i][1][j][k] = 0.0
    for i in range(len(refined_list)):
        refined_list[i] = (refined_list[i][0],
                           refined_list[i][1][np.argsort(refined_list[i][1][:, 0])],
                           refined_list[i][2][np.argsort(refined_list[i][1][:, 0])])
    uniq_list = np.arange(0, len(spg_cursor))
    same_list = []
    shift_list = []
    for i in range(len(spg_cursor)):
        for j in range(i+1, len(spg_cursor)):
            if cursor[i]['stoichiometry'] == cursor[j]['stoichiometry']:
                if np.allclose(refined_list[i][0], refined_list[j][0], atol=latvecprec, rtol=0):
                    if np.allclose(refined_list[i][1], refined_list[j][1], atol=posprec, rtol=0):
                        same_list.append((i, j))
                    else:
                        for dim in range(3):
                            if not rigid_shift(refined_list[i], refined_list[j], dim, posprec):
                                    break
                            elif dim == 3:
                                same_list.append((i, j))
                                shift_list.append((i, j, dim))
                            break
    dupes = list(set([pair[1] for pair in same_list]))
    uniq_list = np.delete(uniq_list, dupes)
    print(len(dupes), 'duplicates found and removed.')
    print(len(shift_list), 'of which were shifted cells.')
    return uniq_list


def rigid_shift(structA, structB, dim, posprec):
    assert(len(structA[2]) == len(structB[2]))
    shift_array = structA[1][:, dim] - structB[1][:, dim]
    # if trivial zero shift, return True
    if np.all((np.abs(shift_array)) < 1e-4):
        return True
    shift_array[np.where(shift_array < 0)] += 1
    # if nontrivial shift, return True
    return np.all((np.abs(np.diff(shift_array)) < 1e-4))


def filter_cursor(cursor, key, min, max):
    """ Returns a cursor obeying the filter on the given key. """
    filtered_cursor = list()
    print('Filtering', key, min, max)
    for doc in cursor:
        try:
            if doc[key] < max and doc[key] >= min:
                filtered_cursor.append(doc)
        except:
            pass
    return filtered_cursor
