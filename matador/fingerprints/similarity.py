# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements filtering based on Fingerprint objects,
although only PDF has been implemented so far.

"""

import copy
from collections import defaultdict
import numpy as np
from matador.fingerprints.pdf import PDF, PDFFactory
from matador.utils.cursor_utils import get_guess_doc_provenance


def get_uniq_cursor(cursor, sim_tol=0.1, energy_tol=1e-2,
                    enforce_same_stoich=True, fingerprint=PDF,
                    hierarchy_order=None, hierarchy_values=None,
                    debug=False, **fingerprint_calc_args):
    """ Uses fingerprint to filter cursor into unique structures to some
    tolerance sim_tol, additionally returning a dict of duplicates and the
    correlation matrix.

    The choice of which of the dulpicates is kept in the unique cursor is
    defined by the "hierarchy". By default, this will guess the provenance
    of a document and prefer structures from "primary sources", i.e.
    ICSD -> OQMD -> Materials Project -> SWAPS -> AIRSS -> GA. A custom hiearchy
    can be provided through `hierarchy_order`, which must be accompanied by a list
    of values per structure to check against that hierarchy.

    Parameters:
        cursor (list) : matador cursor to be filtered

    Keyword Arguments:
        fingerprint (Fingerprint): fingerprint object type to compare
            (DEFAULT: PDF)
        sim_tol (float/bool): tolerance in similarity distance for
            duplicates (if True, default value of 0.1 used)
        energy_tol (float): compare only structures within a certain
            energy tolerance (1e20 if enforce_same_stoich is False)
        enforce_same_stoich (bool): compare only structures of the same
            stoichiometry
        debug (bool): print timings and list similarities
        fingerprint_calc_args (dict): kwargs to pass to fingerprint

    Returns:
        distinct_set (list): ordered list indices of unique documents
        dupe_dict (dict): a dict with keys from distinct_set, listing duplicates
        fingerprint_list (list): a list of <Fingerprint> objects
        sim_mat (np.ndarray): the correlation matrix of pair similarity distances

    """

    if isinstance(sim_tol, bool):
        sim_tol = 0.1

    if not cursor:
        raise RuntimeError("No structures provided to compare.")

    fingerprint_list = []
    if not enforce_same_stoich:
        energy_tol = 1e20
    print('Calculating fingerprints...')

    fingerprint_list = [None for doc in cursor]
    required_inds = set()

    # scipy sparse matrices dont seem to allow non-zero default values, so we'll use a defaultdict
    sim_mat = defaultdict(lambda: 1e10)
    print('Assessing similarities...')
    for i in range(len(fingerprint_list)):
        for j in range(i+1, len(fingerprint_list)):
            # are we checking stoichiometries, if so, ensure they're the same
            if (enforce_same_stoich is False or
                    (sorted(cursor[j]['stoichiometry']) == sorted(cursor[i]['stoichiometry']) and
                     np.abs(cursor[j].get('enthalpy_per_atom', 0) - cursor[i].get('enthalpy_per_atom', 0)) < energy_tol)):
                # need to set both to zero so we can iterate over the dict later
                sim_mat[i, j] = None
                sim_mat[j, i] = None
                required_inds.add(i)
                required_inds.add(j)

    factory = PDFFactory(cursor, required_inds=list(required_inds), **fingerprint_calc_args)

    for i, j in sim_mat:
        if sim_mat[i, j] is None:
            sim = cursor[i][factory.default_key].get_sim_distance(cursor[j][factory.default_key])
            sim_mat[i, j] = sim
            sim_mat[j, i] = sim

    distinct_set = set()
    dupe_set = set()
    dupe_dict = dict()
    for i in range(len(cursor)):
        distinct_set.add(i)
        dupe_dict[i] = set()

    # loop over the similarity matrix and construct the set of "unique" structures
    # and a dictionary containing their duplicates
    for i, j in sim_mat:
        if sim_mat[i, j] <= sim_tol:
            if i not in dupe_set:
                if j in distinct_set:
                    distinct_set.remove(j)
                    del dupe_dict[j]
                dupe_set.add(j)
                dupe_dict[i].add(j)

    total_dupes = len(set(list(dupe_dict.keys()) + [item for key in dupe_dict for item in dupe_dict[key]]))
    if len(cursor) != total_dupes:
        raise RuntimeError("Something went wrong: dupe dict had wrong size {} compared to cursor {}!\nFull output: {}"
                           .format(total_dupes, len(cursor), dupe_dict))

    if hierarchy_order is None:
        hierarchy_order = ['ICSD', 'DOI', 'OQMD', 'MP', 'PF', 'SWAPS', 'AIRSS', 'GA']
    if hierarchy_values is None:
        hierarchy_values = [get_guess_doc_provenance(doc['source']) for doc in cursor]

    print('Applying hierarchy of structures with order: {}'.format(hierarchy_order))
    dupe_dict = _enforce_hierarchy(dupe_dict, hierarchy_values, hierarchy_order)

    all_structures = set(list(dupe_dict.keys()) + [item for key in dupe_dict for item in dupe_dict[key]])
    if len(cursor) != len(all_structures):
        raise RuntimeError("Something went wrong: dupe dict had wrong size {} compared to cursor {}!\nDifference: {}"
                           .format(len(all_structures),
                                   len(cursor),
                                   all_structures.symmetric_difference({i for i in range(len(cursor))})))

    print('Done!')
    return sorted(list(dupe_dict.keys())), dupe_dict, fingerprint_list, sim_mat


def _enforce_hierarchy(dupe_dict, values, hierarchy):
    """ Enforce a general hierarchy of which structures to keep, based
    on the list of values and their importance.

    Parameters:
        dupe_dict (dict): the dictionary keyed by the index of unique structures
            that holds lists of duplicates for that structure.
        values (list): the list of values for each structure on which to enforce
            the hierarchy.
        hierarchy (list): the order in which to consider the values, e.g.
            `['ICSD', 'OQMD']` will promote ICSD structures over OQMD.

    Returns:
        dict: the reshuffled dictionary of duplicates.

    """
    max_val = max(list(dupe_dict.keys()) + [val for t in dupe_dict.values() for val in t])

    if len(values) - 1 != max_val:
        raise RuntimeError("Number of hierarchy values does not much number of items: {} vs {}"
                           .format(len(values)-1, max_val))

    new_dupe_dict = copy.deepcopy(dupe_dict)

    swapped = []
    for i in new_dupe_dict:
        if not list(new_dupe_dict[i]):
            continue
        for value in hierarchy:
            found = False
            for k in [i] + list(new_dupe_dict[i]):
                if values[k] == value:
                    swapped.append((i, k))
                    found = True
                    break
            if found:
                break

    for i, k in swapped:
        if i != k:
            if k in new_dupe_dict:
                new_dupe_dict[k].update([ind for ind in new_dupe_dict[i] if ind != k] + [i])
            else:
                new_dupe_dict[k] = set([ind for ind in new_dupe_dict[i] if ind != k] + [i])
            del new_dupe_dict[i]

    return new_dupe_dict
