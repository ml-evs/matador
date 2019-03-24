# coding: utf-8
""" This file implements the general wrapper to
structural similarity metrics.
"""
# matador modules
from matador.similarity.fingerprint import FingerprintFactory
from matador.similarity.pdf_similarity import PDF
from matador.utils.cursor_utils import get_array_from_cursor, get_guess_doc_provenance
# external libraries
import numpy as np
# standard library
from collections import defaultdict


def get_uniq_cursor(cursor, fingerprint=PDF, sim_tol=0.1, energy_tol=1e-2,
                    enforce_same_stoich=True,
                    debug=False, **fingerprint_calc_args):
    """ Uses fingerprint to filter cursor into unique structures to some
    tolerance sim_tol, additionally returning a dict of duplicates and the
    correlation matrix.

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
        distinct_set (set): a set of indices of unique documents
        dupe_dict (dict): a dict with keys from distinct_set, listing duplicates
        fingerprint_list (list): a list of <Fingerprint> objects
        sim_mat (np.ndarray): the correlation matrix of pair similarity distances

    """
    fingerprint_list = []
    if not enforce_same_stoich:
        energy_tol = 1e20
    if projected:
        fingerprint_calc_args['projected'] = True
    print('Calculating fingerprints...')
    if debug:
        import time
        start = time.time()
    factory = FingerprintFactory(cursor, **fingerprint_calc_args)
    if debug:
        completed = time.time() - start
        print('{} of {} structures completed in {:0.1f} s'.format(fingerprint, len(cursor), completed))

    fingerprint_list = [doc[fingerprint.key] for doc in cursor]
    sim_mat = np.ones((len(fingerprint_list), len(fingerprint_list)))
    print('Assessing similarities...')
    for i in range(len(fingerprint_list)):
        sim_mat[i, i] = 0
        for j in range(i+1, len(fingerprint_list)):
            # are we checking stoichiometries, if so, ensure they're the same
            if (enforce_same_stoich is False or
                (sorted(cursor[j]['stoichiometry']) == sorted(cursor[i]['stoichiometry']) and
                 np.abs(cursor[j].get('enthalpy_per_atom', 0) - cursor[i].get('enthalpy_per_atom', 0)) < energy_tol)):
                sim = fingerprint_list[i].get_sim_distance(fingerprint_list[j], projected=projected)
                sim_mat[i, j] = sim
                sim_mat[j, i] = sim
            else:
                sim = 1e10
                sim_mat[i, j] = sim
                sim_mat[i, j] = sim
    rows, cols = np.where(sim_mat <= sim_tol)
    distinct_set = set()
    dupe_set = set()
    dupe_dict = defaultdict(list)
    prov = [get_guess_doc_provenance(doc['source']) for doc in cursor]
    for i in range(len(sim_mat)):
        distinct_set.add(i)
        dupe_dict[i] = []
    for i, j in zip(rows, cols):
        if i == j:
            continue
        elif i not in dupe_set:
            if j in distinct_set:
                distinct_set.remove(j)
                del dupe_dict[j]
            dupe_set.add(j)
            dupe_dict[i].append(j)

    assert len(cursor) == len(set([key for key in dupe_dict] + [item for key in dupe_dict for item in dupe_dict[key]]))

    for i in dupe_dict:
        to_compare = [i]
        to_compare.extend(dupe_dict[i])
        hierarchy = ['ICSD', 'OQMD', 'SWAPS', 'AIRSS', 'GA']
        swapped = []
        for provenance in hierarchy:
            found = False
            for k in to_compare:
                if prov[k] is provenance:
                    distinct_set.remove(i)
                    distinct_set.add(k)
                    swapped.append((i, k))
                    found = True
                    break
            if found:
                break

    for pair in swapped:
        i, k = pair
        dupe_dict[k] = [ind for ind in dupe_dict[i] if ind != k] + [i]
        del dupe_dict[i]

    print('Done!')
    return distinct_set, dupe_dict, fingerprint_list, sim_mat


def plot_similarity_energy_correlation_matrix(cursor, sim_mat, sim_vmin=0.05, sim_vmax=0.5):
    """ Plot a correlation heatmap where the upper triangular displays
    relative energy differences and the lower triangular displays
    structural similarity distance.

    TO-DO:
    * guidelines for different stoichiometries

    Parameters:
        cursor (list): a matador cursor
        sim_mat (np.array): matrix where S_{ij} = similarity distance(i, j)

    Keyword arguments:
        sim_vmin (float): sim distance at which to show minimum colour
        sim_vmax (float): "---------------------------" maximum "----"

    """
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
    except ImportError as exc:
        raise exc('This is the one function that needs seaborn throughout matador... please install it.')
    assert(len(cursor) == len(sim_mat))
    energy_list = get_array_from_cursor(cursor, 'enthalpy_per_atom')
    energy_mat = np.zeros_like(sim_mat)
    for i in range(len(sim_mat)):
        for j in range(i+1, len(sim_mat)):
            energy_mat[i][j] = (energy_list[j] - energy_list[i])
    _sim_mat = np.tril(sim_mat)
    sim_mask = np.zeros_like(sim_mat, dtype=np.bool)
    sim_mask[np.triu_indices(len(sim_mask), k=1)] = True

    _energy_mat = np.triu(energy_mat)
    energy_mask = np.zeros_like(sim_mat, dtype=np.bool)
    energy_mask[np.tril_indices(len(sim_mask), k=1)] = True

    _, axarr = plt.subplots(ncols=2, nrows=2, figsize=(10, 10),
                            gridspec_kw={'height_ratios': [20, 1],
                                         'width_ratios': [20, 1]})
    cmap_energy = sns.cubehelix_palette(as_cmap=True, rot=-.4, light=1, start=1.5)
    cmap_sim = sns.cubehelix_palette(as_cmap=True, rot=.4, light=1, start=1.5)

    ax = axarr[0][0]
    axarr[1][1].axis('off')

    sns.heatmap(_sim_mat, mask=sim_mask, cmap=cmap_sim,
                vmin=sim_vmin, vmax=sim_vmax, ax=ax,
                xticklabels=False, yticklabels=False, square=True,
                cbar_kws={'label': 'structural similarity distance',
                          'orientation': 'horizontal'},
                cbar_ax=axarr[1][0])
    sns.heatmap(_energy_mat, mask=energy_mask, cmap=cmap_energy,
                vmax=0.05, vmin=0.0, ax=ax,
                xticklabels=False, yticklabels=False, square=True,
                cbar_kws={'label': 'enthalpy difference (eV/atom)'},
                cbar_ax=axarr[0][1])
    return
