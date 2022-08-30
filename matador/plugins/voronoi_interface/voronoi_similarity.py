# encoding: utf-8
""" This file implements lattice-site level measures of similarity, using
Voronoi decompositions.

Sketch of methodology:

* For a list of structures, take each in turn, compute the Voronoi substructure
and create normalised padded arrays so that they can be compared.

* For each structure, find the *unique* sites, to some definition of unique, initially
just a simple np.isclose() on the site arrays.

* Now do an all-to-all comparison of the unique sites in each structure, yielding
an overall list of unique substructures. This step should have a dial that can be turned
such that all sites fold onto each other, or all sites become distinct, i.e. sensitivity
vs specificity.

* [OPTIONAL] The remaining unique sites across the whole list can now be clustered, if desired,
using e.g. hierarchical clustering.

* Every site in every structure can now be assigned to one of the unique sites above, or
one of the unique clusters.

* Finally, the structures themselves can be clustered by the sites that are present, if desired.

"""
from collections import defaultdict
import numpy as np

from matador.utils.cell_utils import frac2cart


def are_sites_the_same(site_A, site_B, rtol=1e-2, atol=1e-2):
    """Simple check for uniqueness of sites based on their
    Voronoi substructure.

    Input:

        | site_A : dict, with elements as keys, containing numpy arrays
                   of normalised solid angles.
        | site_B : dict, as A.

    Args:

        | rtol : relative tolerance in solid angle for np.allclose,
        | atol : absolute tolerance in solid angle for np.allclose.

    Returns:

        | True if sites have approximately the same values for each species,
          else False.

    """
    # if a key is missing from A or B, check that its not some extraneous distant atom
    for key in site_B:
        if key not in site_A:
            site_A_missing_key = np.zeros_like(site_B[key])
            if not (
                np.allclose(site_A_missing_key, site_B[key], atol=atol, rtol=rtol)
                or np.allclose(site_B[key], site_A_missing_key, atol=atol, rtol=rtol)
            ):
                return False
    # now check remaining keys
    for key in site_A:
        if key not in site_B:
            site_B_missing_key = np.zeros_like(site_A[key])
            if not (
                np.allclose(site_B_missing_key, site_A[key], atol=atol, rtol=rtol)
                or np.allclose(site_A[key], site_B_missing_key, atol=atol, rtol=rtol)
            ):
                return False
        else:
            if len(site_A[key]) != len(site_B[key]):
                if len(site_A[key]) > len(site_B[key]):
                    padded_B = np.append(
                        site_B[key], np.zeros((len(site_A[key]) - len(site_B[key])))
                    )
                    if not (
                        np.allclose(site_A[key], padded_B, atol=atol, rtol=rtol)
                        or np.allclose(site_A[key], padded_B, atol=atol, rtol=rtol)
                    ):
                        return False
                else:
                    padded_A = np.append(
                        site_A[key], np.zeros((len(site_B[key]) - len(site_A[key])))
                    )
                    assert np.allclose(
                        padded_A, site_B[key], atol=atol, rtol=rtol
                    ) is np.allclose(site_B[key], padded_A, atol=atol, rtol=rtol)
                    if not (
                        np.allclose(padded_A, site_B[key], atol=atol, rtol=rtol)
                        or np.allclose(site_B[key], padded_A, atol=atol, rtol=rtol)
                    ):
                        return False
            else:
                if not (
                    np.allclose(site_A[key], site_B[key], atol=atol, rtol=rtol)
                    or np.allclose(site_B[key], site_A[key], rtol=rtol, atol=atol)
                ):
                    return False

    return True


def collect_unique_sites(cursor):
    """Collect unique substrucs from each structure into a single dict.

    Input:

        | cursor: list(dict), list of structures with pre-computed unique sites.

    Returns:

        | unique_environments: dict(list), dict with element symbol keys containing
                               a list of each unique substructure.

    """
    unique_environments = defaultdict(list)
    for doc in cursor:
        for elem in doc["unique_substrucs"]:
            for substruc in doc["unique_substrucs"][elem]:
                unique_environments[elem].append(substruc)
    return unique_environments


def get_max_coordination_of_elem(single_elem_environments):
    """For a given set of elemental environments, find the greatest
    number of each element in a site.

    Input:

        | single_elem_environments: list(dict), list of voronoi substructures,

    Returns:

        | max_num_elem: dict(int), greatest number of each element in a substruc.

    """
    max_num_elem = {}
    for site_ind, site in enumerate(single_elem_environments):
        for elem, value in site:
            elem_len = len([val for val in site if val[0] == elem])
            if elem not in max_num_elem:
                max_num_elem[elem] = 0
            if elem_len > max_num_elem[elem]:
                max_num_elem[elem] = elem_len

    return max_num_elem


def create_site_array(
    unique_environments, max_num_elems=None, elems=None, normalise=False
):
    """Create padded numpy arrays based on unique environments provided.

    Input:

        | unique_environments: dict(list), dict with element keys full of
                               local substructures.

    Args:

        | max_num_elems : dict(dict(int)), dict with element keys, containing
                          sub-dict with element keys that contain the largest
                          number of each type of element contributing to a site.
                          If None, this will be calculated based on the input
                          environments. The site array is padded to this size.
        | elems         : list(str), custom list of elements to loop over.
        | normalise     : bool, whether to normalise site array to 1.

    Returns:

        | site_array    : dict(np.ndarray), dict with element keys full of
                          normalised site arrays.
        | max_num_elems : dict(dict(int)), dict with element keys containing the
                          highest number of atoms of an element contributing
                          to a particular site, e.g. {'P': {'K': 10}}.

    """
    site_array = defaultdict(list)
    if elems is None:
        elems = [elem for elem in unique_environments]
    if max_num_elems is None:
        max_num_elems = dict()
        for elem in elems:
            max_num_elems[elem] = get_max_coordination_of_elem(
                unique_environments[elem]
            )
    for elem in elems:
        for site_ind, site in enumerate(unique_environments[elem]):
            site_array[elem].append(dict())
            total_angle = 0
            for _elem in elems:
                if _elem not in max_num_elems[elem]:
                    max_num_elems[elem][_elem] = 0
                site_array[elem][-1][_elem] = np.zeros((max_num_elems[elem][_elem]))
                count = 0
                for species, angle in site:
                    if species == _elem:
                        site_array[elem][-1][_elem][count] = angle
                        count += 1
                        total_angle += angle
            if normalise:
                for _elem in elems:
                    site_array[elem][-1][_elem] /= total_angle

    return site_array, max_num_elems


def set_substruc_dict(doc):
    """Compute voronoi substructure and collect into dict. Sets the
    'voronoi_substruc' and 'substruc_dict' entries in the doc.

    Input:

        | doc: dict, structure to compute substructure of.

    """
    if "voronoi_substruc" not in doc:
        from matador.plugins.voronoi_interface.voronoi_interface import (
            get_voronoi_substructure,
        )

        doc["voronoi_substruc"] = get_voronoi_substructure(doc)
    voronoi_substruc_dict = dict()
    elems = set(doc["atom_types"])
    for elem in elems:
        voronoi_substruc_dict[elem] = [
            substruc[1] for substruc in doc["voronoi_substruc"] if substruc[0] == elem
        ]
    doc["substruc_dict"] = voronoi_substruc_dict


def set_site_array(doc, normalise=False):
    """Set the 'site_array' entry in the chosen document. Creates
    a dict of numpy arrays of normalised solid angles with element keys
    from the import calculated substructure.

    Input:

        | doc: dict, matador doc containing substruc_dict.

    Args:

        | normalise: bool, whether to normalise site_arrays to total angle.

    """
    if "substruc_dict" not in doc:
        set_substruc_dict(doc)
    doc["site_array"], doc["max_num_elems"] = create_site_array(
        doc["substruc_dict"], normalise=normalise
    )


def get_unique_sites(doc, atol=1e-2, rtol=1e-2):
    if "substruc_dict" not in doc:
        set_substruc_dict(doc)
    elems = set(doc["atom_types"])
    substruc_dict = doc["substruc_dict"]
    if "site_array" not in doc:
        set_site_array(doc)
    site_array = doc["site_array"]
    unique_sites = dict()
    degeneracies = dict()
    similar_sites = dict()
    for elem in elems:
        unique_sites[elem] = set()
        similar_sites[elem] = []
        for i in range(len(site_array[elem])):
            for j in range(i + 1, len(site_array[elem])):
                same = are_sites_the_same(
                    site_array[elem][i], site_array[elem][j], atol=atol, rtol=rtol
                )
                if same:
                    added = False
                    for ind, site in enumerate(similar_sites[elem]):
                        if i in site:
                            similar_sites[elem][ind].add(j)
                            added = True
                            break
                        if j in site:
                            similar_sites[elem][ind].add(i)
                            added = True
                            break
                    if not added:
                        similar_sites[elem].append(set([i, j]))
            if not any([i in _set for _set in similar_sites[elem]]):
                similar_sites[elem].append(set([i]))

        # for env in site_array[elem]:
        # print(env)
        from copy import deepcopy

        temp_similar_sites = deepcopy(similar_sites[elem])
        valid = False
        _iter = 0
        while not valid:
            scrubbed = False
            _iter += 1
            if _iter > 100:
                raise RuntimeError
            for index in range(len(site_array[elem])):
                index_sites = []
                for ind, site in enumerate(temp_similar_sites):
                    if index in site:
                        index_sites.append(ind)
                # if index is in multiple sites, combine all sites into the first
                if len(index_sites) > 1:
                    # construct mean of each set of sites, and if they are too far apart, get worried
                    site_mean = []
                    for j in range(len(index_sites)):
                        site_mean.append(defaultdict(list))
                        for _elem in site_array[elem][i]:
                            site_mean[-1][_elem].append(
                                np.mean(
                                    np.asarray(
                                        [
                                            site_array[elem][i][_elem]
                                            for i in temp_similar_sites[index_sites[j]]
                                        ]
                                    ),
                                    axis=0,
                                )
                            )
                            assert np.shape(site_mean[-1][_elem][-1]) == np.shape(
                                site_array[elem][0][_elem]
                            )
                    for j in range(len(index_sites) - 1, 0, -1):
                        # let means deviate by twice the amount
                        if not are_sites_the_same(
                            site_mean[j], site_mean[0], atol=5 * atol, rtol=5 * rtol
                        ):
                            raise RuntimeError(
                                "Numerical instability: site groups with distinct means were clustered."
                            )
                        [
                            temp_similar_sites[index_sites[0]].add(i)
                            for i in list(temp_similar_sites[index_sites[j]])
                        ]
                        del temp_similar_sites[index_sites[j]]
                    scrubbed = True
                    break
                elif len(index_sites) == 0:
                    raise RuntimeError
            if not scrubbed:
                valid = True
        similar_sites[elem] = temp_similar_sites

        unique_sites[elem] = [next(iter(site)) for site in similar_sites[elem]]
        degeneracies[elem] = [len(_set) for _set in similar_sites[elem]]

    doc["unique_site_inds"] = unique_sites
    doc["unique_substrucs"] = dict()
    doc["unique_site_array"] = dict()
    doc["unique_site_stddev"] = dict()
    doc["similar_sites"] = similar_sites
    for elem in elems:
        doc["unique_substrucs"][elem] = [
            substruc_dict[elem][ind] for ind in unique_sites[elem]
        ]
        doc["unique_site_array"][elem] = []
        doc["unique_site_stddev"][elem] = []
        for site in similar_sites[elem]:
            doc["unique_site_array"][elem].append(dict())
            doc["unique_site_stddev"][elem].append(dict())
            for _elem in site_array[elem][0]:
                doc["unique_site_array"][elem][-1][_elem] = np.mean(
                    np.asarray([site_array[elem][i][_elem] for i in site]), axis=0
                )
                doc["unique_site_stddev"][elem][-1][_elem] = np.std(
                    np.asarray([site_array[elem][i][_elem] for i in site]), axis=0
                )
    doc["site_degeneracies"] = degeneracies


def compare_docs(docA, docB, elems):
    matching = 0
    print("COMPARING", " ".join(docA["text_id"]), "vs", " ".join(docB["text_id"]))
    for elem in elems:
        if elem not in docA["unique_substrucs"] or elem not in docB["unique_substrucs"]:
            print("Missing elements.")
            return False
        for i, site in enumerate(docA["unique_substrucs"][elem]):
            for other_site in docB["unique_substrucs"][elem]:
                if are_sites_the_same(site, other_site, ["K", "P"]):
                    matching += docA["site_degeneracies"][elem][i]
    print(matching, "/", sum(docA["site_degeneracies"][elem]))
    if matching / sum(docA["site_degeneracies"][elem]) > 0.5:
        return True


def cluster(
    X,
    hull,
    max_of_elems,
    elems=None,
    method="KMeans",
    elem="P",
    quantile=None,
    n_clusters=None,
    cmap="tab10",
):
    if method == "KMeans":
        from sklearn.cluster import KMeans

        inertia = []
        elem = "P"
        est = KMeans(n_clusters=n_clusters, n_jobs=-2, precompute_distances=True)
        est.fit(X)
        inertia.append(est.inertia_)
        labels = est.labels_
        # means = est.cluster_centers_
        print("Number of clusters: {}".format(n_clusters))
    elif method == "MeanShift":
        from sklearn.cluster import MeanShift, estimate_bandwidth

        bandwidth = estimate_bandwidth(
            X, quantile=quantile if quantile is not None else 0.2
        )
        est = MeanShift(bandwidth=bandwidth, cluster_all=False)
        est.fit(X)
        labels = est.labels_
        # means = est.cluster_centers_
        n_clusters = len(set(labels))
        print("Number of estimate clusters:", n_clusters)
    """
    means_of_point_class = np.asarray([means[label] for label in labels])
    labels = labels[np.argsort(means_of_point_class[:, 0])]
    X = X[np.argsort(means_of_point_class[:, 0])]
    temp_labels = np.zeros_like(labels)
    sorted_labels = np.zeros_like(labels)
    for i in range(len(labels)-1):
       if labels[i] != labels[i+1]:
           temp_labels[i+1] += 1
    current_group = 0
    for i in range(len(temp_labels)):
       if temp_labels[i] == 1:
           current_group += 1
       sorted_labels[i] = current_group
    """

    #    fig, axarr = plt.subplots(2, 2, figsize=(10,10))
    #    sns.set_palette(sns.color_palette(cmap, n_colors=n_clusters))
    #    colours = [sns.palettes.color_palette(cmap, n_colors=n_clusters)[label] for label in labels]
    #
    #    sns.stripplot(x=labels, y=plot_X[:, 0], ax=axarr[0, 0], jitter=0.1)
    #    axarr[0, 0].set_ylabel('Fraction of solid angle seen as P')
    #
    #    axarr[0, 1].scatter(plot_X[:, 1], plot_X[:, 3], c=colours, lw=0, alpha=0.3)
    #    axarr[0, 1].set_xlabel('Fraction of solid angle seen as P')
    #    axarr[0, 1].set_ylabel('Number of contrib P atoms')
    #
    #    axarr[1, 0].scatter(plot_X[:, 0], plot_X[:, 2], c=colours, lw=0, alpha=0.3)
    #    axarr[1, 0].set_xlabel('Fraction of solid angle seen as K')
    #    axarr[1, 0].set_ylabel('Number of contrib K atoms')
    #
    #    sns.stripplot(x=labels, y=plot_X[:, 3], ax=axarr[1, 1], jitter=0.1)
    #    axarr[1, 1].set_ylabel('Number of contrib Pf atoms')
    #    fig.suptitle('Number of clusters: {}'.format(n_clusters))

    for doc in hull.hull_cursor:
        site_classes = []
        if elem not in doc["unique_substrucs"]:
            doc["class"] = -1
            doc["class_err"] = 0
            doc["class_sites"] = [-1]
            continue
        for ind, substruc in enumerate(doc["unique_substrucs"][elem]):
            site_array = np.zeros((1, int(np.sum(max_of_elems[elem]))))
            elem_counters = np.zeros((len(elems)), dtype=int)
            for elem_ind, _elem in enumerate(elems):
                if elem_ind != 0:
                    elem_counters[elem_ind] += max_of_elems[elem][elem_ind - 1]
            for elem_ind, _elem in enumerate(elems):
                for atom in substruc:
                    if atom[0] == _elem:
                        site_array[0][elem_counters[elem_ind]] = atom[1]
                        elem_counters[elem_ind] += 1
            site_classes.append(est.predict(site_array))
        site_classes = [val[0] for val in site_classes]
        doc["class"] = np.median(site_classes)
        doc["class_sites"] = [val for val in site_classes]
        doc["class_err"] = (len(set(doc["class_sites"])) - 1) / n_clusters
        class_vector = np.zeros((n_clusters))
        for site in site_classes:
            class_vector[int(site)] += 1
        doc["class_vector"] = class_vector

    return est


def plot_class_hull(hull, n_clusters, plot_class=None):
    from random import sample
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12, 20))
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax2 = hull.plot_2d_hull(show=False, ax=ax2)
    colours = plt.cm.Dark2(np.linspace(0, 1, n_clusters * 10))
    class_concs = defaultdict(list)
    for doc in hull.hull_cursor:
        class_concs[int(doc["class"])].append(doc["concentration"])
    #        ax.errorbar(doc['concentration'], doc['class'], yerr=doc['class_err'], fmt='o', c=colours[int(10*doc['class'])], alpha=0.01, zorder=100)
    for _class in class_concs:
        ax.errorbar(
            np.mean(class_concs[_class]),
            _class,
            xerr=np.std(class_concs[_class]),
            fmt="o",
            c=colours[int(10 * _class)],
            lw=2,
            zorder=1000,
        )
    bounds = [0.32, 0.57, 0.74]
    for bound in bounds:
        ax.axvline(bound, ls="--", lw=1)
    ax.set_ylim(-1.5, n_clusters)
    ax.set_xlim(-0.1, 1.1)
    for doc in hull.hull_cursor:
        ax2.scatter(
            doc["concentration"],
            doc["formation_enthalpy_per_atom"],
            c="#eeeeee",
            s=60,
            zorder=9999,
            lw=0,
        )
    for doc in sample(hull.hull_cursor, len(hull.hull_cursor)):
        if plot_class is not None:
            if doc["class"] not in plot_class:
                continue
            else:
                ax2.scatter(
                    doc["concentration"],
                    doc["formation_enthalpy_per_atom"],
                    c=colours[int(10 * (doc["class"]))],
                    s=50,
                    alpha=max(1 - 2 * doc["class_err"], 0.1),
                    zorder=99999,
                    lw=0,
                )
        else:
            if doc["hull_distance"] <= 1e-12:
                ax2.scatter(
                    doc["concentration"],
                    doc["formation_enthalpy_per_atom"],
                    c=colours[int(10 * (doc["class"]))],
                    s=75,
                    zorder=99999999,
                    lw=1,
                )
            else:
                ax2.scatter(
                    doc["concentration"],
                    doc["formation_enthalpy_per_atom"],
                    c=colours[int(10 * (doc["class"]))],
                    s=50,
                    alpha=max(1 - 2 * doc["class_err"], 0.1),
                    zorder=99999,
                    lw=0,
                )


def find_site_index(doc, elem, elem_ind):
    elem_count = 0
    for idx, atom in enumerate(doc["atom_types"]):
        if atom == elem:
            if elem_count == elem_ind:
                print("Found substruc at {}".format(elem_count))
                return elem_count
        elem_count += 1
    return False


def viz_site(doc, targ_substruc, elem, rmax=6):
    if "positions_abs" not in doc:
        doc["positions_abs"] = frac2cart(doc["lattice_cart"], doc["positions_frac"])
    for ind, substruc in enumerate(doc["substruc_dict"][elem]):
        if substruc == targ_substruc:
            elem_ind = ind
            elem_count = find_site_index(doc, elem, elem_ind)
            if elem_count is not False:
                break
    targ_site = doc["positions_frac"][elem_count]
    targ_pos = doc["positions_abs"][elem_count]
    print(elem_count)
    from itertools import product
    from collections import defaultdict

    doc["lattice_cart"] = np.asarray(doc["lattice_cart"])
    neighbour_doc = defaultdict(list)
    neighbour_doc["positions_abs"].append(targ_pos)
    neighbour_doc["positions_frac"].append(targ_site)
    neighbour_doc["atom_types"].append("U")
    neighbour_doc["lattice_cart"] = doc["lattice_cart"]
    for j, pos in enumerate(doc["positions_abs"]):
        for prod in product(range(-1, 2), repeat=3):
            trans = np.zeros((3))
            for ind, multi in enumerate(prod):
                trans += doc["lattice_cart"][ind] * multi
            dist2 = 0
            for i in range(3):
                dist2 += (targ_pos[i] - (pos + trans)[i]) ** 2
            if dist2 < rmax**2:
                neighbour_doc["positions_abs"].append(pos + trans)
                neighbour_doc["positions_frac"].append(
                    (np.asarray(doc["positions_frac"][j]) + prod).tolist()
                )
                neighbour_doc["atom_types"].append(doc["atom_types"][j])
    print(targ_site)
    return elem_count, neighbour_doc
