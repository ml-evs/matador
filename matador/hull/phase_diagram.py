# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements the base PhaseDiagram creator that interfaces
with QueryConvexHull and EnsembleHull.

"""


from traceback import print_exc
import bisect

import scipy.spatial
import numpy as np

from matador.utils.hull_utils import (
    barycentric2cart, vertices2plane, vertices2line, FakeHull, is_point_in_triangle
)
from matador.utils.chem_utils import get_formula_from_stoich
from matador.utils.cursor_utils import get_array_from_cursor, display_results, set_cursor_from_array

EPS = 1e-12


class PhaseDiagram:
    """ This class encapsulates the actual phase data, e.g. the actual
    energy and compositions found to be stable.

    Attributes:
        structures (numpy.ndarray): the array passed to init used to
            make the hull, with the first (num_species-1) columns
            containing normalised concentrations, and the final column
            containing formation energy.
        convex_hull (scipy.spatial.ConvexHull): the actual convex hull
            returned by SciPy.
        formation_key (list): index/key specification of formation energy
            per atom from top level of each document.

    """
    def __init__(self, cursor, formation_key, dimension):
        """ Compute the convex hull of data passed, to retrieve hull
        distances and thus stable structures.

        Parameters:
            cursor (list[dict]): list of matador documents to make
                phase diagram from.
            formation_key (str or list):  location of the formation energy
                inside each document, either a single key or iterable of
                keys to use with `recursive_get`.

        """
        self._dimension = dimension
        self.cursor = cursor
        self.formation_key = formation_key

        structures = np.hstack((
            get_array_from_cursor(cursor, 'concentration').reshape(len(cursor), dimension-1),
            get_array_from_cursor(cursor, self.formation_key).reshape(len(cursor), 1)))

        # define self._structure_slice as the filtered array of points actually used to create the convex hull
        # which can include/exclude points from the passed structures. This array is the one indexed by
        # vertices/simplices in ConvexHull

        if self._dimension == 3:
            # add a point "above" the hull
            # for simple removal of extraneous vertices (e.g. top of 2D hull)
            dummy_point = [0.333, 0.333, 1e5]
            # if ternary, use all structures, not just those with negative eform for compatibility reasons
            self._structure_slice = np.vstack((structures, dummy_point))
        else:
            # filter out those with positive formation energy, to reduce expense computing hull
            self._structure_slice = structures[np.where(structures[:, -1] <= 0 + EPS)]

        # filter out "duplicates" in _structure_slice
        # this prevents breakages if no structures are on the hull and chempots are duplicated
        # but it might be faster to hardcode this case individually
        self._structure_slice = np.unique(self._structure_slice, axis=0)

        # if we only have the chempots (or worse) with negative formation energy, don't even make the hull
        if len(self._structure_slice) <= dimension:
            if len(self._structure_slice) < dimension:
                raise RuntimeError('No chemical potentials on hull... either mysterious use of custom chempots, or worry!')
            self.convex_hull = FakeHull()
        else:
            try:
                self.convex_hull = scipy.spatial.ConvexHull(self._structure_slice)
            except scipy.spatial.qhull.QhullError:
                print(self._structure_slice)
                print('Error with QHull, plotting formation energies only...')
                print_exc()
                self.convex_hull = FakeHull()

        # remove vertices that have positive formation energy
        filtered_vertices = [vertex for vertex in self.convex_hull.vertices if self._structure_slice[vertex, -1] <= 0 + EPS]
        bad_simplices = set()
        for ind, simplex in enumerate(self.convex_hull.simplices):
            for vertex in simplex:
                if vertex not in filtered_vertices:
                    bad_simplices.add(ind)

        filtered_simplices = [simplex for ind, simplex in enumerate(self.convex_hull.simplices) if ind not in bad_simplices]

        self.convex_hull = FakeHull()
        self.convex_hull.points = self._structure_slice
        self.convex_hull.vertices = list(filtered_vertices)
        self.convex_hull.simplices = list(filtered_simplices)

        self.hull_dist = self.get_hull_distances(structures, precompute=True)
        set_cursor_from_array(self.cursor, self.hull_dist, 'hull_distance')
        self.structures = structures
        self.stable_structures = [doc for doc in self.cursor if doc['hull_distance'] < EPS]

    def __str__(self):
        """ Print underlying phase diagram. """
        return display_results(self.cursor,
                               hull=True,
                               colour=False,
                               energy_key=self.formation_key,
                               sort=False,
                               return_str=True)

    def get_hull_distances(self, structures, precompute=False, **kwargs):
        """ Returns array of distances to pre-computed binary or ternary
        hull, from array containing concentrations and energies.

        Parameters:
            structures (numpy.ndarray): N x n array of concentrations and
                enthalpies for N structures, with up to 2 columns of
                concentrations and the last column containing the
                structure's formation enthalpy.

        Keyword arguments:
            precompute (bool): whether or not to bootstrap hull
                distances from previously computed values at the same
                stoichiometry.

        Returns:
            numpy.ndarray: N-dim array storing distances to
                the hull for N structures,

        """

        if precompute:
            # dict with formula keys, containing tuple of pre-computed enthalpy/atom and hull distance
            cached_formula_dists = dict()
            cache_hits = 0
            cache_misses = 0

        if isinstance(structures, list):
            structures = np.asarray(structures)

        # if only chem pots on hull, dist = energy
        if len(self._structure_slice) == self._dimension:
            hull_dist = np.ones((len(structures)))
            hull_dist = structures[:, -1]

        # if binary hull, do binary search
        elif self._dimension == 2:
            tie_line_comp = self._structure_slice[self.convex_hull.vertices, 0]
            tie_line_energy = self._structure_slice[self.convex_hull.vertices, -1]
            tie_line_comp = np.asarray(tie_line_comp)
            tie_line_energy = tie_line_energy[np.argsort(tie_line_comp)]
            tie_line_comp = tie_line_comp[np.argsort(tie_line_comp)]

            hull_dist = np.empty((len(structures)))
            hull_dist.fill(np.nan)
            if precompute:
                for ind, _ in enumerate(structures):
                    formula = get_formula_from_stoich(self.cursor[ind]['stoichiometry'], sort=True, tex=False)
                    if formula in cached_formula_dists:
                        hull_dist[ind] = (structures[ind, -1] - cached_formula_dists[formula][0] +
                                          cached_formula_dists[formula][1])
                        cache_hits += 1
                    else:
                        i = bisect.bisect_left(tie_line_comp, structures[ind, 0])
                        gradient, intercept = vertices2line([[tie_line_comp[i-1], tie_line_energy[i-1]],
                                                             [tie_line_comp[i], tie_line_energy[i]]])
                        # calculate hull_dist
                        hull_dist[ind] = structures[ind, -1] - (gradient * structures[ind, 0] + intercept)
                        cached_formula_dists[formula] = (structures[ind, -1], hull_dist[ind])
                        cache_misses += 1
            else:
                for ind, _ in enumerate(structures):
                    i = bisect.bisect_left(tie_line_comp, structures[ind, 0])
                    gradient, intercept = vertices2line([[tie_line_comp[i-1], tie_line_energy[i-1]],
                                                         [tie_line_comp[i], tie_line_energy[i]]])
                    # calculate hull_dist
                    hull_dist[ind] = structures[ind, -1] - (gradient * structures[ind, 0] + intercept)

        # if ternary, use barycentric coords
        elif self._dimension == 3:
            # loop through structures and find which plane they correspond to
            # using barycentric coordinates, if a formula has already been
            # computed then calculate delta relative to that and skip
            self.convex_hull.planes = [[self._structure_slice[vertex] for vertex in simplex]
                                       for simplex in self.convex_hull.simplices]
            structures_finished = [False] * len(structures)
            hull_dist = np.empty(len(structures))
            hull_dist.fill(np.nan)
            cart_planes_inv = []
            planes_height_fn = []
            for ind, plane in enumerate(self.convex_hull.planes):
                cart_planes = barycentric2cart(plane).T
                cart_planes[-1, :] = 1
                # if projection of triangle in 2D is a line, do binary search
                if np.linalg.det(cart_planes) == 0:
                    cart_planes_inv.append(None)
                    planes_height_fn.append(None)
                else:
                    cart_planes_inv.append(np.linalg.inv(cart_planes))
                    planes_height_fn.append(vertices2plane(plane))
            for idx, structure in enumerate(structures):
                for ind, plane in enumerate(self.convex_hull.planes):
                    if cart_planes_inv[ind] is None:
                        continue
                    if precompute and get_formula_from_stoich(self.cursor[idx]['stoichiometry'], sort=True,
                                                              tex=False) in cached_formula_dists:
                        formula = get_formula_from_stoich(self.cursor[idx]['stoichiometry'], sort=True, tex=False)
                        if formula in cached_formula_dists:
                            cache_hits += 1
                            hull_dist[idx] = (structures[idx, -1] - cached_formula_dists[formula][0] +
                                              cached_formula_dists[formula][1])
                            structures_finished[idx] = True

                    elif is_point_in_triangle(structure, cart_planes_inv[ind], preprocessed_triangle=True):
                        structures_finished[idx] = True
                        hull_dist[idx] = planes_height_fn[ind](structure)
                        if precompute:
                            cached_formula_dists[
                                get_formula_from_stoich(self.cursor[idx]['stoichiometry'], sort=True,
                                                        tex=False)] = (structure[-1], hull_dist[idx])
                            cache_misses += 1
                        break

            # mask values very close to 0 with 0
            hull_dist[np.where(np.abs(hull_dist) < EPS)] = 0

            failed_structures = []
            for ind, structure in enumerate(structures_finished):
                if not structure:
                    failed_structures.append(ind)

            if failed_structures:
                raise RuntimeError('There were issues calculating the hull distance for {} structures.'
                                   .format(len(failed_structures)))

        # otherwise, set to zero until proper N-d distance can be implemented
        else:
            raise NotImplementedError
            # self.hull.planes = [[self._structure_slice[vertex] for vertex in simplex]
            #                     for simplex in self.hull.simplices]
            # for idx, structure in enumerate(structures):
            #     if precompute and get_formula_from_stoich(self.cursor[idx]['stoichiometry'],
            #                                               tex=False) in cached_formula_dists:
            #         formula = get_formula_from_stoich(self.cursor[idx]['stoichiometry'], tex=False)
            #         if formula in cached_formula_dists:
            #             cache_hits += 1
            #             hull_dist[idx] = (structures[idx, -1] - cached_formula_dists[formula][0] +
            #                               cached_formula_dists[formula][1])
            #         structures_finished[idx] = True
            #     else:
            #         if point_in_polygon(structure[:-1], plane):
            #             hull_dist[idx] = 0
            #             continue

            # loop over compositions
            #   loop over planes
            #       check if point is interior to plane
            #           if so, compute "height" above plane, using equation of hyperplane from scipy
        if np.isnan(hull_dist).any():
            raise RuntimeError("Some hull distances failed")

        return hull_dist

    def point_in_polygon(self):
        """ Use "ray-tracing" or winding number approach
        to compute whether a point lies inside a polygon.
        """
        raise NotImplementedError
