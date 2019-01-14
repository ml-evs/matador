# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements convex hull functionality from database
queries.

"""


from traceback import print_exc
from bisect import bisect_left
from copy import deepcopy
import sys
import os
import re

from scipy.spatial import ConvexHull
from scipy.spatial.qhull import QhullError
from bson.son import SON
import pymongo as pm
import numpy as np

from matador.utils.print_utils import print_failure, print_notify, print_warning
from matador.utils.hull_utils import barycentric2cart, vertices2plane, vertices2line, FakeHull
from matador.utils.chem_utils import parse_element_string, get_padded_composition, get_num_intercalated
from matador.utils.chem_utils import get_generic_grav_capacity, get_formula_from_stoich, get_stoich_from_formula
from matador.utils.chem_utils import get_formation_energy, KELVIN_TO_EV, get_concentration
from matador.utils.cursor_utils import set_cursor_from_array, get_array_from_cursor
from matador.utils.cursor_utils import display_results
from matador.export import generate_hash

EPS = 1e-12


class QueryConvexHull:
    """ Construct a binary or ternary phase diagram from a
    matador.query.DBQuery object, or a list of structures.

    Attributes:
        cursor (list): list of all structures used to create phase diagram.
        hull_cursor (list): list of all documents within hull_cutoff.
        chempot_cursor (list): list of chemical potential documents.
        structures (numpy.ndarray): all structures used to create hull.
        structure_slice (numpy.ndarray): array of concentrations and
            formation energies for each structure with E_F <= 0. The indices
            in hull correspond to entries of this array.
        hull_dist (np.ndarray): array of distances from hull for each structure.
        species (list): list of chemical potential symbols.
        num_elements (int): number of elements present in the chemical potentials.
        elements (list): the elements present in the convex hull.
        voltage_data (dict): if voltage_curve() has been called, then
            this is a dictionary containing Q, x, V and reaction pathways.
        volume_data (dict): if volume_curve() has been called, then
            this is a dictionary containing x and volumes per B (in AxBy).

    """

    def __init__(self, query=None, cursor=None, elements=None, species=None, subcmd='hull', quiet=False,
                 plot_kwargs=None, **kwargs):
        """ Initialise the class from either a DBQuery or a cursor (list
        of matador dicts) and construct the appropriate phase diagram.

        Keyword arguments:
            query (matador.query.DBQuery): object containing structures,
            cursor (list(dict)): alternatively specify list of matador docs.
            species (list(str)): list of elements/chempots to use, used to provide a useful order,
            elements (list(str)): deprecated form of the above.
            subcmd (str): either 'hull' or 'voltage',
            kwargs (dict): mostly CLI arguments, see matador hull --help for full options.
            plot_kwargs (dict): arguments to pass to plot_hull function

        """
        self.args = dict()
        if query is not None:
            self.args.update(query.args)
        self.args.update(kwargs)

        self.devel = True
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd
        if plot_kwargs is None:
            plot_kwargs = {}
        if plot_kwargs.get('show') is None:
            plot_kwargs['show'] = True
        self._query = query
        self.from_cursor = False
        self.plot_params = False

        if self._query is not None:
            self.cursor = list(query.cursor)
            use_source = False
            if self._query.args['subcmd'] not in ['hull', 'voltage', 'hulldiff']:
                raise RuntimeError('Query was not prepared with subcmd=hull, cannot make a hull...')
        else:
            self.cursor = cursor
            self.from_cursor = True
            use_source = True

        # set up attributes for later
        self.chempot_cursor = None
        self.hull_cursor = None
        self.structure_slice = None
        self.hull_dist = None
        self.species = None
        self.voltage_data = {}
        self.volume_data = {}
        self.elements = []
        self.num_elements = 0

        # set some hull options
        if quiet:
            cached_stdout = sys.stdout
            f = open(os.devnull, 'w')
            sys.stdout = f

        if self.args.get('energy_key') is not None:
            # consider if temperature present and looking for thermo_* keys
            if self.args.get('temperature') is not None:
                if self.args.get('energy_key').startswith('thermo_'):
                    prefix = ''
                else:
                    prefix = 'thermo_'
            else:
                prefix = ''

            if self.args.get('energy_key').endswith('_per_atom'):
                self.energy_key = prefix + self.args.get('energy_key')
            else:
                self.energy_key = prefix + self.args.get('energy_key') + '_per_atom'
        else:
            self.energy_key = 'enthalpy_per_atom'

        self.temperature = self.args.get('temperature')
        if self.args.get('hull_cutoff') is not None:
            self.hull_cutoff = float(self.args['hull_cutoff'])
        else:
            self.hull_cutoff = 0.0

        if self.cursor is None:
            raise RuntimeError('Failed to find structures to create hull!')

        if species is None:
            species = elements

        self._non_elemental = False
        if species is None:
            if elements is None:
                if isinstance(self.args.get('composition'), list):
                    species = self.args.get('composition')[0]
                else:
                    species = self.args.get('composition')
                if ':' in species:
                    species = species.split(':')
            else:
                species = elements

        if isinstance(species, str):
            if ':' in species:
                species = species.split(':')
            else:
                species = [spec for spec in re.split(r'([A-Z][a-z]*)', species) if spec]

        self.species = species
        assert isinstance(self.species, list)
        for _species in self.species:
            if len(parse_element_string(_species, stoich=True)) > 1:
                self._non_elemental = True

        self._dimension = len(self.species)
        if self._dimension > 2 and self._query is not None and not self._query.args.get('intersection'):
            raise SystemExit('Please query with -int/--intersection when creating ternary+ hulls.')

        self.construct_phase_diagram()

        if not self.hull_cursor:
            print_warning('No structures on hull with chosen chemical potentials.')
        else:
            print_notify('{} structures found within {} eV of the hull, including chemical potentials.'
                         .format(len(self.hull_cursor), self.hull_cutoff))

        display_results(self.hull_cursor, self.args, hull=True, use_source=use_source)

        if not self.args.get('no_plot'):
            from matador import plotting

        if self.args['subcmd'] == 'voltage':
            self.voltage_curve([doc for doc in self.hull_cursor if doc['hull_distance'] <= 1e-9])
            if not self.args.get('no_plot'):
                plotting.plot_voltage_curve(self)
                self.plot_hull(**plot_kwargs, debug=self.args.get('debug'))

        if self.args.get('volume'):
            self.volume_curve()
            if not self.args.get('no_plot'):
                plotting.plot_volume_curve(self)

        if self.args['subcmd'] == 'hull' and not self.args.get('no_plot'):
            self.plot_hull(**plot_kwargs, debug=self.args.get('debug'))

        if self.args.get('uniq'):
            from matador.similarity.similarity import get_uniq_cursor
            if self.args.get('uniq') is True:
                sim_tol = 0.1
            else:
                sim_tol = self.args.get('uniq')
            print_notify('Filtering for unique structures...')
            unique_set, _, _, _ = get_uniq_cursor(self.hull_cursor, debug=self.args.get('debug'), sim_tol=sim_tol)
            old_cursor_len = len(self.hull_cursor)
            self.hull_cursor = [self.hull_cursor[ind] for ind in unique_set]
            display_results(self.hull_cursor, hull=True, args=self.args)
            print('Filtered {} down to {}'.format(old_cursor_len, len(unique_set)))

        if quiet:
            f.close()
            sys.stdout = cached_stdout

    @property
    def savefig(self):
        """ True if any figure type argument was passed. """
        return any([self.args.get('pdf'), self.args.get('png'), self.args.get('svg')])

    def plot_hull(self, **kwargs):
        """ Hull plot helper function. """
        from matador import plotting
        if self._dimension == 3:
            ax = plotting.plot_ternary_hull(self, **kwargs)
        elif self._dimension == 2:
            ax = plotting.plot_2d_hull(self, **kwargs)
        else:
            print_notify('Unable to plot phase diagram of dimension {}.'.format(self._dimension))
        return ax

    def fake_chempots(self, custom_elem=None):
        """ Spoof documents for command-line chemical potentials.

        Keyword arguments:
            custom_elem (list(str)): list of element symbols to generate chempots for.

        """
        self.chempot_cursor = []

        if custom_elem is None:
            custom_elem = self.species

        if len(custom_elem) != len(self.species):
            raise RuntimeError('Wrong number of compounds/chemical potentials specified: {} vs {}'.format(custom_elem, self.args.get('chempots')))
        for i, _ in enumerate(self.args.get('chempots')):
            self.chempot_cursor.append(dict())
            self.chempot_cursor[i]['stoichiometry'] = get_stoich_from_formula(custom_elem[i])
            self.chempot_cursor[i]['enthalpy_per_atom'] = -1*abs(self.args.get('chempots')[i])
            self.chempot_cursor[i]['enthalpy'] = self.chempot_cursor[i]['enthalpy_per_atom'] * sum(elem[1] for elem in self.chempot_cursor[i]['stoichiometry'])
            self.chempot_cursor[i]['num_fu'] = 1
            self.chempot_cursor[i]['text_id'] = ['command', 'line']
            self.chempot_cursor[i]['_id'] = generate_hash(hash_len=10)
            self.chempot_cursor[i]['source'] = ['command_line']
            self.chempot_cursor[i]['space_group'] = 'xxx'
            self.chempot_cursor[i]['enthalpy_per_b'] = self.chempot_cursor[i]['enthalpy_per_atom']
            self.chempot_cursor[i]['num_a'] = 0
            self.chempot_cursor[i]['cell_volume'] = 1
            self.chempot_cursor[i]['concentration'] = [1 if i == ind else 0 for ind in range(self._dimension-1)]
        self.chempot_cursor[0]['num_a'] = float('inf')
        notify = 'Custom chempots:'
        for chempot in self.chempot_cursor:
            notify += '{:3} = {} eV/fu, '.format(get_formula_from_stoich(chempot['stoichiometry'], sort=False),
                                                 chempot['enthalpy'])

        if self.args.get('debug'):
            for match in self.chempot_cursor:
                print(match)
        print(len(notify) * '─')
        print(notify)
        print(len(notify) * '─')

    def get_hull_distances(self, structures, precompute=False):
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
        if len(self.structure_slice) == self._dimension:
            hull_dist = np.ones((len(structures)))
            hull_dist = structures[:, -1]

        # if binary hull, do binary search
        elif self._dimension == 2:
            tie_line_comp = self.structure_slice[self.convex_hull.vertices, 0]
            tie_line_energy = self.structure_slice[self.convex_hull.vertices, -1]
            tie_line_comp = np.asarray(tie_line_comp)
            tie_line_energy = tie_line_energy[np.argsort(tie_line_comp)]
            tie_line_comp = tie_line_comp[np.argsort(tie_line_comp)]

            hull_dist = np.ones((len(structures)))
            if precompute:
                for ind, _ in enumerate(structures):
                    formula = get_formula_from_stoich(self.cursor[ind]['stoichiometry'], sort=True, tex=False)
                    if formula in cached_formula_dists:
                        hull_dist[ind] = (structures[ind, -1] - cached_formula_dists[formula][0] +
                                          cached_formula_dists[formula][1])
                        cache_hits += 1
                    else:
                        i = bisect_left(tie_line_comp, structures[ind, 0])
                        gradient, intercept = vertices2line([[tie_line_comp[i-1], tie_line_energy[i-1]],
                                                             [tie_line_comp[i], tie_line_energy[i]]])
                        # calculate hull_dist
                        hull_dist[ind] = structures[ind, -1] - (gradient * structures[ind, 0] + intercept)
                        cached_formula_dists[formula] = (structures[ind, -1], hull_dist[ind])
                        cache_misses += 1
            else:
                for ind, _ in enumerate(structures):
                    i = bisect_left(tie_line_comp, structures[ind, 0])
                    gradient, intercept = vertices2line([[tie_line_comp[i-1], tie_line_energy[i-1]],
                                                         [tie_line_comp[i], tie_line_energy[i]]])
                    # calculate hull_dist
                    hull_dist[ind] = structures[ind, -1] - (gradient * structures[ind, 0] + intercept)

        # if ternary, use barycentric coords
        elif self._dimension == 3:
            # for each plane, convert each point into barycentric coordinates
            # for that plane and test for negative values
            self.convex_hull.planes = [[self.structure_slice[vertex] for vertex in simplex]
                                       for simplex in self.convex_hull.simplices]
            structures_finished = [False] * len(structures)
            hull_dist = np.ones((len(structures) + 1))
            planes_R_inv = []
            planes_height_fn = []
            for ind, plane in enumerate(self.convex_hull.planes):
                R = barycentric2cart(plane).T
                R[-1, :] = 1
                # if projection of triangle in 2D is a line, do binary search
                if np.linalg.det(R) == 0:
                    planes_R_inv.append(None)
                    planes_height_fn.append(None)
                else:
                    planes_R_inv.append(np.linalg.inv(R))
                    planes_height_fn.append(vertices2plane(plane))
            for idx, structure in enumerate(structures):
                for ind, plane in enumerate(self.convex_hull.planes):
                    if structures_finished[idx] or planes_R_inv[ind] is None:
                        continue
                    if precompute and get_formula_from_stoich(self.cursor[idx]['stoichiometry'], sort=True,
                                                              tex=False) in cached_formula_dists:
                        formula = get_formula_from_stoich(self.cursor[idx]['stoichiometry'], sort=True, tex=False)
                        if formula in cached_formula_dists:
                            cache_hits += 1
                            hull_dist[idx] = (structures[idx, -1] - cached_formula_dists[formula][0] +
                                              cached_formula_dists[formula][1])
                            structures_finished[idx] = True
                    else:
                        barycentric_structure = barycentric2cart(structure.reshape(1, 3)).T
                        barycentric_structure[-1, :] = 1
                        plane_barycentric_structure = np.matrix(planes_R_inv[ind]) * np.matrix(barycentric_structure)
                        if (plane_barycentric_structure >= 0 - 1e-12).all():
                            structures_finished[idx] = True
                            hull_dist[idx] = planes_height_fn[ind](structure)
                            if precompute:
                                cached_formula_dists[
                                    get_formula_from_stoich(self.cursor[idx]['stoichiometry'], sort=True,
                                                            tex=False)] = (structure[-1], hull_dist[idx])
                                cache_misses += 1

            for idx, dist in enumerate(hull_dist):
                if np.abs(dist) < EPS:
                    hull_dist[idx] = 0

            failed_structures = []
            for ind, structure in enumerate(structures_finished):
                if not structure:
                    failed_structures.append(ind)
            if failed_structures:
                raise RuntimeError('There were issues calculating the hull distance for {} structures.'.format(len(failed_structures)))

            hull_dist = hull_dist[:-1]

        # otherwise, set to zero until proper N-d distance can be implemented
        else:
            raise NotImplementedError
            # self.hull.planes = [[self.structure_slice[vertex] for vertex in simplex]
                                # for simplex in self.hull.simplices]
            # for idx, structure in enumerate(structures):
                # if precompute and get_formula_from_stoich(self.cursor[idx]['stoichiometry'],
                                                          # tex=False) in cached_formula_dists:
                    # formula = get_formula_from_stoich(self.cursor[idx]['stoichiometry'], tex=False)
                    # if formula in cached_formula_dists:
                        # cache_hits += 1
                        # hull_dist[idx] = (structures[idx, -1] - cached_formula_dists[formula][0] +
                                          # cached_formula_dists[formula][1])
                    # structures_finished[idx] = True
                # else:
                    # if point_in_polygon(structure[:-1], plane):
                        # hull_dist[idx] = 0
                        # continue

            # loop over compositions
                # loop over planes
                    # check if point is interior to plane
                        # if so, compute "height" above plane, using equation of hyperplane from scipy

        return hull_dist

    def point_in_polygon(self):
        """ Use "ray-tracing" or winding number approach
        to compute whether a point lies inside a polygon.
        """
        raise NotImplementedError

    def construct_phase_diagram(self):
        """ Create a phase diagram with arbitrary chemical potentials.

        Expects self.chempot_cursor to be set with unique entries.

        """
        # only filter if hull is using non-elemental chempots to save time
        if self._non_elemental or self.from_cursor:
            self.cursor = self.filter_cursor_by_chempots(self.species, self.cursor)
            print('Cursor filtered down to {} structures.'.format(len(self.cursor)))
        else:
            for ind, doc in enumerate(self.cursor):
                self.cursor[ind]['concentration'] = get_concentration(doc, self.species)

        self.set_chempots()

        formation_key = 'formation_{}'.format(self.energy_key)
        for ind, doc in enumerate(self.cursor):
            self.cursor[ind][formation_key] = get_formation_energy(self.chempot_cursor, doc,
                                                                   energy_key=self.energy_key,
                                                                   temperature=self.temperature)

        # create stacked array of hull data
        structures = np.hstack((
            get_array_from_cursor(self.cursor, 'concentration').reshape(len(self.cursor), self._dimension-1),
            get_array_from_cursor(self.cursor, formation_key).reshape(len(self.cursor), 1)))

        if not self._non_elemental:
            self._setup_per_b_fields()
        elif self.args.get('subcmd') in ['voltage', 'volume']:
            raise NotImplementedError('Pseudo-binary/ternary voltages not yet implemented.')

        if self._dimension == 3:
            # add a point "above" the hull
            # for simple removal of extraneous vertices (e.g. top of 2D hull)
            dummy_point = [0.333, 0.333, 1e5]
            # if ternary, use all structures, not just those with negative eform for compatibility reasons
            self.structure_slice = np.vstack((structures, dummy_point))
        else:
            # filter out those with positive formation energy, to reduce expense computing hull
            self.structure_slice = structures[np.where(structures[:, -1] <= 0 + EPS)]

        # if we only have the chempots (or worse) with negative formation energy, don't even make the hull
        if len(self.structure_slice) <= self._dimension:
            if len(self.structure_slice) < self._dimension:
                raise RuntimeError('No chemical potentials on hull... either mysterious use of custom chempots, or worry!')
            self.convex_hull = FakeHull()
        else:
            try:
                self.convex_hull = ConvexHull(self.structure_slice)
            except QhullError:
                print('Error with QHull, plotting formation energies only...')
                print('To view errors, use --debug')
                if self.args.get('debug'):
                    print_exc()
                self.convex_hull = FakeHull()

        # remove vertices that have positive formation energy
        filtered_vertices = [vertex for vertex in self.convex_hull.vertices if self.structure_slice[vertex, -1] <= 0 + EPS]
        temp_simplices = self.convex_hull.simplices
        bad_simplices = []
        for ind, simplex in enumerate(temp_simplices):
            for vertex in simplex:
                if vertex not in filtered_vertices:
                    bad_simplices.append(ind)
                    break
        filtered_simplices = [simplex for ind, simplex in enumerate(temp_simplices) if ind not in bad_simplices]
        del self.convex_hull
        self.convex_hull = FakeHull()
        self.convex_hull.vertices = list(filtered_vertices)
        self.convex_hull.simplices = list(filtered_simplices)

        self.hull_dist = self.get_hull_distances(structures, precompute=True)
        set_cursor_from_array(self.cursor, self.hull_dist, 'hull_distance')

        # ensure hull cursor is sorted by enthalpy_per_atom,
        # then by concentration, as it will be by default if from database
        hull_cursor = [self.cursor[idx] for idx in np.where(self.hull_dist <= self.hull_cutoff + EPS)[0]]
        if self.temperature is not None:
            hull_cursor = sorted(hull_cursor, key=lambda doc: doc[self.energy_key][self.temperature])
        else:
            hull_cursor = sorted(hull_cursor, key=lambda doc: doc[self.energy_key])
        hull_cursor = sorted(hull_cursor, key=lambda k: k['concentration'])

        # if summary requested and we're in hulldiff mode, filter hull_cursor for lowest per stoich
        if self.args.get('summary') and self.args['subcmd'] == 'hulldiff':
            self.hull_cursor = []
            compositions = set()
            for ind, member in enumerate(hull_cursor):
                formula = get_formula_from_stoich(member['stoichiometry'], sort=True)
                if formula not in compositions:
                    compositions.add(formula)
                    self.hull_cursor.append(member)
        # otherwise, hull cursor includes all structures within hull_cutoff
        else:
            self.hull_cursor = hull_cursor

        self.structures = structures

    def set_chempots(self):
        """ Search for chemical potentials that match the structures in
        the query cursor and add them to the cursor.

        """
        query = self._query
        query_dict = dict()
        species_stoich = [sorted(get_stoich_from_formula(spec, sort=False)) for spec in self.species]
        self.chempot_cursor = []

        if self.args.get('chempots') is not None:
            self.fake_chempots(custom_elem=self.species)

        elif self.from_cursor:
            chempot_cursor = sorted([doc for doc in self.cursor if doc['stoichiometry'] in species_stoich],
                                    key=lambda k: (k[self.energy_key] if self.temperature is None
                                                   else k[self.energy_key][self.temperature]))

            for species in species_stoich:
                for doc in chempot_cursor:
                    if doc['stoichiometry'] == species:
                        self.chempot_cursor.append(doc)
                        break

            if len(self.chempot_cursor) != len(self.species):
                raise RuntimeError('Found {} of {} required chemical potentials'.format(len(self.chempot_cursor), len(self.species)))

            for ind, doc in enumerate(self.chempot_cursor):
                self.chempot_cursor[ind]['enthalpy_per_b'] = doc[self.energy_key]
                self.chempot_cursor[ind]['cell_volume_per_b'] = float('inf') if ind == 0 else 0.0

            if self.args.get('debug'):
                print([mu['stoichiometry'] for mu in self.chempot_cursor])

        else:
            print(60 * '─')
            self.chempot_cursor = len(self.species) * [None]
            # scan for suitable chem pots in database
            for ind, elem in enumerate(self.species):

                print('Scanning for suitable', elem, 'chemical potential...')
                query_dict['$and'] = deepcopy(list(query.calc_dict['$and']))

                if not self.args.get('ignore_warnings'):
                    query_dict['$and'].append(query._query_quality())

                if len(species_stoich[ind]) == 1:
                    query_dict['$and'].append(query._query_composition(custom_elem=[elem]))
                else:
                    query_dict['$and'].append(query._query_stoichiometry(custom_stoich=[elem]))

                # if oqmd, only query composition, not parameters
                if query.args.get('tags') is not None:
                    query_dict['$and'].append(query._query_tags())

                mu_cursor = query.repo.find(SON(query_dict)).sort(self.energy_key, pm.ASCENDING)
                if mu_cursor.count() == 0:
                    print_notify('Failed... searching without spin polarization field...')
                    scanned = False
                    while not scanned:
                        for idx, dicts in enumerate(query_dict['$and']):
                            for key in dicts:
                                if key == 'spin_polarized':
                                    del query_dict['$and'][idx][key]
                                    break
                            if idx == len(query_dict['$and']) - 1:
                                scanned = True
                    mu_cursor = query.repo.find(SON(query_dict)).sort(self.energy_key, pm.ASCENDING)

                if mu_cursor.count() == 0:
                    raise RuntimeError('No chemical potentials found for {}...'.format(elem))

                self.chempot_cursor[ind] = mu_cursor[0]
                if self.chempot_cursor[ind] is not None:
                    print('Using', ''.join([self.chempot_cursor[ind]['text_id'][0], ' ',
                                            self.chempot_cursor[ind]['text_id'][1]]), 'as chem pot for', elem)
                    print(60 * '─')
                else:
                    print_failure('No possible chem pots available for {}.'.format(elem))
                    raise RuntimeError('Exiting...')

            for i, mu in enumerate(self.chempot_cursor):
                self.chempot_cursor[i]['enthalpy_per_b'] = mu[self.energy_key]
                self.chempot_cursor[i]['num_a'] = 0

            self.chempot_cursor[0]['num_a'] = float('inf')

        # don't check for IDs if we're loading from cursor
        if not self.from_cursor:
            ids = [doc['_id'] for doc in self.cursor]
            if self.chempot_cursor[0]['_id'] not in ids:
                self.cursor.insert(0, self.chempot_cursor[0])
            for match in self.chempot_cursor[1:]:
                if match['_id'] not in ids:
                    self.cursor.append(match)

        # add faked chempots to overall cursor
        elif self.args.get('chempots') is not None:
            self.cursor.insert(0, self.chempot_cursor[0])
            for match in self.chempot_cursor[1:]:
                self.cursor.append(match)

        # find all elements present in the chemical potentials
        elements = []
        for mu in self.chempot_cursor:
            for elem, _ in mu['stoichiometry']:
                if elem not in elements:
                    elements.append(elem)
        self.elements = elements
        self.num_elements = len(elements)

    @staticmethod
    def filter_cursor_by_chempots(species, cursor):
        """ For the desired chemical potentials, remove any incompatible structures
        from cursor.

        Parameters:
            species (list): list of chemical potential formulae.
            cursor (list): list of matador documents to filter.

        Returns:
            list: the filtered cursor.
        """

        from matador.utils.cursor_utils import filter_cursor_by_chempots
        return filter_cursor_by_chempots(species, cursor)

    def _setup_per_b_fields(self):
        """ Calculate the enthalpy and volume per "B" in A_x B. """
        for ind, doc in enumerate(self.cursor):
            nums_b = len(self.species[1:]) * [0]
            for elem in doc['stoichiometry']:
                for chem_pot_ind, chem_pot in enumerate(self.species[1:]):
                    if elem[0] == chem_pot:
                        nums_b[chem_pot_ind] += elem[1]

            num_b = sum(nums_b)
            num_fu = doc['num_fu']
            if num_b == 0:
                self.cursor[ind]['enthalpy_per_b'] = 12345e5
                self.cursor[ind]['cell_volume_per_b'] = 12345e5
            else:
                if self.energy_key == 'enthalpy_per_atom':
                    self.cursor[ind]['enthalpy_per_b'] = doc['enthalpy'] / (num_b * num_fu)
                else:

                    if self.temperature is not None:
                        key = self.energy_key.split("_per_atom")[0]
                        self.cursor[ind]['enthalpy_per_b'] = doc[key][self.temperature] / (num_b * num_fu)
                    else:
                        self.cursor[ind]['enthalpy_per_b'] = doc[self.energy_key.split("_per_atom")[0]] / (num_b * num_fu)

                self.cursor[ind]['cell_volume_per_b'] = doc['cell_volume'] / (num_b * num_fu)

        Q = np.zeros((len(self.cursor)))
        for i, doc in enumerate(self.cursor):
            concs = self.cursor[i]['concentration']
            concs = get_padded_composition(doc['stoichiometry'], self.elements)
            Q[i] = get_generic_grav_capacity(concs, self.elements)
        set_cursor_from_array(self.cursor, Q, 'gravimetric_capacity')

    def voltage_curve(self, hull_cursor, quiet=False):
        """ Take a computed convex hull and calculate voltages for either binary or ternary
        systems. Sets the self.voltage_data attribute with various fields.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        Keyword arguments:
            quiet (bool): if False, print voltage data

        """
        if self._dimension == 2:
            self._calculate_binary_voltage_curve(hull_cursor, quiet=quiet)
        elif self._dimension == 3:
            self._calculate_ternary_voltage_curve(hull_cursor, quiet=quiet)
        else:
            raise RuntimeError('Unable to calculate voltage curve for hull of dimension {}'.format(self._dimension))

        data_str = ''
        for ind, path in enumerate(self.voltage_data['Q']):
            if ind != 0:
                data_str += '\n'
            if self._dimension == 3:
                data_str += '# ' + get_formula_from_stoich(self.voltage_data['endstoichs'][ind]) + '\n'
            else:
                data_str += '# ' + ''.join(self.species) + '\n'
            data_str += '# {:^10},\t{:^10},\t{:^10}\n'.format('x', 'Q (mAh/g)', 'Voltage (V)')
            for idx, _ in enumerate(path):
                data_str += '{:>10.2f},\t{:>10.2f},\t{:>10.4f}'.format(self.voltage_data['x'][ind][idx],
                                                                       self.voltage_data['Q'][ind][idx],
                                                                       self.voltage_data['voltages'][ind][idx])
                if idx != len(path) - 1:
                    data_str += '\n'
        if self.args.get('csv'):
            with open(''.join(self.species) + '_voltage.csv', 'w') as f:
                f.write(data_str)
        if not quiet:
            print('\nVoltage data:')
            print('\n' + data_str)

    def volume_curve(self, quiet=False):
        """ Take stable compositions and volume and calculate
        volume expansion per "B" in AB binary.

        Keyword arguments:
            quiet (bool): if False, print voltage data

        """

        if self._dimension == 2:
            self._calculate_binary_volume_curve()
        else:
            raise NotImplementedError('Volume curves have only been implemented for binary phase diagrams.')

        data_str = ''
        data_str += '# ' + ''.join(self.species) + '\n'
        data_str += '# {:>10},\t{:>10}\n'.format('x in {d[0]}x{d[1]}'.format(d=self.species),
                                                 'Volume per {} (Ang^3)'.format(self.species[1]))
        for idx, _ in enumerate(self.volume_data['x']):
            data_str += '{:>10.2f},\t{:>10.4f}'.format(self.volume_data['x'][idx],
                                                       self.volume_data['vol_per_y'][idx])
            if idx != len(self.volume_data['x']) - 1:
                data_str += '\n'
        if self.args.get('csv'):
            with open(''.join(self.species) + '_volume.csv', 'w') as f:
                f.write(data_str)
        if not quiet:
            print('\nVolume data:')
            print('\n' + data_str)

    def _calculate_binary_voltage_curve(self, hull_cursor, quiet=False):
        """ Generate binary voltage curve, setting the self.voltage_data
        dictionary.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        Keyword arguments:
            quiet (bool): if False, print voltage data

        """

        if not quiet:
            print('Generating voltage curve...')
        if self.temperature is not None and self.energy_key is not None:
            mu_enthalpy = []
            for doc in self.chempot_cursor:
                mu_enthalpy.append(doc[self.energy_key][self.temperature])
            mu_enthalpy = np.asarray(mu_enthalpy)
            if len(mu_enthalpy) != len(self.chempot_cursor):
                raise RuntimeError('Some keys were missing.')
        elif self.energy_key is not None:
            mu_enthalpy = get_array_from_cursor(self.chempot_cursor, self.energy_key)
        else:
            mu_enthalpy = get_array_from_cursor(self.chempot_cursor, 'enthalpy_per_atom')
        x = get_num_intercalated(hull_cursor)
        # sort for voltage calculation
        Q = get_array_from_cursor(hull_cursor, 'gravimetric_capacity')
        Q = Q[np.argsort(x)]
        stable_enthalpy_per_b = get_array_from_cursor(hull_cursor, 'enthalpy_per_b')[np.argsort(x)]
        x = np.sort(x)
        x, uniq_idxs = np.unique(x, return_index=True)
        stable_enthalpy_per_b = stable_enthalpy_per_b[uniq_idxs]
        Q = Q[uniq_idxs]
        V = []
        for i, _ in enumerate(x):
            V.append(
                -(stable_enthalpy_per_b[i] - stable_enthalpy_per_b[i-1]) / (x[i] - x[i-1]) + mu_enthalpy[0])
        V[0] = V[1]
        V[-1] = 0
        # make V, Q and x available for plotting
        self.voltage_data['voltages'] = []
        self.voltage_data['voltages'].append(V)
        self.voltage_data['Q'] = []
        self.voltage_data['Q'].append(Q)
        self.voltage_data['x'] = []
        self.voltage_data['x'].append(x)

    def _calculate_ternary_voltage_curve(self, hull_cursor, quiet=False):
        """ Calculate tenary voltage curve, setting self.voltage_data.
        First pass written by James Darby, jpd47@cam.ac.uk.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        Keyword arguments:
            quiet (bool): if False, print voltage data

        """
        points = np.hstack((
            get_array_from_cursor(hull_cursor, 'concentration'),
            get_array_from_cursor(hull_cursor, 'enthalpy_per_atom').reshape(len(hull_cursor), 1)))
        stoichs = get_array_from_cursor(hull_cursor, 'stoichiometry')
        mu_enthalpy = get_array_from_cursor(self.chempot_cursor, 'enthalpy_per_atom')
        enthalpy_active_ion = mu_enthalpy[0]
        # do another convex hull on just the known hull points, to allow access to useful indices
        hull = ConvexHull(points)

        endpoints = []
        endstoichs = []
        for ind, point in enumerate(points):
            if point[0] == 0 and point[1] != 0 and point[1] != 1:
                if not any([point.tolist() == test_point.tolist() for test_point in endpoints]):
                    endpoints.append(point)
                    endstoichs.append(stoichs[ind])
        if not quiet:
            print('{} starting point(s) found.'.format(len(endstoichs)))
            for endstoich in endstoichs:
                print(get_formula_from_stoich(endstoich), end=' ')
            print('\n')

        # iterate over possible endpoints of delithiation
        _reactions = []
        _voltages = []
        _Q = []
        for reaction_ind, endpoint in enumerate(endpoints):
            ratio = endpoint[1] / (1 - endpoint[0] - endpoint[1])
            if not quiet:
                print(30 * '-')
                print('Reaction {}, {}:'.format(reaction_ind, get_formula_from_stoich(endstoichs[reaction_ind])))
            y0 = endpoint[1] / (1 - endpoint[0])
            simp_in = 0
            intersections = []
            crossover = []
            # put starting point into crossover just in case it is not detected
            sum_conc = sum([int(species[1]) for species in endstoichs[reaction_ind]])
            conc = [int(species[1]) / sum_conc for species in endstoichs[reaction_ind]]
            if len(conc) == 2:
                conc.insert(0, 0.0)
            crossover.append(conc)
            for simplex in hull.simplices:
                tints = []
                for i in range(3):
                    j = (i + 1) % 3
                    e = points[simplex[i], 0]
                    f = points[simplex[i], 1]
                    g = points[simplex[j], 0] - points[simplex[i], 0]
                    h = points[simplex[j], 1] - points[simplex[i], 1]

                    x1 = e
                    y1 = f
                    z1 = 1 - x1 - y1
                    x2 = points[simplex[j], 0]
                    y2 = points[simplex[j], 1]
                    z2 = 1 - x2 - y2

                    if np.abs(h + g * y0) > EPS:
                        tin = (e * h + g * y0 - f * g) / (h + g * y0)
                        s2 = (y0 - e * y0 - f) / (h + g * y0)
                        if tin >= 0 and tin <= 1 and s2 >= 0 and s2 <= 1:
                            tints = np.append(tints, tin)
                            a = 1
                            # x1-x2 never == 0 on points we care about
                            if np.abs(x1 - x2) > EPS:
                                b = (y1 - y2) / (x1 - x2)
                                c = (z1 - z2) / (x1 - x2)
                                x_cross = tin
                                y_cross = b * (tin - x1) / a + y1
                                z_cross = c * (tin - x1) / a + z1
                                # only append unique points
                                if (len(crossover) == 0 or not np.any([np.isclose([x_cross, y_cross, z_cross], val)
                                                                       for val in crossover])):
                                    if y1 != 0 and y2 != 0 and round(float(z1 / y1), 5) == round(float(
                                            z2 / y2), 5) and round(float(z1 / y1), 5) == round(ratio, 5):
                                        pass
                                    else:
                                        crossover.append([x_cross, y_cross, z_cross])
                if len(tints) != 0:
                    temp = [simp_in, np.amin(tints), np.amax(tints)]
                    # condition removes the big triangle and the points which only graze the line of interest
                    if all([temp[2] > EPS, temp[1] < 1, temp[2] - temp[1] > EPS, temp[2] - temp[1] < 1,
                            temp[1] != temp[2]]):
                        intersections = np.append(intersections, temp)
                simp_in += 1

            # if tie line runs from fully de-lithiated to pure lithium (for example), then print and skip
            if not quiet:
                if len(intersections) == 0:
                    print_notify('No intermediate structures found for starting point {}.'.format(get_formula_from_stoich(endstoichs[reaction_ind])))
                    continue

            intersections = np.asarray(intersections)
            intersections = intersections.reshape(-1, 3)
            intersections = intersections[intersections[:, 1].argsort()]
            ends_of_rows = []
            min_values = []
            rows_to_keep = []
            # remove row corresponding to largest triangle, i.e. chempots only, and near duplicates (i.e. points with multiple tie-lines)
            for ind, row in enumerate(intersections):
                if not (row[1:].tolist() == [0, 1] or row[1:].tolist() in ends_of_rows or
                        np.any(np.isclose(row.tolist()[1], [val for val in min_values]))):
                    rows_to_keep.append(ind)
                    ends_of_rows.append(row[1:].tolist())
                    min_values.append(row.tolist()[1])
            intersections = intersections[rows_to_keep]

            voltages = []
            crossover = sorted(crossover)
            Q = sorted([get_generic_grav_capacity(point, self.species) for point in crossover])
            reactions = []
            reaction = [get_formula_from_stoich(endstoichs[reaction_ind])]
            reactions.append(reaction)
            for ind, face in enumerate(intersections):
                simplex_index = int(face[0])
                reaction = []
                reaction = [get_formula_from_stoich(hull_cursor[idx]['stoichiometry'])
                            for idx in hull.simplices[simplex_index]
                            if get_formula_from_stoich(hull_cursor[idx]['stoichiometry']) not in reaction]
                reactions.append(reaction)
                if not quiet:
                    print('{d[0]} + {d[1]} + {d[2]}'.format(d=reaction))
                Evec = points[hull.simplices[simplex_index], 2]
                Comp = points[hull.simplices[simplex_index], :]
                Comp[:, 2] = 1 - Comp[:, 0] - Comp[:, 1]

                Comp = Comp.T
                Compinv = np.linalg.inv(Comp)

                X = [1, 0, 0]

                V = -(Compinv.dot(X)).dot(Evec)
                V = V + enthalpy_active_ion
                # double up on first voltage
                if ind == 0:
                    voltages.append(V)
                if not quiet:
                    if ind != len(intersections) - 1:
                        print(5 * (ind + 1) * ' ' + ' ---> ', end='')
                voltages.append(V)

            _reactions.append(reactions)
            _Q.append(Q)
            _voltages.append(voltages)
            if not quiet:
                print('\n')
        assert len(_Q) == len(_voltages)

        self.voltage_data['x'] = _Q
        self.voltage_data['Q'] = _Q
        self.voltage_data['voltages'] = _voltages
        self.voltage_data['reactions'] = _reactions
        self.voltage_data['endstoichs'] = endstoichs

    def _calculate_binary_volume_curve(self):
        """ Take stable compositions and volume and calculate volume
        expansion per "B" in AB binary.

        """
        stable_comp = get_array_from_cursor(self.hull_cursor, 'concentration')
        for doc in self.hull_cursor:
            if 'cell_volume_per_b' not in doc:
                print(doc)
        stable_vol = get_array_from_cursor(self.hull_cursor, 'cell_volume_per_b')
        self.volume_data['x'] = np.asarray([comp / (1 - comp) for comp in stable_comp[:-1]]).flatten()
        self.volume_data['vol_per_y'] = np.asarray([vol for vol in stable_vol[:-1]])
        self.volume_data['bulk_volume'] = self.volume_data['vol_per_y'][0]
        self.volume_data['bulk_species'] = self.species[-1]
        self.volume_data['volume_ratio_with_bulk'] = self.volume_data['vol_per_y'] / self.volume_data['bulk_volume']
