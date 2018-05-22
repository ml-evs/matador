# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements convex hull functionality from database
queries.

"""


from traceback import print_exc
from bisect import bisect_left
import sys
import re
from os import devnull
# external libraries
from scipy.spatial import ConvexHull
from scipy.spatial.qhull import QhullError
from bson.son import SON
import pymongo as pm
import numpy as np
# matador modules
from matador.utils.print_utils import print_failure, print_notify, print_warning
from matador.utils.hull_utils import barycentric2cart, vertices2plane, vertices2line, FakeHull
from matador.utils.chem_utils import get_binary_grav_capacities, get_molar_mass, get_num_intercalated
from matador.utils.chem_utils import get_generic_grav_capacity, get_formula_from_stoich
from matador.utils.chem_utils import get_formation_energy, get_concentration, KELVIN_TO_EV
from matador.utils.cursor_utils import set_cursor_from_array, get_array_from_cursor
from matador.utils.cursor_utils import display_results

EPS = 1e-12


class QueryConvexHull(object):
    """ Construct a binary or ternary phase diagram from a
    matador.query.DBQuery object, or a list of structures.

    Attributes:
        cursor (list): list of all structures used to create phase diagram.
        hull_cursor (list): list of all documents within hull_cutoff.
        chempot_cursor (list): list of chemical potential documents.
        structure_slice (numpy.ndarray): array of concentrations and
            formation energies for each structure.
        hull_dist (np.ndarray): array of distances from hull for each structure.
        elements (list): list of chemical potential symbols.
        voltage_data (dict): if voltage_curve() has been called, then
            this is a dictionary containing Q, x, V and reaction pathways.
        voltage_data (dict): if volume_curve() has been called, then
            this is a dictionary containing x and volumes per B (in AxBy).

    """

    def __init__(self, query=None, cursor=None, elements=None, subcmd='hull', quiet=False,
                 plot_kwargs=None, **kwargs):
        """ Initialise the class from either a DBQuery or a cursor (list
        of matador dicts) and construct the appropriate phase diagram.

        Keyword arguments:
            query (matador.query.DBQuery): object containing structures,
            cursor (list(dict)): alternatively specify list of matador docs.
            elements (list(str)): list of elements to use, used to provide a useful order,
            subcmd (str): either 'hull' or 'voltage',
            kwargs (dict): mostly CLI arguments, see matador hull --help for full options.
            plot_kwargs (dict): arguments to pass to plot_hull function

        """
        self.args = kwargs
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd
        if plot_kwargs is None:
            plot_kwargs = {}
        self._query = query
        self.from_cursor = False
        self.plot_params = False
        if self._query is not None:
            self.cursor = list(query.cursor)
            use_source = False
        else:
            self.cursor = cursor
            self.from_cursor = True
            use_source = True

        # set up attributes for later
        self.chempot_cursor = None
        self.hull_cursor = None
        self.structure_slice = None
        self.hull_dist = None
        self.voltage_data = {}
        self.volume_data = {}

        if self.cursor is None:
            raise RuntimeError('Failed to find structures to create hull!')
        if elements is None:
            elements = set()
            for doc in self.cursor:
                for species, _ in doc['stoichiometry']:
                    elements.add(species)
            self.elements = list(elements)
        else:
            if isinstance(elements, str):
                self.elements = list(elements)
            else:
                self.elements = elements

            # filter out structures with any elements with missing chem pots
            self.cursor = [doc for doc in self.cursor if
                           all([atom in self.elements for atom, num in doc['stoichiometry']])]

        if quiet:
            cached_stdout = sys.stdout
            f = open(devnull, 'w')
            sys.stdout = f

        if self.args.get('energy_key') is not None:
            self._energy_key = self.args.get('energy_key')
        else:
            self._energy_key = 'enthalpy_per_atom'

        self.temperature = self.args.get('temperature')

        if self.args.get('hull_temp') is not None:
            self.hull_cutoff = float(self.args['hull_temp'] * KELVIN_TO_EV)
        elif self.args.get('hull_cutoff') is not None:
            self.hull_cutoff = float(self.args['hull_cutoff'])
        else:
            self.hull_cutoff = 0.0

        self.hull_2d()

        if not self.hull_cursor:
            print_warning('No structures on hull with chosen chemical potentials.')
        else:
            if self.args.get('hull_temp'):
                print_notify(
                    str(len(self.hull_cursor)) + ' structures within ' + str(self.args.get('hull_temp')) +
                    ' K of the hull with chosen chemical potentials.')
            else:
                print_notify(
                    str(len(self.hull_cursor)) + ' structures within ' + str(self.hull_cutoff) +
                    ' eV of the hull with chosen chemical potentials.')

        display_results(self.hull_cursor, self.args, hull=True, use_source=use_source)

        if not self.args.get('no_plot'):
            from matador import plotting

        if self.args['subcmd'] == 'voltage':
            self.voltage_curve([doc for doc in self.hull_cursor if doc['hull_distance'] <= 1e-9])
            if not self.args.get('no_plot'):
                plotting.plot_voltage_curve(self)
                self.plot_hull(**plot_kwargs)

        if self.args.get('volume'):
            self.volume_curve()
            if not self.args.get('no_plot'):
                plotting.plot_volume_curve(self)

        if self.args['subcmd'] == 'hull' and not self.args.get('no_plot'):
            self.plot_hull(**plot_kwargs)

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
        if self._ternary:
            plotting.plot_ternary_hull(self, **kwargs)
        else:
            plotting.plot_2d_hull(self, **kwargs)
        return

    def set_plot_param(self):
        """ Set some plotting options global to voltage and hull plots. """
        from matador import plotting
        import matplotlib.pyplot as plt

        if self.savefig:
            plt.style.use('article')
        import seaborn as sns
        sns.set(font_scale=1.2)
        sns.set_style('ticks')
        sns.set_style({
            'axes.facecolor': 'white', 'figure.facecolor': 'white',
            'axes.linewidth': 0.5,
            'axes.grid': False,
            'legend.frameon': False,
            'axes.axisbelow': True})
        dark2_8 = plt.cm.get_cmap('Dark2').colors
        self.default_cmap_list = plotting.get_linear_cmap(dark2_8[1:4], list_only=True)
        self.default_cmap = plotting.get_linear_cmap(dark2_8[1:4], list_only=False)
        # first colour reserved for hull
        # penultimate colour reserved for off hull above cutoff
        # last colour reserved for OQMD
        dark2_8_hex = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666']
        self.colours = dark2_8_hex
        self.colours.append('#bc80bd')
        self.plot_params = True
        return

    def get_chempots(self):
        """ Search for chemical potentials that match the structures in
        the query cursor and add them to the cursor.

        """
        query = self._query
        query_dict = dict()
        if not self._non_binary:
            elements = self.elements
        else:
            elements = self.chempot_search
        if self.args.get('chempots') is not None:
            self.fake_chempots(custom_elem=elements)
        elif self.from_cursor:
            chempot_cursor = sorted([doc for doc in self.cursor if len(doc['stoichiometry']) == 1],
                                    key=lambda k: k[self._energy_key])
            self.chempot_cursor = []
            for elem in elements:
                for doc in chempot_cursor:
                    if doc['stoichiometry'][0][0] == elem:
                        self.chempot_cursor.append(doc)
                        break
            if len(self.chempot_cursor) != len(elements):
                raise RuntimeError('Found {} of {} required chemical potentials'.format(len(self.chempot_cursor), len(elements)))
            for ind, doc in enumerate(self.chempot_cursor):
                self.chempot_cursor[ind]['hull_distance'] = 0
                self.chempot_cursor[ind]['enthalpy_per_b'] = doc[self._energy_key]
        else:
            print(60 * '─')
            self.chempot_cursor = len(elements) * [None]
            # scan for suitable chem pots in database
            for ind, elem in enumerate(elements):
                print('Scanning for suitable', elem, 'chemical potential...')
                from copy import deepcopy
                query_dict['$and'] = deepcopy(list(query.calc_dict['$and']))
                if not self.args.get('ignore_warnings'):
                    query_dict['$and'].append(query._query_quality())
                if not self._non_binary or ind == 0:
                    query_dict['$and'].append(query._query_composition(custom_elem=[elem]))
                else:
                    query_dict['$and'].append(query._query_stoichiometry(custom_stoich=[elem]))
                # if oqmd, only query composition, not parameters
                if query.args.get('tags') is not None:
                    query_dict['$and'].append(query._query_tags())
                mu_cursor = query.repo.find(SON(query_dict)).sort(self._energy_key, pm.ASCENDING)
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
                    mu_cursor = query.repo.find(SON(query_dict)).sort(self._energy_key, pm.ASCENDING)
                    if mu_cursor.count() == 0:
                        raise RuntimeError('No chemical potentials found...')

                self.chempot_cursor[ind] = mu_cursor[0]
                if self.chempot_cursor[ind] is not None:
                    print('Using', ''.join([self.chempot_cursor[ind]['text_id'][0], ' ',
                          self.chempot_cursor[ind]['text_id'][1]]), 'as chem pot for', elem)
                    print(60 * '─')
                else:
                    print_failure('No possible chem pots found for ' + elem + '.')
                    raise SystemExit('Exiting...')
            for i, mu in enumerate(self.chempot_cursor):
                self.chempot_cursor[i]['hull_distance'] = 0.0
                self.chempot_cursor[i]['enthalpy_per_b'] = mu[self._energy_key]
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

    def fake_chempots(self, custom_elem=None):
        """ Spoof documents for command-line chemical potentials.

        Keyword arguments:
            custom_elem (list(str)): list of element symbols to generate chempots for.

        """
        from matador.utils.chem_utils import get_stoich_from_formula
        from matador.export import generate_hash
        self.chempot_cursor = []

        if custom_elem is None:
            custom_elem = self.elements
        for i, _ in enumerate(self.args.get('chempots')):
            self.chempot_cursor.append(dict())
            self.chempot_cursor[i]['enthalpy_per_atom'] = -1*abs(self.args.get('chempots')[i])
            self.chempot_cursor[i]['enthalpy'] = -1*abs(self.args.get('chempots')[i])
            self.chempot_cursor[i]['num_fu'] = 1
            self.chempot_cursor[i]['text_id'] = ['command', 'line']
            self.chempot_cursor[i]['_id'] = generate_hash(hash_len=10)
            self.chempot_cursor[i]['source'] = ['command_line']
            if self._non_binary and i == len(self.chempot_cursor) - 1:
                self.chempot_cursor[i]['stoichiometry'] = get_stoich_from_formula(custom_elem)
            else:
                self.chempot_cursor[i]['atom_types'] = [custom_elem[i]]
                self.chempot_cursor[i]['stoichiometry'] = [[custom_elem[i], 1]]
            self.chempot_cursor[i]['space_group'] = 'xxx'
            self.chempot_cursor[i]['hull_distance'] = 0.0
            self.chempot_cursor[i]['enthalpy_per_b'] = self.chempot_cursor[i]['enthalpy_per_atom']
            self.chempot_cursor[i]['num_a'] = 0
            self.chempot_cursor[i]['cell_volume'] = 1
        self.chempot_cursor[0]['num_a'] = float('inf')
        notify = 'Custom chempots:'
        for chempot in self.chempot_cursor:
            notify += '{:3} = {} eV/atom, '.format(get_formula_from_stoich(chempot['stoichiometry']),
                                                   chempot['enthalpy'])

        if self.args.get('debug'):
            for match in self.chempot_cursor:
                print(match)
        print(len(notify) * '─')
        print(notify)
        print(len(notify) * '─')

    def get_hull_distances(self, structures, precompute=True):
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
        tie_line_comp = self.structure_slice[self.hull.vertices, 0]
        tie_line_energy = self.structure_slice[self.hull.vertices, -1]
        tie_line_comp = np.asarray(tie_line_comp)
        tie_line_energy = tie_line_energy[np.argsort(tie_line_comp)]
        tie_line_comp = tie_line_comp[np.argsort(tie_line_comp)]
        if precompute:
            # dict with formula keys, containing tuple of pre-computed enthalpy/atom and hull distance
            cached_formula_dists = dict()
            cache_hits = 0
            cache_misses = 0
        # if only chem pots on hull, dist = energy
        if len(self.structure_slice) == 2:
            hull_dist = np.ones((len(structures)))
            hull_dist = structures[:, -1]
        # if binary hull, do binary search
        elif len(self.structure_slice[0]) == 2:
            hull_dist = np.ones((len(structures)))
            if precompute:
                for ind, _ in enumerate(structures):
                    formula = get_formula_from_stoich(self.cursor[ind]['stoichiometry'], tex=False)
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
        elif len(self.structure_slice[0]) == 3:
            # for each plane, convert each point into barycentric coordinates
            # for that plane and test for negative values
            self.hull.planes = [[self.structure_slice[vertex] for vertex in simplex]
                                for simplex in self.hull.simplices]
            structures_finished = [False] * len(structures)
            hull_dist = np.ones((len(structures) + 1))
            planes_R_inv = []
            planes_height_fn = []
            for ind, plane in enumerate(self.hull.planes):
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
                for ind, plane in enumerate(self.hull.planes):
                    if structures_finished[idx] or planes_R_inv[ind] is None:
                        continue
                    if precompute and get_formula_from_stoich(self.cursor[idx]['stoichiometry'],
                                                              tex=False) in cached_formula_dists:
                        formula = get_formula_from_stoich(self.cursor[idx]['stoichiometry'], tex=False)
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
                                    get_formula_from_stoich(self.cursor[idx]['stoichiometry'],
                                                            tex=False)] = (structure[-1], hull_dist[idx])
                                cache_misses += 1

            failed_structures = []
            for ind in range(len(structures_finished)):
                if not structures_finished[ind]:
                    failed_structures.append(ind)
            if failed_structures:
                raise RuntimeError('There were issues calculating the hull distance for {} structures.'.format(len(failed_structures)))

            hull_dist = hull_dist[:-1]

        # otherwise, set to zero until proper N-d distance can be implemented
        else:
            for ind in self.hull.vertices:
                hull_dist[ind] = 0.0

        return hull_dist

    def hull_2d(self):
        """ Create a convex hull for a binary or ternary system. Sets
        several pieces of member data, most importantly self.hull and
        self.hull_cursor, as well as adding hull distances to self.cursor.

        """
        self._non_binary = False
        if self._query is not None:
            query = self._query
            elements_str = ''.join(query.args.get('composition'))
            if ':' in elements_str:
                self._non_binary = True
                self.chempot_search = elements_str.split(':')
                if query.args.get('intersection'):
                    print_failure('Please disable intersection when creating a non-binary hull.')
                    raise SystemExit('Exiting...')
            self.elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements_str) if elem.isalpha()]
        assert (len(self.elements) < 4 and len(self.elements) > 1)
        self._ternary = False
        if len(self.elements) == 3 and not self._non_binary:
            self._ternary = True
        self.get_chempots()
        if self._non_binary:
            print('Contructing hull with non-elemental chemical potentials...')
        elif self._ternary:
            print('Constructing ternary hull...')
            if self._query is not None and not self.args.get('intersection'):
                print_warning('Please query with -int/--intersection when creating ternary hulls.')
                raise SystemExit('Exiting...')
        else:
            print('Constructing binary hull...')
        # define hull by order in command-line arguments
        self.x_elem = [self.elements[0]]
        self.one_minus_x_elem = list(self.elements[1:])
        one_minus_x_elem = self.one_minus_x_elem
        formation_key = 'formation_' + self._energy_key
        # grab relevant information from query results; also make function?
        for ind, doc in enumerate(self.cursor):
            if not self._ternary:
                # calculate number of atoms of type B per formula unit
                nums_b = len(one_minus_x_elem) * [0]
                for elem in doc['stoichiometry']:
                    for chem_pot_ind, chem_pot in enumerate(one_minus_x_elem):
                        if elem[0] == chem_pot:
                            nums_b[chem_pot_ind] += elem[1]
                num_b = sum(nums_b)
                num_fu = doc['num_fu']
                # get enthalpy and volume per unit B: TODO - generalise this
                if num_b == 0:
                    self.cursor[ind]['enthalpy_per_b'] = 12345e5
                    self.cursor[ind]['cell_volume_per_b'] = 12345e5
                else:
                    self.cursor[ind]['enthalpy_per_b'] = doc['enthalpy'] / (num_b * num_fu)
                    self.cursor[ind]['cell_volume_per_b'] = doc['cell_volume'] / (num_b * num_fu)
            self.cursor[ind][formation_key] = get_formation_energy(self.chempot_cursor, doc, energy_key=self._energy_key, temperature=self.temperature)
            self.cursor[ind]['concentration'] = get_concentration(doc, self.elements)
        # create stacked array of hull data
        structures = np.hstack((
            get_array_from_cursor(self.cursor, 'concentration'),
            get_array_from_cursor(self.cursor, formation_key).reshape(len(self.cursor), 1)))
        if not self._ternary and not self._non_binary:
            Q = get_binary_grav_capacities(get_num_intercalated(self.cursor), get_molar_mass(self.elements[1]))
            set_cursor_from_array(self.cursor, Q, 'gravimetric_capacity')
        else:
            Q = np.zeros((len(self.cursor)))
            for i in range(len(self.cursor)):
                concs = structures[i, 0:-1].tolist()
                concs.append(1 - concs[0] - concs[1])
                Q[i] = get_generic_grav_capacity(concs, self.elements)
            set_cursor_from_array(self.cursor, Q, 'gravimetric_capacity')
        # create hull with SciPy routine, including only points with formation energy < 0
        if self._ternary:
            self.structure_slice = structures
            self.structure_slice = np.vstack((self.structure_slice, np.array([0, 0, 1e5])))
        elif self._non_binary:
            # if non-binary hull, remove middle concentration
            structures = structures[:, [0, -1]]
            self.structure_slice = structures[np.where(structures[:, -1] <= 0 + 1e-9)]
        else:
            self.structure_slice = structures[np.where(structures[:, -1] <= 0 + 1e-9)]

        if len(self.structure_slice) <= 2:
            if len(self.structure_slice) < 2:
                print_warning('No chemical potentials on hull... either mysterious use of custom chempots, or worry!')
            self.hull = FakeHull()

        else:
            try:
                self.hull = ConvexHull(self.structure_slice)
            except QhullError:
                print_exc()
                print('Error with QHull, plotting points only...')
            # filter out top of hull - ugly
            if self._ternary:
                filtered_vertices = [vertex for vertex in self.hull.vertices if self.structure_slice[vertex, -1] <= 0 + 1e-9]
                temp_simplices = self.hull.simplices
                bad_simplices = []
                for ind, simplex in enumerate(temp_simplices):
                    for vertex in simplex:
                        if vertex not in filtered_vertices:
                            bad_simplices.append(ind)
                            break
                filtered_simplices = [simplex for ind, simplex in enumerate(temp_simplices) if ind not in bad_simplices]
                del self.hull
                self.hull = FakeHull()
                self.hull.vertices = list(filtered_vertices)
                self.hull.simplices = list(filtered_simplices)

        self.hull_dist = self.get_hull_distances(structures)
        set_cursor_from_array(self.cursor, self.hull_dist, 'hull_distance')

        # ensure hull cursor is sorted by enthalpy_per_atom, then by concentration, as it will be by default if from database
        hull_cursor = [self.cursor[idx] for idx in np.where(self.hull_dist <= self.hull_cutoff + 1e-12)[0]]
        hull_cursor = sorted(hull_cursor, key=lambda doc: doc[self._energy_key])
        hull_cursor = sorted(hull_cursor, key=lambda k: k['concentration'])

        # if summary requested and we're in hulldiff mode, filter hull_cursor for lowest per stoich
        if self.args.get('summary') and self.args['subcmd'] == 'hulldiff':
            self.hull_cursor = []
            compositions = set()
            for ind, member in enumerate(hull_cursor):
                formula = get_formula_from_stoich(sorted(member['stoichiometry']))
                if formula not in compositions:
                    compositions.add(formula)
                    self.hull_cursor.append(member)

        # otherwise, hull cursor includes all structures within hull_cutoff
        else:
            self.hull_cursor = hull_cursor

        self.structures = structures

    def voltage_curve(self, hull_cursor, quiet=False):
        """ Take a computed convex hull and calculate voltages for either binary or ternary
        systems. Sets the self.voltage_data attribute with various fields.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        Keyword arguments:
            quiet (bool): if False, print voltage data

        """
        if not self._ternary:
            self._calculate_binary_voltage_curve(hull_cursor, quiet=quiet)
        elif self._ternary:
            self._calculate_ternary_voltage_curve(hull_cursor, quiet=quiet)

        data_str = ''
        for ind, path in enumerate(self.voltage_data['Q']):
            if ind != 0:
                data_str += '\n'
            if self._ternary:
                data_str += '# ' + get_formula_from_stoich(self.voltage_data['endstoichs'][ind]) + '\n'
            else:
                data_str += '# ' + ''.join(self.elements) + '\n'
            data_str += '# {:>10},\t{:>10}\n'.format('Q (mAh/g)', 'Voltage (V)')
            for idx, _ in enumerate(path):
                data_str += '{:>10.2f},\t{:>10.4f}'.format(self.voltage_data['Q'][ind][idx],
                                                           self.voltage_data['voltages'][ind][idx])
                if idx != len(path) - 1:
                    data_str += '\n'
        if self.args.get('csv'):
            with open(''.join(self.elements) + '_voltage.csv', 'w') as f:
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

        if not self._ternary and not self._non_binary:
            self._calculate_binary_volume_curve()
        else:
            raise NotImplementedError('Volume curves have only been implemented for binary phase diagrams.')

        data_str = ''
        data_str += '# ' + ''.join(self.elements) + '\n'
        data_str += '# {:>10},\t{:>10}\n'.format('x in {d[0]}x{d[1]}'.format(d=self.elements),
                                                 'Volume per {} (Ang^3)'.format(self.elements[1]))
        for idx, _ in enumerate(self.volume_data['x']):
            data_str += '{:>10.2f},\t{:>10.4f}'.format(self.volume_data['x'][idx],
                                                       self.volume_data['vol_per_y'][idx])
            if idx != len(self.volume_data['x']) - 1:
                data_str += '\n'
        if self.args.get('csv'):
            with open(''.join(self.elements) + '_volume.csv', 'w') as f:
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
        for i in range(len(x)):
            V.append(
                -(stable_enthalpy_per_b[i] - stable_enthalpy_per_b[i-1]) / (x[i] - x[i-1]) + (mu_enthalpy[0]))
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
        _x = []
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
            Q = sorted([get_generic_grav_capacity(point, self.elements) for point in crossover])
            x = []
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
            _x.append(x)
            _voltages.append(voltages)
            if not quiet:
                print('\n')
        assert len(_Q) == len(_voltages)

        self.voltage_data['x'] = _x
        self.voltage_data['Q'] = _Q
        self.voltage_data['voltages'] = _voltages
        self.voltage_data['reactions'] = _reactions
        self.voltage_data['endstoichs'] = endstoichs

    def _calculate_binary_volume_curve(self):
        """ Take stable compositions and volume and calculate volume
        expansion per "B" in AB binary.

        """
        stable_comp = get_array_from_cursor(self.hull_cursor, 'concentration')
        stable_vol = get_array_from_cursor(self.hull_cursor, 'cell_volume_per_b')
        self.volume_data['x'] = np.asarray([comp / (1 - comp) for comp in stable_comp[:-1]]).flatten()
        self.volume_data['vol_per_y'] = np.asarray([vol for vol in stable_vol[:-1]])
        self.volume_data['bulk_volume'] = self.volume_data['vol_per_y'][0]
        self.volume_data['bulk_species'] = self.elements[-1]
        self.volume_data['volume_ratio_with_bulk'] = self.volume_data['vol_per_y'] / self.volume_data['bulk_volume']
