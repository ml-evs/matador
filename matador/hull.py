# coding: utf-8
""" This file implements convex hull functionality
from database queries.
"""

from __future__ import print_function
# standard library
from traceback import print_exc
from bisect import bisect_left
from sys import exit
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
from matador.utils.chem_utils import get_formation_energy, get_concentration
from matador.utils.cursor_utils import set_cursor_from_array, get_array_from_cursor
from matador.utils.cursor_utils import display_results

EPS = 1e-12


class QueryConvexHull(object):
    """ Construct a binary or ternary phase diagram
    from matador.DBQuery object.
    """
    def __init__(self, query=None, cursor=None, elements=None, subcmd='hull', quiet=False, **kwargs):
        """ Initialise the class from either a DBQuery or a cursor (list of matador dicts)
        and construct the appropriate phase diagram.

        Args:

            | query    : matador.DBQuery, object containing structures,
            | cursor   : list(dict), alternatively specify list of matador docs.
            | elements : list(str), list of elements to use, used to provide a useful order,
            | subcmd   : either 'hull' or 'voltage',
            | kwargs   : mostly CLI arguments, see matador hull --help for full options.

        """
        self.args = kwargs
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd
        self.query = query
        self.from_cursor = False
        self.plot_param = False
        if self.query is not None:
            self.cursor = list(query.cursor)
            use_source = False
        else:
            self.cursor = cursor
            self.from_cursor = True
            use_source = True
        if self.cursor is None:
            raise RuntimeError('Failed to find structures to create hull!')
        if elements is None:
            elements = set()
            for doc in self.cursor:
                for species, _ in doc['stoichiometry']:
                    elements.add(species)
            self.elements = list(elements)
        else:
            self.elements = elements
            # filter out structures with any elements with missing chem pots
            self.cursor = [doc for doc in self.cursor if all([atom in self.elements for atom, num in doc['stoichiometry']])]

        if quiet:
            f = open(devnull, 'w')
            sys.stdout = f

        K2eV = 8.61733e-5
        if self.args.get('hull_temp') is not None:
            self.hull_cutoff = float(self.args['hull_temp']*K2eV)
        elif self.args.get('hull_cutoff') is not None:
            self.hull_cutoff = float(self.args['hull_cutoff'])
        else:
            self.hull_cutoff = 0.0

        if self.args.get('chempots') is not None:
            self.chem_pots = self.args.get('chempots')
            for ind, pot in enumerate(self.chem_pots):
                if pot > 0:
                    self.chem_pots[ind] = -1*self.chem_pots[ind]
        else:
            self.chem_pots = None

        self.hull_2d()

        if len(self.hull_cursor) == 0:
            print_warning('No structures on hull with chosen chemical potentials.')
        else:
            if self.args.get('hull_temp'):
                print_notify(str(len(self.hull_cursor)) + ' structures within ' +
                             str(self.args.get('hull_temp')) +
                             ' K of the hull with chosen chemical potentials.')
            else:
                print_notify(str(len(self.hull_cursor)) + ' structures within ' +
                             str(self.hull_cutoff) +
                             ' eV of the hull with chosen chemical potentials.')

        display_results(self.hull_cursor, self.args, hull=True, use_source=use_source)

        if not self.args.get('no_plot'):
            from matador import plotting
            self.set_plot_param()

        if self.args['subcmd'] == 'voltage':
            self.voltage_curve([doc for doc in self.hull_cursor if doc['hull_distance'] <= 1e-9])
            if not self.args.get('no_plot'):
                plotting.plot_voltage_curve(self)
                self.plot_hull()

        if self.args.get('volume'):
            self.volume_curve()
            if not self.args.get('no_plot'):
                plotting.plot_volume_curve(self, save=self.savefig)

        if self.args['subcmd'] == 'hull' and not self.args.get('no_plot'):
            if self.args.get('bokeh'):
                plotting.plot_2d_hull_bokeh(self)
            else:
                self.plot_hull()

        if not self.args.get('no_plot') and not self.savefig:
            import matplotlib.pyplot as plt
            plt.show()

        if quiet:
            f.close()
            sys.stdout = sys.__stdout__

    @property
    def savefig(self):
        return any([self.args.get('pdf'), self.args.get('png'), self.args.get('svg')])

    def plot_hull(self):
        """ Hull plot helper function. """
        from matador import plotting
        if self.ternary:
            plotting.plot_ternary_hull(self)
        else:
            plotting.plot_2d_hull(self)
        return

    def set_plot_param(self):
        """ Set some plotting options global to
        voltage and hull plots.
        """
        from matador import plotting
        import matplotlib.pyplot as plt

        if self.savefig:
            try:
                plt.style.use('article')
            except:
                print_exc()
                pass
        try:
            import seaborn as sns
            sns.set(font_scale=1.2)
            sns.set_style('ticks')
            sns.set_style({
                'axes.facecolor': 'white', 'figure.facecolor': 'white',
                'font.sans-serif': ['Linux Biolinum O', 'Helvetica', 'Arial'],
                'axes.linewidth': 0.5,
                'axes.grid': False,
                'legend.frameon': False,
                'axes.axisbelow': True})
        except:
            print_exc()
            pass
        self.scale = 1
        try:
            c = plt.cm.viridis(np.linspace(0, 1, 100))
            del c
            self.mpl_new_ver = True
        except:
            print_exc()
            self.mpl_new_ver = False
        Dark2_8 = plt.cm.get_cmap('Dark2').colors
        self.default_cmap_list = plotting.get_linear_cmap(Dark2_8[1:4], list_only=True)
        self.default_cmap = plotting.get_linear_cmap(Dark2_8[1:4], list_only=False)
        # first colour reserved for hull
        # penultimate colour reserved for off hull above cutoff
        # last colour reserved for OQMD
        Dark2_8_hex = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                       '#66a61e', '#e6ab02', '#a6761d', '#666666']
        self.colours = Dark2_8_hex
        self.colours.append('#bc80bd')
        self.plot_params = True
        return

    def get_chempots(self):
        """ Search for chemical potentials that match the structures in the query cursor,
        and add them to the cursor.
        """
        query = self.query
        self.mu_enthalpy = np.zeros((2))
        self.mu_volume = np.zeros((2))
        query_dict = dict()
        if not self.non_binary:
            elements = self.elements
        else:
            elements = self.chempot_search
        if self.chem_pots is not None:
            self.fake_chempots(custom_elem=elements)
        elif self.from_cursor:
            chempot_cursor = sorted([doc for doc in self.cursor if len(doc['stoichiometry']) == 1],
                                    key=lambda k: k['enthalpy_per_atom'])
            self.match = []
            for elem in elements:
                for doc in chempot_cursor:
                    if doc['stoichiometry'][0][0] == elem:
                        self.match.append(doc)
                        break
            if len(self.match) != len(elements):
                raise RuntimeError('Found {} of {} required chemical potentials'.format(len(self.match), len(elements)))
            for ind, doc in enumerate(self.match):
                self.match[ind]['hull_distance'] = 0
                self.match[ind]['enthalpy_per_b'] = doc['enthalpy_per_atom']
        else:
            print(60*'─')
            self.match = len(elements)*[None]
            # scan for suitable chem pots in database
            for ind, elem in enumerate(elements):
                print('Scanning for suitable', elem, 'chemical potential...')
                query_dict['$and'] = list(query.calc_dict['$and'])
                if not self.args.get('ignore_warnings'):
                    query_dict['$and'].append(query.query_quality())
                if not self.non_binary or ind == 0:
                    query_dict['$and'].append(query.query_composition(custom_elem=[elem]))
                else:
                    query_dict['$and'].append(query.query_stoichiometry(custom_stoich=[elem]))
                # if oqmd, only query composition, not parameters
                if query.args.get('tags') is not None:
                    query_dict['$and'].append(query.query_tags())
                if self.args.get('debug'):
                    print('Chemical potential query:')
                    from json import dumps
                    print(dumps(query_dict, indent=2))
                mu_cursor = query.repo.find(SON(query_dict)).sort('enthalpy_per_atom',
                                                                  pm.ASCENDING)
                if mu_cursor.count() == 0:
                    print('Failed... searching without spin polarization field...')
                    scanned = False
                    while not scanned:
                        for idx, dicts in enumerate(query_dict['$and']):
                            for key in dicts:
                                if key == 'spin_polarized':
                                    del query_dict['$and'][idx][key]
                                    break
                            if idx == len(query_dict['$and'])-1:
                                scanned = True
                    mu_cursor = query.repo.find(SON(query_dict)).sort('enthalpy_per_atom',
                                                                      pm.ASCENDING)

                try:
                    self.match[ind] = mu_cursor[0]
                except:
                    self.match[ind] = None
                if self.match[ind] is not None:
                    if ind == 0:
                        self.mu_enthalpy[ind] = float(self.match[ind]['enthalpy_per_atom'])
                        self.mu_volume[ind] = float(self.match[ind]['cell_volume'] /
                                                    self.match[ind]['num_atoms'])
                    else:
                        self.mu_enthalpy[1] += float(self.match[ind]['enthalpy_per_atom'])
                        self.mu_volume[1] = float(self.match[ind]['cell_volume'] /
                                                  self.match[ind]['num_atoms'])
                    print('Using', ''.join([self.match[ind]['text_id'][0], ' ',
                          self.match[ind]['text_id'][1]]), 'as chem pot for', elem)
                    print(60*'─')
                else:
                    print_failure('No possible chem pots found for ' + elem + '.')
                    exit()
            for i, mu in enumerate(self.match):
                self.match[i]['hull_distance'] = 0.0
                self.match[i]['enthalpy_per_b'] = mu['enthalpy_per_atom']
                self.match[i]['num_a'] = 0
            self.match[0]['num_a'] = float('inf')
        if not self.ternary and not self.from_cursor and not self.args.get('intersection'):
            self.cursor.insert(0, self.match[0])
            for match in self.match[1:]:
                self.cursor.append(match)
        return

    def fake_chempots(self, custom_elem=None):
        """ Spoof documents for command-line chemical potentials.

        Args:

            | custom_elem : list(str), list of element symbols to generate chempots for.

        """
        from matador.utils.chem_utils import get_stoich_from_formula
        self.match = [dict(), dict()]
        if custom_elem is None:
            custom_elem = self.elements
        for i, mu in enumerate(self.match):
            self.mu_enthalpy[i] = self.chem_pots[i]
            self.match[i]['enthalpy_per_atom'] = self.mu_enthalpy[i]
            self.match[i]['enthalpy'] = self.mu_enthalpy[i]
            self.match[i]['num_fu'] = 1
            self.match[i]['text_id'] = ['command', 'line']
            if self.non_binary and i == len(self.match)-1:
                self.match[i]['stoichiometry'] = get_stoich_from_formula(custom_elem)
            else:
                self.match[i]['stoichiometry'] = [[custom_elem[i], 1]]
            self.match[i]['space_group'] = 'xxx'
            self.match[i]['hull_distance'] = 0.0
            self.match[i]['enthalpy_per_b'] = self.match[i]['enthalpy_per_atom']
            self.match[i]['num_a'] = 0
            self.match[i]['cell_volume'] = 1
        self.match[0]['num_a'] = float('inf')
        notify = ('Using custom energies of ' + str(self.mu_enthalpy[0]) + ' eV/atom ' +
                  'and ' + str(self.mu_enthalpy[1]) + ' eV/atom as chemical potentials.')
        if self.args.get('debug'):
            for match in self.match:
                print(match)
        print(len(notify)*'─')
        print(notify)
        print(len(notify)*'─')

    def get_hull_distances(self, structures, precompute=True):
        """ Returns array of distances to pre-computed binary or ternary hull, from array
        containing concentrations and energies.

        Input:

            | structures : [N x n] np.ndarray, concentrations and enthalpies for N structures,
                           with up to 2 columns of concentrations and the last column containing
                           the structure's formation enthalpy.

        Args:

            | precompute: bool, whether or not to bootstrap hull distances from previously computed
                          values at the same stoichiometry.

        Returns:

            | hull_dist       : [N x 0] np.ndarray, distances to the hull for N structures,
            | tie_line_energy : [M x 1] np.ndarray, energies for structures on the precomputed hull,
                                sorted by concentration of the first element (the active ion).
            | tie_line_comp   : [M x 1] np.ndarray, sorted concentrations of first element in
                                structures on the precomputed hull.

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
                        hull_dist[ind] = structures[ind, -1] - cached_formula_dists[formula][0] + cached_formula_dists[formula][1]
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
            self.hull.planes = [[self.structure_slice[vertex] for vertex in simplex] for simplex in self.hull.simplices]
            structures_finished = [False]*len(structures)
            hull_dist = np.ones((len(structures)+1))
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
                    if precompute and get_formula_from_stoich(self.cursor[idx]['stoichiometry'], tex=False) in cached_formula_dists:
                        formula = get_formula_from_stoich(self.cursor[idx]['stoichiometry'], tex=False)
                        if formula in cached_formula_dists:
                            cache_hits += 1
                            hull_dist[idx] = structures[idx, -1] - cached_formula_dists[formula][0] + cached_formula_dists[formula][1]
                            structures_finished[idx] = True
                    else:
                        barycentric_structure = barycentric2cart(structure.reshape(1, 3)).T
                        barycentric_structure[-1, :] = 1
                        plane_barycentric_structure = np.matrix(planes_R_inv[ind]) * np.matrix(barycentric_structure)
                        if (plane_barycentric_structure >= 0-1e-12).all():
                            structures_finished[idx] = True
                            hull_dist[idx] = planes_height_fn[ind](structure)
                            if precompute:
                                cached_formula_dists[get_formula_from_stoich(self.cursor[idx]['stoichiometry'], tex=False)] = (structure[-1], hull_dist[idx])
                                cache_misses += 1
            self.failed_structures = []
            for ind in range(len(structures_finished)):
                if not structures_finished[ind]:
                    self.failed_structures.append(ind)
            self.failed_structures = np.asarray(self.failed_structures)
            hull_dist = hull_dist[:-1]

        # otherwise, set to zero until proper N-d distance can be implemented
        else:
            for ind in self.hull.vertices:
                hull_dist[ind] = 0.0

        return hull_dist, tie_line_energy, tie_line_comp

    def hull_2d(self, dis=False):
        """ Create a convex hull for a binary system. Sets several pieces of member data,
        most importantly self.hull and self.hull_cursor, as well as adding hull distances to
        self.cursor.
        """
        self.non_binary = False
        if self.query is not None:
            query = self.query
            self.elements = query.args.get('composition')
            if ':' in self.elements[0]:
                self.non_binary = True
                self.chempot_search = self.elements[0].split(':')
                if query.args.get('intersection'):
                    print_failure('Please disable intersection when creating a non-binary hull.')
                    exit()
            self.elements = [elem for elem in re.split(r'([A-Z][a-z]*)', self.elements[0]) if elem.isalpha()]
        assert(len(self.elements) < 4 and len(self.elements) > 1)
        self.ternary = False
        if len(self.elements) == 3 and not self.non_binary:
            self.ternary = True
        self.get_chempots()
        if self.non_binary:
            print('Contructing hull with non-elemental chemical potentials...')
        elif self.ternary:
            print('Constructing ternary hull...')
            if self.query is not None and not self.args.get('intersection'):
                print_warning('Please query with -int/--intersection when creating ternary hulls.')
                exit('Exiting...')
        else:
            print('Constructing binary hull...')
        # define hull by order in command-line arguments
        self.x_elem = [self.elements[0]]
        self.one_minus_x_elem = list(self.elements[1:])
        one_minus_x_elem = self.one_minus_x_elem
        # grab relevant information from query results; also make function?
        for ind, doc in enumerate(self.cursor):
            if not self.ternary:
                # calculate number of atoms of type B per formula unit
                nums_b = len(one_minus_x_elem)*[0]
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
                    self.cursor[ind]['enthalpy_per_b'] = doc['enthalpy'] / (num_b*num_fu)
                    self.cursor[ind]['cell_volume_per_b'] = doc['cell_volume'] / (num_b*num_fu)
            self.cursor[ind]['formation_enthalpy_per_atom'] = get_formation_energy(self.match, doc)
            self.cursor[ind]['concentration'] = get_concentration(doc, self.elements)
        # create stacked array of hull data
        structures = np.hstack((get_array_from_cursor(self.cursor, 'concentration'),
                                get_array_from_cursor(self.cursor, 'formation_enthalpy_per_atom').reshape(len(self.cursor), 1)))
        if not self.ternary and not self.non_binary:
            Q = get_binary_grav_capacities(get_num_intercalated(self.cursor), get_molar_mass(self.elements[1]))
            set_cursor_from_array(self.cursor, Q, 'gravimetric_capacity')
        else:
            Q = np.zeros((len(self.cursor)))
            for i in range(len(self.cursor)):
                concs = structures[i, 0:-1].tolist()
                concs.append(1-concs[0]-concs[1])
                Q[i] = get_generic_grav_capacity(concs, self.elements)
            set_cursor_from_array(self.cursor, Q, 'gravimetric_capacity')
        # create hull with SciPy routine, including only points with formation energy < 0
        if self.ternary:
            self.structure_slice = structures
            self.structure_slice = np.vstack((self.structure_slice, np.array([0, 0, 1e5])))
        elif self.non_binary:
            # if non-binary hull, remove middle concentration
            structures = structures[:, [0, -1]]
            self.structure_slice = structures[np.where(structures[:, -1] <= 0 + 1e-9)]
        else:
            self.structure_slice = structures[np.where(structures[:, -1] <= 0 + 1e-9)]
        if len(self.structure_slice) <= 2:
            if len(self.structure_slice) < 2:
                print_warning('No chemical potentials on hull... either mysterious use of custom chempots, or worry!')
            self.hull = FakeHull()
            self.hull_dist, self.hull_energy, self.hull_comp = self.get_hull_distances(structures)
            # should add chempots only to hull_cursor
            set_cursor_from_array(self.cursor, self.hull_dist, 'hull_distance')
        else:
            try:
                self.hull = ConvexHull(self.structure_slice)
            except QhullError:
                print_exc()
                print('Error with QHull, plotting points only...')
            # filter out top of hull - ugly
            if self.ternary:
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

            self.hull_dist, self.hull_energy, self.hull_comp = self.get_hull_distances(structures)
            set_cursor_from_array(self.cursor, self.hull_dist, 'hull_distance')

        # ensure hull cursor is sorted by enthalpy_per_atom, as it will be by default if from database
        hull_cursor = sorted([self.cursor[idx] for idx in np.where(self.hull_dist <= self.hull_cutoff + 1e-12)[0]],
                             key=lambda doc: doc['enthalpy_per_atom'])
        # if summary requested, filter for lowest per stoich
        if self.args.get('summary'):
            self.hull_cursor = []
            compositions = set()
            for ind, member in enumerate(hull_cursor):
                formula = get_formula_from_stoich(sorted(member['stoichiometry']))
                if formula not in compositions:
                    compositions.add(formula)
                    self.hull_cursor.append(member)
        else:
            self.hull_cursor = hull_cursor
        self.hull_cursor = sorted(self.hull_cursor, key=lambda k: k['concentration'])
        self.structures = structures

    def voltage_curve(self, hull_cursor):
        """ Take a computed convex hull and calculate voltages for either binary or ternary
        systems. Sets the self.x, self.Q and self.V member data for plotting.

        Input:

            | hull_cursor : list(dict), list of structures to include in the voltage curve.

        """

        if not self.ternary:
            print('Generating voltage curve...')
            mu_enthalpy = get_array_from_cursor(self.match, 'enthalpy_per_atom')
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
                V.append(-(stable_enthalpy_per_b[i] - stable_enthalpy_per_b[i-1]) /
                          (x[i] - x[i-1]) +
                          (mu_enthalpy[0]))
            V[0] = V[1]
            V[-1] = 0
            # make V, Q and x available for plotting
            self.voltages = []
            self.voltages.append(V)
            self.Q = []
            self.Q.append(Q)
            self.x = []
            self.x.append(x)

        elif self.ternary:
            """ Written by James Darby, jpd47@cam.ac.uk. """
            points = np.hstack((get_array_from_cursor(hull_cursor, 'concentration'),
                                get_array_from_cursor(hull_cursor, 'enthalpy_per_atom').reshape(len(hull_cursor), 1)))
            stoichs = get_array_from_cursor(hull_cursor, 'stoichiometry')
            mu_enthalpy = get_array_from_cursor(self.match, 'enthalpy_per_atom')
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
            print('{} starting point(s) found.'.format(len(endstoichs)))
            for endstoich in endstoichs:
                print(get_formula_from_stoich(endstoich), end=' ')
            print('\n')
            self.endstoichs = endstoichs

            # iterate over possible endpoints of delithiation
            self.voltages = []
            self.Q = []
            self.x = []
            for reaction_ind, endpoint in enumerate(endpoints):
                ratio = endpoint[1] / (1 - endpoint[0] - endpoint[1])
                print(30*'-')
                print('Reaction {}, {}:'.format(reaction_ind, get_formula_from_stoich(endstoichs[reaction_ind])))
                y0 = endpoint[1] / (1 - endpoint[0])
                simp_in = 0
                intersections = []
                crossover = []
                # put starting point into crossover just in case it is not detected
                sum_conc = sum([int(species[1]) for species in endstoichs[reaction_ind]])
                conc = [int(species[1])/sum_conc for species in endstoichs[reaction_ind]]
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

                        if np.abs(h + g*y0) > EPS:
                            tin = (e*h + g*y0 - f*g)/(h + g*y0)
                            s2 = (y0 - e*y0 - f) / (h + g*y0)
                            if tin >= 0 and tin <= 1 and s2 >= 0 and s2 <= 1:
                                tints = np.append(tints, tin)
                                a = 1
                                # x1-x2 never == 0 on points we care about
                                if np.abs(x1 - x2) > EPS:
                                    b = (y1 - y2)/(x1 - x2)
                                    c = (z1 - z2)/(x1 - x2)
                                    x_cross = tin
                                    y_cross = b * (tin-x1)/a + y1
                                    z_cross = c * (tin-x1)/a + z1
                                    # only append unique points
                                    if len(crossover) == 0 or not np.any([np.isclose([x_cross, y_cross, z_cross], crossover[i]) for i in range(len(crossover))]):
                                        if y1 != 0 and y2 != 0 and round(float(z1/y1), 5) == round(float(z2/y2), 5) and round(float(z1/y1), 5) == round(ratio, 5):
                                            pass
                                        else:
                                            crossover.append([x_cross, y_cross, z_cross])
                    if len(tints) != 0:
                        temp = [simp_in, np.amin(tints), np.amax(tints)]
                        # condition removes the big triangle and the points which only graze the line of interest
                        if all([temp[2] > EPS, temp[1] < 1, temp[2] - temp[1] > EPS, temp[2] - temp[1] < 1, temp[1] != temp[2]]):
                            intersections = np.append(intersections, temp)
                    simp_in += 1

                # if tie line runs from fully de-lithiated to pure lithium (for example), then print and skip
                if len(intersections) == 0:
                    print_warning('No intermediate structures found for starting point {}.'.format(get_formula_from_stoich(endstoichs[reaction_ind])))
                    continue

                intersections = np.asarray(intersections)
                intersections = intersections.reshape(-1, 3)
                intersections = intersections[intersections[:, 1].argsort()]
                ends_of_rows = []
                min_values = []
                rows_to_keep = []
                # remove row corresponding to largest triangle, i.e. chempots only, and near duplicates (i.e. points with multiple tie-lines)
                for ind, row in enumerate(intersections):
                    if not (row[1:].tolist() == [0, 1] or row[1:].tolist() in ends_of_rows or np.any(np.isclose(row.tolist()[1], [val for val in min_values]))):
                        rows_to_keep.append(ind)
                        ends_of_rows.append(row[1:].tolist())
                        min_values.append(row.tolist()[1])
                intersections = intersections[rows_to_keep]

                voltages = []
                crossover = sorted(crossover)
                Q = sorted([get_generic_grav_capacity(point, self.elements) for point in crossover])
                x = []
                reaction = [get_formula_from_stoich(endstoichs[reaction_ind])]
                for ind, face in enumerate(intersections):
                    simplex_index = int(face[0])
                    reaction = []
                    reaction = [get_formula_from_stoich(hull_cursor[idx]['stoichiometry'])
                                for idx in hull.simplices[simplex_index]
                                if get_formula_from_stoich(hull_cursor[idx]['stoichiometry']) not in reaction]
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
                    if ind != len(intersections)-1:
                        print(5*(ind+1)*' ' + ' ---> ', end='')
                    voltages.append(V)

                self.Q.append(Q)
                self.x.append(x)
                self.voltages.append(voltages)
                print('\n')
            assert len(self.Q) == len(self.voltages)

        print('Voltage data:')
        data_str = ''
        for ind, path in enumerate(self.Q):
            if ind != 0:
                data_str += '\n'
            if self.ternary:
                data_str += '# ' + get_formula_from_stoich(endstoichs[ind]) + '\n'
            else:
                data_str += '# ' + ''.join(self.elements) + '\n'
            data_str += '# {:>10},\t{:>10}\n'.format('Q (mAh/g)', 'Voltage (V)')
            for idx, _ in enumerate(path):
                data_str += '{:>10.2f},\t{:>10.4f}'.format(self.Q[ind][idx],
                                                           self.voltages[ind][idx])
                if idx != len(path) - 1:
                    data_str += '\n'
        if self.args.get('csv'):
            with open(''.join(self.elements) + '_voltage.csv', 'w') as f:
                f.write(data_str)
        print('\n' + data_str)

        return

    def volume_curve(self):
        """ Take stable compositions and volume and calculate
        volume expansion per "B" in AB binary.
        """
        stable_comp = get_array_from_cursor(self.hull_cursor, 'concentration')
        stable_vol = get_array_from_cursor(self.hull_cursor, 'cell_volume_per_b')
        # here, in A_x B_y
        # and v is the volume per x atom
        self.x = [comp/(1-comp) for comp in stable_comp]
        self.vol_per_y = [vol for vol in stable_vol]
        return
