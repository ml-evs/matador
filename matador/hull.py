# coding: utf-8
""" This file implements convex hull functionality
from database queries.
"""

from __future__ import print_function
# matador modules
from .utils.print_utils import print_failure, print_notify, print_warning
from .utils.hull_utils import barycentric2cart, vertices2plane, vertices2line
from .utils.chem_utils import get_binary_grav_capacities, get_molar_mass, get_num_intercalated
from .utils.chem_utils import get_generic_grav_capacity, get_formula_from_stoich
from .utils.chem_utils import get_formation_energy, get_concentration
from .utils.cursor_utils import set_cursor_from_array, get_array_from_cursor, filter_cursor
from .utils.cursor_utils import display_results
from .utils.glmol_wrapper import get_glmol_placeholder_string
from .export import generate_hash, generate_relevant_path
# external libraries
from scipy.spatial import ConvexHull
from bson.son import SON
import pymongo as pm
import numpy as np
# standard library
from traceback import print_exc
from bisect import bisect_left
from sys import exit
import re


class QueryConvexHull(object):
    """ Construct a binary or ternary phase diagram
    from matador.DBQuery object.
    """
    def __init__(self, query=None, cursor=None, elements=None, subcmd='hull', **kwargs):
        """

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
        if self.query is not None:
            self.cursor = list(query.cursor)
        else:
            self.cursor = cursor
            self.from_cursor = True
            for ind, doc in enumerate(self.cursor):
                self.cursor[ind]['text_id'] = ['xxx', 'yyy']
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

        display_results(self.hull_cursor, self.args, hull=True)

        if not self.args.get('no_plot'):
            self.set_plot_param()

        if self.args['subcmd'] == 'voltage':
            self.voltage_curve(self.hull_cursor)
            if not self.args.get('no_plot'):
                if self.args.get('subplot'):
                    self._subplot_voltage_hull()
                else:
                    self.plot_voltage_curve()
                self.plot_hull()

        if self.args.get('volume'):
            self.volume_curve()
            if not self.args.get('no_plot'):
                self.plot_volume_curve()

        if self.args['subcmd'] == 'hull' and not self.args.get('no_plot'):
            if self.args.get('bokeh'):
                self.plot_2d_hull_bokeh()
            else:
                if self.args.get('debug') and self.ternary:
                    self.plot_3d_ternary_hull()
                if self.ternary:
                    self.plot_ternary_hull()
                else:
                    self.plot_2d_hull()

        self.savefig = any([self.args.get('pdf'), self.args.get('png')])
        if not self.args.get('no_plot') and not self.savefig:
            import matplotlib.pyplot as plt
            plt.show()

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
            self.match = [doc for doc in self.cursor if len(doc['stoichiometry']) == 1 and doc['stoichiometry'][0][0] in self.elements]
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
        if not self.ternary and not self.from_cursor:
            self.cursor.insert(0, self.match[0])
            for match in self.match[1:]:
                self.cursor.append(match)
        return

    def fake_chempots(self, custom_elem=None):
        """ Spoof documents for command-line chemical potentials.

        Args:

            | custom_elem : list(str), list of element symbols to generate chempots for.

        """
        self.match = [dict(), dict()]
        if custom_elem is None:
            custom_elem = self.elements
        for i, mu in enumerate(self.match):
            self.mu_enthalpy[i] = self.chem_pots[i]
            self.match[i]['enthalpy_per_atom'] = self.mu_enthalpy[i]
            self.match[i]['enthalpy'] = self.mu_enthalpy[i]
            self.match[i]['num_fu'] = 1
            self.match[i]['text_id'] = ['command', 'line']
            # vomit-inducing cludge so that this works for custom chemical potentials for binarys that don't exist, provided
            # the first element has a two character symbol and the second has a one character symbol, e.g. TiP4...
            if self.non_binary and i == len(self.match)-1:
                self.match[i]['stoichiometry'] = [[custom_elem[i][:2], 1], [custom_elem[i][2:3], int(custom_elem[i][-1])]]
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

    def get_hull_distances(self, structures):
        """ Returns array of distances to pre-computed binary or ternary hull, from array
        containing concentrations and energies.

        Input:

            | structures : [N x n] np.ndarray, concentrations and enthalpies for N structures,
                           with up to 2 columns of concentrations and the last column containing
                           the structure's formation enthalpy.

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
        # if only chem pots on hull, dist = energy
        if len(self.structure_slice) == 2:
            hull_dist = np.ones((len(structures)))
            hull_dist = structures[:, -1]
        # if binary hull, do binary search
        elif len(self.structure_slice[0]) == 2:
            hull_dist = np.ones((len(structures)))
            for ind in range(len(structures)):
                i = bisect_left(tie_line_comp, structures[ind, 0])
                gradient, intercept = vertices2line([[tie_line_comp[i-1], tie_line_energy[i-1]],
                                                     [tie_line_comp[i], tie_line_energy[i]]])
                # calculate hull_dist
                hull_dist[ind] = structures[ind, -1] - (gradient * structures[ind, 0] + intercept)
        # otherwise, set to zero until proper N-d distance can be implemented
        else:
            # for each plane, convert each point into barycentric coordinates
            # for that plane and test for negative values
            self.hull.planes = [[self.structure_slice[vertex] for vertex in simplex] for simplex in self.hull.simplices]
            self.plane_points = []
            structures_sorted = [False]*len(structures)
            hull_dist = np.ones((len(structures)+1))
            for ind, plane in enumerate(self.hull.planes):
                self.plane_points.append([])
                R = barycentric2cart(plane).T
                R[-1, :] = 1
                # if projection of triangle in 2D is a line, do binary search
                if np.linalg.det(R) == 0:
                    if self.args.get('debug'):
                        print('TRANSFORMATION MATRIX IS SINGULAR')
                    continue
                else:
                    get_height_above_plane = vertices2plane(plane)
                    R_inv = np.linalg.inv(R)
                    for idx, structure in enumerate(structures):
                        if not structures_sorted[idx]:
                            barycentric_structure = barycentric2cart(structure.reshape(1, 3)).T
                            barycentric_structure[-1, :] = 1
                            plane_barycentric_structure = np.matrix(R_inv) * np.matrix(barycentric_structure)
                            if (plane_barycentric_structure >= 0-1e-12).all():
                                self.plane_points[-1].append(idx)
                                structures_sorted[idx] = True
                                hull_dist[idx] = get_height_above_plane(structure)

            # for ind in self.hull.vertices:
                # hull_dist[ind] = 0.0
            self.failed_structures = []
            for ind in range(len(structures_sorted)):
                if not structures_sorted[ind]:
                    self.failed_structures.append(ind)
            self.failed_structures = np.asarray(self.failed_structures)
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
                if self.ternary:
                    self.hull_dist = self.hull_dist[:-1]
                set_cursor_from_array(self.cursor, self.hull_dist, 'hull_distance')
            except:
                print_exc()
                print('Error with QHull, plotting points only...')

        hull_cursor = [self.cursor[idx] for idx in np.where(self.hull_dist <= self.hull_cutoff + 1e-12)[0]]
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
        self.structures = structures
        try:
            self.info = self.get_text_info(html=self.args.get('bokeh'))
            self.hull_info = self.get_text_info(cursor=self.cursor, hull=True, html=self.args.get('bokeh'))
        except:
            print_exc()
            pass

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

                        if h + g*y0 != 0:
                            tin = (e*h + g*y0 - f*g)/(h + g*y0)
                            s2 = (y0 - e*y0 - f) / (h + g*y0)
                            if tin >= 0 and tin <= 1 and s2 >= 0 and s2 <= 1:
                                tints = np.append(tints, tin)
                                a = 1
                                # x1-x2 never == 0 on points we care about
                                if x1 - x2 != 0:
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
                        # print('tints:', tints)
                        temp = [simp_in, np.amin(tints), np.amax(tints)]
                        # condition removes the big triangle and the points which only graze the line of interest
                        if temp[2] > 0 and temp[1] < 1 and temp[2] > temp[1] and temp[2] - temp[1] < 1:
                            intersections = np.append(intersections, temp)
                    simp_in += 1

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
                data_str += '{:>10.2f},\t{:>10.4f}'.format(self.Q[ind][idx], self.voltages[ind][idx])
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
        x = []
        # and v is the volume per x atom
        v = []
        for i in range(len(stable_comp)):
            x.append(stable_comp[i]/(1-stable_comp[i]))
            v.append(stable_vol[i])
        self.x = x
        self.vol_per_y = v
        return

    def plot_2d_hull(self, ax=None, dis=False, show=False, plot_points=True, plot_hull_points=True):
        """ Plot calculated hull, returning ax and fig objects for further editing.

        Args:

            | ax               : matplotlib axis object, an existing axis on which to plot,
            | show             : bool, whether or not to display the plot in an X window,
            | plot_points      : bool, whether or not to display off-hull structures.
            | plot_hull_points : bool, whether or not to display on-hull structures.

        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as colours
        if ax is None:
            if self.args.get('pdf') or self.args.get('png') or self.args.get('svg'):
                fig = plt.figure(facecolor=None, figsize=(8, 6))
            else:
                fig = plt.figure(facecolor=None)
            ax = fig.add_subplot(111)
        scatter = []
        x_elem = [self.elements[0]]
        one_minus_x_elem = list(self.elements[1:])
        tie_line = self.structure_slice[self.hull.vertices]
        plt.draw()
        # star structures on hull
        if len(self.structure_slice) != 2:
            if plot_hull_points:
                ax.scatter(tie_line[:, 0], tie_line[:, 1],
                           c=self.colours[1],
                           marker='o', zorder=99999, edgecolor='k',
                           s=self.scale*40, lw=1.5, alpha=1)
            ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1],
                    c=self.colours[0], lw=2, alpha=1, zorder=1000)
            if self.hull_cutoff > 0:
                ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1] + self.hull_cutoff,
                        '--', c=self.colours[1], lw=1, alpha=0.5, zorder=1000, label='')
            # annotate hull structures
            if self.args.get('labels'):
                ion = self.hull_cursor[0]['stoichiometry'][0][0]
                for elem in self.hull_cursor[1]['stoichiometry']:
                    if elem[0] == ion:
                        num_ion = elem[1]
                    else:
                        num_b = elem[1]
                ratio_A = num_ion / num_b
                for elem in self.hull_cursor[2]['stoichiometry']:
                    if elem[0] == ion:
                        num_ion = elem[1]
                    else:
                        num_b = elem[1]
                ratio_B = num_ion / num_b
                if ratio_A > ratio_B:
                    self.label_cursor = reversed(self.hull_cursor[1:-1])
                else:
                    self.label_cursor = self.hull_cursor[1:-1]
                for ind, doc in enumerate(self.label_cursor):
                    arrowprops = dict(arrowstyle="-|>", color='k')
                    if (ind+2) < np.argmin(tie_line[:, 1]):
                        position = (0.8*tie_line[ind+2, 0], 1.15*(tie_line[ind+2, 1])-0.05)
                    elif (ind+2) == np.argmin(tie_line[:, 1]):
                        position = (tie_line[ind+2, 0], 1.15*(tie_line[ind+2, 1])-0.05)
                    else:
                        position = (min(1.1*tie_line[ind+2, 0]+0.15, 0.95), 1.15*(tie_line[ind+2, 1])-0.05)
                    ax.annotate(get_formula_from_stoich(doc['stoichiometry'], elements=self.elements, tex=True),
                                xy=(tie_line[ind+2, 0], tie_line[ind+2, 1]),
                                xytext=position,
                                textcoords='data',
                                ha='right',
                                va='bottom',
                                arrowprops=arrowprops,
                                zorder=1)
            lw = self.scale * 0 if self.mpl_new_ver else 1
            # points for off hull structures
            if self.hull_cutoff == 0:
                # if no specified hull cutoff, ignore labels and colour
                # by distance from hull
                cmap_full = plt.cm.get_cmap('Dark2')
                cmap = colours.LinearSegmentedColormap.from_list(
                    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_full.name, a=0, b=1),
                    cmap_full(np.linspace(0.15, 0.4, 100)))
                if plot_points:
                    scatter = ax.scatter(self.structures[np.argsort(self.hull_dist), 0][::-1],
                                         self.structures[np.argsort(self.hull_dist), -1][::-1],
                                         s=self.scale*40, lw=lw, alpha=1, c=np.sort(self.hull_dist)[::-1],
                                         edgecolor='k', zorder=10000, cmap=cmap, norm=colours.LogNorm(0.02, 2))
                    cbar = plt.colorbar(scatter, aspect=30, pad=0.02, ticks=[0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                    cbar.ax.tick_params(length=0)
                    cbar.ax.set_yticklabels([0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                    cbar.ax.yaxis.set_ticks_position('right')
                    cbar.ax.set_frame_on(False)
                    cbar.outline.set_visible(False)
                    cbar.set_label('Distance from hull (eV)')
            if self.hull_cutoff != 0:
                # if specified hull cutoff, label and colour those below
                c = self.colours[1]
                for ind in range(len(self.structures)):
                    if self.hull_dist[ind] <= self.hull_cutoff or self.hull_cutoff == 0:
                        if plot_points:
                            scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1],
                                           s=self.scale*40, lw=lw, alpha=0.9, c=c, edgecolor='k',
                                           label=self.info[ind], zorder=300))
                if plot_points:
                    ax.scatter(self.structures[1:-1, 0], self.structures[1:-1, 1], s=self.scale*30, lw=lw,
                               alpha=0.3, c=self.colours[-2],
                               edgecolor='k', zorder=10)
            # tie lines
            ax.set_ylim(-0.1 if np.min(self.structure_slice[self.hull.vertices, 1]) > 0
                        else np.min(self.structure_slice[self.hull.vertices, 1])-0.15,
                        0.1 if np.max(self.structure_slice[self.hull.vertices, 1]) > 1
                        else np.max(self.structure_slice[self.hull.vertices, 1])+0.1)
        else:
            scatter = []
            print_exc()
            c = self.colours[1]
            lw = self.scale * 0 if self.mpl_new_ver else 1
            for ind in range(len(self.hull_cursor)):
                if type(self.hull_cursor[ind]['concentration']) is list:
                    if plot_points:
                        scatter.append(ax.scatter(self.hull_cursor[ind]['concentration'][0], self.hull_cursor[ind]['formation_enthalpy_per_atom'],
                                       s=self.scale*40, lw=1.5, alpha=1, c=c, edgecolor='k',
                                       zorder=1000))
                else:
                    if plot_points:
                        scatter.append(ax.scatter(self.hull_cursor[ind]['concentration'], self.hull_cursor[ind]['formation_enthalpy_per_atom'],
                                       s=self.scale*40, lw=1.5, alpha=1, c=c, edgecolor='k',
                                       zorder=1000))
                ax.plot([0, 1], [0, 0], lw=2, c=self.colours[0], zorder=900)
            for ind in range(len(self.structures)):
                if plot_points:
                    scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1],
                                   s=self.scale*40, lw=lw, alpha=0.9, c=c, edgecolor='k',
                                   zorder=300))

        if len(one_minus_x_elem) == 1:
            ax.set_title(x_elem[0] + '$_\mathrm{x}$' + one_minus_x_elem[0] + '$_\mathrm{1-x}$')
        if self.non_binary:
            ax.set_title(self.chempot_search[0] + '$_\mathrm{x}$(' + self.chempot_search[1] + ')$_\mathrm{1-x}$')
        plt.locator_params(nbins=3)
        ax.set_xlabel('x in {}$_\mathrm{{x}}${}$_\mathrm{{1-x}}$'.format(x_elem[0], one_minus_x_elem[0]))
        ax.grid(False)
        ax.set_xlim(-0.05, 1.05)
        ax.set_xticks([0, 0.33, 0.5, 0.66, 1])
        ax.set_xticklabels(ax.get_xticks())
        ax.set_yticks([0, -0.2, -0.4])
        ax.set_yticklabels(ax.get_yticks())
        ax.set_ylabel('Formation energy (eV/atom)')
        try:
            import seaborn as sns
            sns.despine(ax=ax, left=False, bottom=False)
        except:
            pass
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull.pdf',
                        dpi=500, bbox_inches='tight')
        if self.args.get('svg'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull.svg',
                        dpi=500, bbox_inches='tight')
        if self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull.png',
                        dpi=500, bbox_inches='tight')
        elif show:
            plt.show()

        return ax

    def plot_2d_hull_bokeh(self):
        """ Plot interactive hull with Bokeh. """
        from bokeh.plotting import figure, save, output_file
        from bokeh.models import ColumnDataSource, HoverTool, Range1d
        import matplotlib.pyplot as plt
        # grab tie-line structures
        tie_line_data = dict()
        tie_line_data['composition'] = list()
        tie_line_data['energy'] = list()
        for ind in range(len(self.hull.vertices)):
            if self.structure_slice[self.hull.vertices[ind], 1] <= 0:
                tie_line_data['composition'].append(self.structure_slice[self.hull.vertices[ind], 0])
                tie_line_data['energy'].append(self.structure_slice[self.hull.vertices[ind], 1])
        tie_line_data['energy'] = np.asarray(tie_line_data['energy'])
        tie_line_data['composition'] = np.asarray(tie_line_data['composition'])
        tie_line_data['energy'] = tie_line_data['energy'][np.argsort(tie_line_data['composition'])]
        tie_line_data['composition'] = np.sort(tie_line_data['composition'])

        # points for off hull structures
        hull_data = dict()
        hull_data['composition'] = self.structures[:, 0]
        hull_data['energy'] = self.structures[:, 1]
        hull_data['hull_distance'] = self.hull_dist
        hull_data['formula'], hull_data['text_id'] = [], []
        hull_data['space_group'], hull_data['hull_dist_string'] = [], []
        for structure in self.info:
            hull_data['formula'].append(structure[0])
            hull_data['text_id'].append(structure[1])
            hull_data['space_group'].append(structure[2])
            hull_data['hull_dist_string'].append(structure[3])
        cmap_limits = [0, 0.5]
        colormap = plt.cm.get_cmap('Dark2')
        cmap_input = np.interp(hull_data['hull_distance'], cmap_limits, [0.15, 0.4], left=0.15, right=0.4)
        colours = colormap(cmap_input, 1, True)
        bokeh_colours = ["#%02x%02x%02x" % (r, g, b) for r, g, b in colours[:, 0:3]]
        fixed_colours = colormap([0.0, 0.15], 1, True)
        tie_line_colour, on_hull_colour = ["#%02x%02x%02x" % (r, g, b) for r, g, b in fixed_colours[:, 0:3]]

        tie_line_source = ColumnDataSource(data=tie_line_data)
        hull_source = ColumnDataSource(data=hull_data)

        hover = HoverTool(tooltips="""
                          <div>
                              <div>
                                  <span style="font-size: 16px; font-family: "Fira Sans", sans-serif">
                                      Formula: @formula <br>
                                      ID: @text_id <br>
                                      Space group: @space_group <br>
                                      Distance from hull: @hull_dist_string
                                  </span>
                              </div>
                          </div>
                          """)

        tools = ['pan', 'wheel_zoom', 'reset', 'save']
        tools.append(hover)

        fig = figure(tools=tools)

        fig.xaxis.axis_label = 'x'
        fig.yaxis.axis_label = 'Formation energy (eV/atom)'
        fig.xaxis.axis_label_text_font_size = '20pt'
        fig.xaxis.axis_label_text_font = "Fira Sans, sans-serif"
        fig.yaxis.axis_label_text_font_size = '20pt'
        fig.yaxis.axis_label_text_font = "Fira Sans, sans-serif"
        fig.yaxis.axis_label_text_font_style = 'normal'
        fig.xaxis.axis_label_text_font_style = 'normal'
        fig.background_fill_alpha = 0
        fig.border_fill_alpha = 0
        fig.title.text_font_size = '20pt'
        fig.title.align = 'center'

        ylim = [-0.1 if np.min(self.structure_slice[self.hull.vertices, 1]) > 0
                else np.min(self.structure_slice[self.hull.vertices, 1])-0.15,
                0.1 if np.max(self.structure_slice[self.hull.vertices, 1]) > 1
                else np.max(self.structure_slice[self.hull.vertices, 1])+0.1]
        fig.x_range = Range1d(-0.1, 1.1)
        fig.y_range = Range1d(ylim[0], ylim[1])

        fig.line('composition', 'energy',
                 source=tie_line_source,
                 line_width=4,
                 line_color=tie_line_colour)
        hull_scatter = fig.scatter('composition', 'energy',
                                   source=hull_source,
                                   alpha=1,
                                   size=10,
                                   fill_color=bokeh_colours,
                                   line_color=None)
        fig.tools[0].renderers.append(hull_scatter)
        fig.square('composition', 'energy',
                   source=tie_line_source,
                   line_color='black',
                   color=on_hull_colour,
                   line_width=2,
                   alpha=1,
                   size=10)
        fig.plot_width = 800
        fig.plot_height = 600
        path = '/u/fs1/me388/data/hulls/'
        fname = generate_relevant_path(self.args) + '_' + generate_hash() + '.html'
        output_file(path+fname, title='Convex hull')
        print('Hull will be available shortly at http://www.tcm.phy.cam.ac.uk/~me388/hulls/' + fname)
        save(fig)
        glmol = False
        if glmol:
            html_string, js_string = get_glmol_placeholder_string()
            with open(path+fname) as f:
                flines = f.readlines()
                for ind, line in enumerate(flines):
                    if "<div class=\"bk-root\">" in line:
                        flines.insert(ind - 1, html_string)
                        break
                flines.append(js_string)
            with open(path+fname, 'w') as f:
                f.write('\n'.join(map(str, flines)))

    def plot_3d_ternary_hull(self):
        """ Plot calculated ternary hull in 3D. """
        from mpl_toolkits.mplot3d import axes3d
        # avoids annoying flake8 warning
        del axes3d
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        coords = barycentric2cart(self.structures)
        stable = coords[np.where(self.hull_dist < 0 + 1e-9)]
        stable = np.asarray(stable)
        ax.plot_trisurf(stable[:, 0], stable[:, 1], stable[:, 2], cmap=plt.cm.gnuplot, linewidth=1, color='grey', alpha=0.2)
        ax.scatter(stable[:, 0], stable[:, 1], stable[:, 2], s=100, c='k', marker='o')
        if len(self.failed_structures) > 0:
            ax.scatter(coords[self.failed_structures, 0], coords[self.failed_structures, 1], coords[self.failed_structures, 2], c='r')
        ax.set_zlim(-1, 1)
        ax.view_init(-90, 90)
        plt.show()

    def plot_hull(self):
        """ Hull plot helper function. """
        if self.ternary:
            self.plot_ternary_hull()
        else:
            self.plot_2d_hull()
        return

    def plot_ternary_hull(self, axis=None, show=False):
        """ Plot calculated ternary hull as a 2D projection.

        Takes optional matplotlib subplot axis as a parameter, and returns
        python-ternary subplot axis object.

        """
        import ternary
        import matplotlib.pyplot as plt
        import matplotlib
        import matplotlib.colors as colours
        try:
            import seaborn as sns
            sns.set_style({
                'axes.facecolor': 'white', 'figure.facecolor': 'white',
                'xtick.major.size': 0, 'xtick.minor.size': 0,
                'ytick.major.size': 0, 'ytick.minor.size': 0,
                'axes.linewidth': 0.0})
        except:
            print_exc()
            pass

        print('Plotting ternary hull...')
        if self.args.get('capmap') or self.args.get('efmap'):
            scale = 100
        elif self.args.get('sampmap'):
            scale = 20
        else:
            scale = 1
        fontsize = matplotlib.rcParams['font.size']

        if axis is not None:
            fig, ax = ternary.figure(scale=scale, ax=axis)
        else:
            fig, ax = ternary.figure(scale=scale)
        if self.args.get('capmap') or self.args.get('efmap') or self.args.get('sampmap'):
            fig.set_size_inches(8, 5)
        else:
            fig.set_size_inches(6.67, 5)
        ax.boundary(linewidth=2.0, zorder=99)
        ax.gridlines(color='black', multiple=scale*0.1, linewidth=0.5)

        ax.clear_matplotlib_ticks()
        if scale == 1:
            ax.ticks(axis='lbr', linewidth=1, multiple=scale*0.2, offset=0.02, fsize=fontsize-2)
        else:
            ax.ticks(axis='lbr', linewidth=1, multiple=scale*0.2, offset=0.02, fsize=fontsize-2,
                     ticks=[str(round(num, 1)) for num in np.linspace(0.0, 1.0, 6)])

        ax.set_title(''.join(self.elements), fontsize=fontsize+2, y=1.02)
        ax.left_axis_label(self.elements[2], fontsize=fontsize+2)
        ax.right_axis_label(self.elements[1], fontsize=fontsize+2)
        ax.bottom_axis_label(self.elements[0], fontsize=fontsize+2)

        concs = np.zeros((len(self.structures), 3))

        concs[:, :-1] = self.structures[:, :-1]
        for i in range(len(concs)):
            # set third triangular coordinate
            concs[i, -1] = 1 - concs[i, 0] - concs[i, 1]

        stable = np.asarray([concs[ind] for ind in self.hull.vertices])

        # sort by hull distances so things are plotting the right order
        concs = concs[np.argsort(self.hull_dist)].tolist()
        hull_dist = np.sort(self.hull_dist)

        filtered_concs = []
        filtered_hull_dists = []
        for ind, conc in enumerate(concs):
            if conc not in filtered_concs:
                if hull_dist[ind] <= self.hull_cutoff:
                    filtered_concs.append(conc)
                    filtered_hull_dists.append(hull_dist[ind])

        if self.hull_cutoff != 0:
            concs = np.asarray(filtered_concs)
            hull_dist = np.asarray(filtered_hull_dists)
            concs = concs[np.where(hull_dist <= self.hull_cutoff)]
            hull_dist = hull_dist[np.where(hull_dist <= self.hull_cutoff)]
        else:
            concs = np.asarray(concs)
            hull_dist = np.asarray(hull_dist)

        Ncolours = 1000
        min_cut = 0.01
        max_cut = 0.2
        min_colour = 0.15
        max_colour = 0.4
        colours_hull = plt.cm.Dark2(np.linspace(min_colour, max_colour, Ncolours))

        for plane in self.hull.planes:
            plane.append(plane[0])
            plane = np.asarray(plane)
            ax.plot(scale*plane, c=self.colours[0], lw=1.5, alpha=1, zorder=98)

        if self.args.get('pathways'):
            for phase in stable:
                if phase[0] == 0 and phase[1] != 0 and phase[2] != 0:
                    ax.plot([scale*phase, [scale, 0, 0]], c='r', alpha=0.2, lw=6, zorder=99)

        cmap_full = plt.cm.get_cmap('Dark2')
        cmap = colours.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_full.name, a=0, b=1),
            cmap_full(np.linspace(min_colour, max_colour, Ncolours)))

        colours_list = []
        colour_metric = hull_dist
        for i in range(len(colour_metric)):
            if colour_metric[i] >= max_cut:
                colours_list.append(Ncolours-1)
            elif colour_metric[i] <= min_cut:
                colours_list.append(0)
            else:
                colours_list.append(int((Ncolours-1)*(colour_metric[i] / max_cut)))
        colours_list = np.asarray(colours_list)
        ax.scatter(scale*concs, colormap=cmap, colorbar=True, cbarlabel='Distance from hull (eV/atom)',
                   c=hull_dist, vmax=max_cut, vmin=min_cut, zorder=1000, s=40, alpha=0)
        ax.scatter(scale*stable, marker='o', color=colours_hull[0], edgecolors='black', zorder=9999999,
                   s=150, lw=1.5)
        for i in range(len(concs)):
            ax.scatter(scale*concs[i].reshape(1, 3),
                       color=colours_hull[colours_list[i]],
                       marker='o',
                       zorder=10000-colours_list[i],
                       # alpha=max(0.1, 1-2*hull_dist[i]),
                       s=70*(1-float(colours_list[i])/Ncolours)+15,
                       lw=1, edgecolors='black')
        if self.args.get('capmap'):
            capacities = dict()
            from ternary.helpers import simplex_iterator
            for (i, j, k) in simplex_iterator(scale):
                capacities[(i, j, k)] = get_generic_grav_capacity([float(i)/scale, float(j)/scale, float(scale-i-j)/scale], self.elements)
            ax.heatmap(capacities, style="hexagonal", cbarlabel='Gravimetric capacity (maH/g)',
                       vmin=0, vmax=3000, cmap='Pastel2')
        elif self.args.get('efmap'):
            energies = dict()
            fake_structures = []
            from ternary.helpers import simplex_iterator
            for (i, j, k) in simplex_iterator(scale):
                fake_structures.append([float(i)/scale, float(j)/scale, 0.0])
            fake_structures = np.asarray(fake_structures)
            plane_energies, _, _ = self.get_hull_distances(fake_structures)
            ind = 0
            for (i, j, k) in simplex_iterator(scale):
                energies[(i, j, k)] = -1*plane_energies[ind]
                ind += 1
            ax.heatmap(energies, style="hexagonal", cbarlabel='Formation energy (eV/atom)',
                       vmax=0, cmap='bone')
        elif self.args.get('sampmap'):
            sampling = dict()
            from ternary.helpers import simplex_iterator
            eps = 1.0/float(scale)
            for (i, j, k) in simplex_iterator(scale):
                sampling[(i, j, k)] = np.size(np.where((concs[:, 0] <= float(i)/scale + eps) *
                                                       (concs[:, 0] >= float(i)/scale - eps) *
                                                       (concs[:, 1] <= float(j)/scale + eps) *
                                                       (concs[:, 1] >= float(j)/scale - eps) *
                                                       (concs[:, 2] <= float(k)/scale + eps) *
                                                       (concs[:, 2] >= float(k)/scale - eps)))
            ax.heatmap(sampling, style="hexagonal", cbarlabel='Number of structures',
                       cmap='afmhot')
        plt.tight_layout()
        if self.args.get('png'):
            plt.savefig(''.join(self.elements) + '.png', dpi=400, transparent=True, bbox_inches='tight')
        if self.args.get('pdf'):
            plt.savefig(''.join(self.elements) + '.pdf', dpi=400, transparent=True, bbox_inches='tight')
        elif show:
            ax.show()
        return ax

    def plot_voltage_curve(self, show=False):
        """ Plot calculated voltage curve. """
        import matplotlib.pyplot as plt
        if self.args.get('pdf') or self.args.get('png'):
            if len(self.voltages) != 1:
                fig = plt.figure(facecolor=None, figsize=(4, 3.5))
            else:
                fig = plt.figure(facecolor=None, figsize=(4, 3.5))
        else:
            fig = plt.figure(facecolor=None)
        axQ = fig.add_subplot(111)
        if self.args.get('expt') is not None:
            try:
                expt_data = np.loadtxt(self.args.get('expt'), delimiter=',')
            except:
                print_exc()
                pass
            if self.args.get('expt_label'):
                axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label=self.args.get('expt_label'))
            else:
                axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label='Experiment')
        for ind, voltage in enumerate(self.voltages):
            for i in range(len(voltage)-1):
                if i == 0 and self.args.get('expt'):
                    axQ.plot([self.Q[ind][i-1], self.Q[ind][i]], [voltage[i], voltage[i]], marker='*',
                             lw=2, c=self.colours[ind], label='DFT (this work)')
                elif i == 0 and len(self.voltages) != 1:
                    axQ.plot([self.Q[ind][i-1], self.Q[ind][i]], [voltage[i], voltage[i]], marker='o',
                             lw=2, c=self.colours[ind], label=get_formula_from_stoich(self.endstoichs[ind], tex=True))
                else:
                    axQ.plot([self.Q[ind][i-1], self.Q[ind][i]], [voltage[i], voltage[i]], marker='o',
                             lw=2, c=self.colours[ind])
                    if i != len(voltage)-2:
                        axQ.plot([self.Q[ind][i], self.Q[ind][i]], [voltage[i], voltage[i+1]], marker='o',
                                 lw=2, c=self.colours[ind])
        if self.args.get('labels'):
            ion = self.hull_cursor[0]['stoichiometry'][0][0]
            for elem in self.hull_cursor[1]['stoichiometry']:
                if elem[0] == ion:
                    num_ion = elem[1]
                else:
                    num_b = elem[1]
            ratio_A = num_ion / num_b
            for elem in self.hull_cursor[2]['stoichiometry']:
                if elem[0] == ion:
                    num_ion = elem[1]
                else:
                    num_b = elem[1]
            ratio_B = num_ion / num_b
            num_labels = len(self.hull_cursor)-2
            if ratio_A > ratio_B:
                self.label_cursor = list(reversed(self.hull_cursor[1:-1]))
            else:
                self.label_cursor = self.hull_cursor[1:-1]
            for i in range(num_labels):
                axQ.annotate(get_formula_from_stoich(self.label_cursor[i]['stoichiometry'], elements=self.elements, tex=True),
                             xy=(self.Q[0][i+1], self.voltages[0][i+1]+0.001),
                             textcoords='data',
                             ha='left',
                             zorder=9999)
        if self.args.get('expt') or len(self.voltages) != 1:
            axQ.legend(loc=1)
        axQ.set_ylabel('Voltage (V) vs {}$^+$/{}'.format(self.elements[0], self.elements[0]))
        axQ.set_xlabel('Gravimetric cap. (mAh/g)')
        start, end = axQ.get_ylim()
        axQ.set_ylim(0, 1.1*end)
        start, end = axQ.get_xlim()
        axQ.set_xlim(0, 1.1*end)
        axQ.grid('off')
        plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
        try:
            import seaborn as sns
            sns.despine()
            dark_grey = '#262626'
            for spine in ['left', 'bottom']:
                axQ.spines[spine].set_linewidth(0.5)
                axQ.spines[spine].set_color(dark_grey)
        except:
            pass
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_voltage.pdf',
                        dpi=500)
        elif self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_voltage.png',
                        dpi=500)
        elif show:
            plt.show()

    def plot_volume_curve(self, show=False):
        """ Plot calculated volume curve. """
        import matplotlib.pyplot as plt
        if self.args.get('pdf') or self.args.get('png'):
            fig = plt.figure(facecolor=None, figsize=(4, 3.5))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
        stable_hull_dist = get_array_from_cursor(self.hull_cursor, 'hull_distance')
        hull_vols = []
        hull_comps = []
        bulk_vol = self.vol_per_y[-1]
        for i in range(len(self.vol_per_y)):
            if stable_hull_dist[i] <= 0 + 1e-16:
                hull_vols.append(self.vol_per_y[i])
                hull_comps.append(self.x[i])
                s = 40
                zorder = 1000
                markeredgewidth = 1.5
                c = self.colours[1]
                alpha = 1
            else:
                s = 30
                zorder = 900
                alpha = 0.3
                markeredgewidth = 0
                c = 'grey'
            ax.scatter(self.x[i]/(1+self.x[i]), self.vol_per_y[i]/bulk_vol, marker='o', s=s, edgecolor='k', lw=markeredgewidth,
                       c=c, zorder=zorder, alpha=alpha)
        hull_comps, hull_vols = np.asarray(hull_comps), np.asarray(hull_vols)
        ax.plot(hull_comps/(1+hull_comps), hull_vols/bulk_vol, marker='o', lw=4,
                c=self.colours[0], zorder=100)
        ax.set_xlabel('$\mathrm{x}$ in $\mathrm{'+self.elements[0]+'_x'+self.elements[1]+'}_{1-x}$')
        ax.set_ylabel('Volume ratio with bulk')
        ax.set_ylim(0, 5*np.sort(hull_vols)[-2]/bulk_vol)
        ax.set_xlim(-0.05, 1.05)
        ax.yaxis.set_label_position('left')
        ax.set_xticklabels(ax.get_xticks())
        ax.grid('off')
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        tick_locs = np.linspace(0, 1, 6, endpoint=True).tolist()
        ax2.set_xticks(tick_locs)
        new_tick_labels = ['{}'.format(int(get_generic_grav_capacity([loc, 1-loc], [self.elements[0], self.elements[1]]))) for loc in tick_locs[:-1]]
        new_tick_labels.append('$\infty$')
        ax2.set_xlabel('Gravimetric capacity (mAh/g)')
        ax2.set_xticklabels(new_tick_labels)
        ax2.grid('off')
        try:
            import seaborn as sns
            sns.despine(top=False, right=False)
            dark_grey = '#262626'
            for spine in ['left', 'top', 'right', 'bottom']:
                ax.spines[spine].set_color(dark_grey)
                ax2.spines[spine].set_color(dark_grey)
                ax.spines[spine].set_linewidth(0.5)
                ax2.spines[spine].set_linewidth(0.5)
        except:
            pass
        # ax.yaxis.set_ticks(range(0, int(end)+1, 5))
        plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_volume.pdf',
                        dpi=300)
        if self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_volume.png',
                        dpi=300, bbox_inches='tight')
        elif show:
            plt.show()

    def _subplot_voltage_hull(self, dis=False):
        """ Plot calculated hull with inset voltage curve.

        DEPRECATED.

        """
        raise DeprecationWarning
        import matplotlib.pyplot as plt
        if self.args.get('pdf') or self.args.get('png'):
            fig = plt.figure(facecolor=None, figsize=(4.5, 1.5))
        else:
            fig = plt.figure(facecolor=None, figsize=(4.5, 1.5))
        ax = plt.subplot2grid((1, 3), (0, 0), colspan=2)
        ax2 = plt.subplot2grid((1, 3), (0, 2))
        scatter = []
        hull_scatter = []
        x_elem = self.elements[0]
        one_minus_x_elem = self.elements[1]
        plt.locator_params(nbins=3)
        # star structures on hull
        for ind in range(len(self.hull.vertices)):
            if self.structure_slice[self.hull.vertices[ind], 1] <= 0:
                hull_scatter.append(ax.scatter(self.structure_slice[self.hull.vertices[ind], 0],
                                               self.structure_slice[self.hull.vertices[ind], 1],
                                               c=self.colours[0],
                                               marker='*', zorder=99999, edgecolor='k',
                                               s=self.scale*150, lw=1, alpha=1,
                                               label=self.info[self.hull.vertices[ind]]))
        lw = self.scale * 0.05 if self.mpl_new_ver else 1
        # points for off hull structures
        for ind in range(len(self.structures)):
            if self.hull_dist[ind] <= self.hull_cutoff or self.hull_cutoff == 0:
                c = self.colours[self.source_ind[ind]] \
                    if self.hull_cutoff == 0 else self.colours[1]
                scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1],
                               s=self.scale*30, lw=lw, alpha=0.9, c=c, edgecolor='k',
                               label=self.info[ind], zorder=100))
        if self.hull_cutoff != 0:
            c = self.colours[self.source_ind[ind]] if self.hull_cutoff == 0 else self.colours[1]
            ax.scatter(self.structures[1:-1, 0], self.structures[1:-1, 1], s=self.scale*30, lw=lw,
                       alpha=0.3, c=self.colours[-2],
                       edgecolor='k', zorder=10)
        # tie lines
        for ind in range(len(self.hull_comp)-1):
            ax.plot([self.hull_comp[ind], self.hull_comp[ind+1]],
                    [self.hull_energy[ind], self.hull_energy[ind+1]],
                    c=self.colours[0], lw=2, alpha=1, zorder=1000, label='')
            if self.hull_cutoff > 0:
                ax.plot([self.hull_comp[ind], self.hull_comp[ind+1]],
                        [self.hull_energy[ind]+self.hull_cutoff,
                         self.hull_energy[ind+1]+self.hull_cutoff],
                        '--', c=self.colours[1], lw=1, alpha=0.5, zorder=1000, label='')
        ax.set_xlim(-0.05, 1.05)
        # data cursor
        ax.set_ylim(-0.1 if np.min(self.structure_slice[self.hull.vertices, 1]) > 0
                    else np.min(self.structure_slice[self.hull.vertices, 1]) - 0.05,
                    0.5 if np.max(self.structures[:, 1]) > 0.5
                    else np.max(self.structures[:, 1]) + 0.1)
        ax.set_title('$\mathrm{'+x_elem+'_x'+one_minus_x_elem+'_{1-x}}$')
        ax.set_xlabel('$x$', labelpad=-3)
        ax.set_xticks([0, 1])
        ax.set_yticks([-0.4, 0, 0.4])
        ax.set_ylabel('$E_\mathrm{F}$ (eV/atom)')
        # plot voltage
        for i in range(2, len(self.voltages)):
            ax2.scatter(self.x[i-1], self.voltages[i-1],
                        marker='*', s=100, edgecolor='k', c=self.colours[0], zorder=1000)
            ax2.plot([self.x[i], self.x[i]], [self.voltages[i], self.voltages[i-1]],
                     lw=2, c=self.colours[0])
            ax2.plot([self.x[i-1], self.x[i]], [self.voltages[i-1], self.voltages[i-1]],
                     lw=2, c=self.colours[0])
        ax2.set_ylabel('Voltage (V)')
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()
        ax2.set_xlim(0, np.max(np.asarray(self.x[1:]))+1)
        ax2.set_ylim(np.min(np.asarray(self.voltages[1:]))-0.1,
                     np.max(np.asarray(self.voltages[1:]))+0.1)
        ax2.set_xlabel('$n_\mathrm{Li}$', labelpad=-3)
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull_voltage.pdf',
                        dpi=300, bbox_inches='tight')
        elif self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull_voltage.png',
                        dpi=300, bbox_inches='tight')
        else:
            fig.show()

    def get_text_info(self, cursor=None, hull=False, html=False):
        """ Grab textual info for Bokeh plot labels. """
        info = []
        if cursor is None:
            cursor = self.cursor
        if hull:
            stoich_strings = []
        for ind, doc in enumerate(cursor):
            stoich_string = ''
            for elem in doc['stoichiometry']:
                stoich_string += elem[0]
                stoich_string += '$_{' + str(elem[1]) + '}$' if elem[1] != 1 else ''
            if hull:
                if stoich_string not in stoich_strings:
                    stoich_strings.append(stoich_string)
            info_string = "{0:^10}\n{1:^24}\n{2:^5s}\n{3:.3f} eV".format(stoich_string,
                                                                         doc['text_id'][0] + ' ' + doc['text_id'][1],
                                                                         doc['space_group'],
                                                                         doc['hull_distance'])
            if html:
                for char in ['$', '_', '{', '}']:
                    info_string = info_string.replace(char, '')
                info_string = info_string.split('\n')
            info.append(info_string)
        if hull:
            info = stoich_strings
        return info

    def set_plot_param(self):
        """ Set some plotting options global to
        voltage and hull plots.
        """
        import matplotlib.pyplot as plt
        try:
            plt.style.use('bmh')
        except:
            print_exc()
            pass
        if self.args.get('pdf') or self.args.get('png'):
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
        Dark2_8 = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a',
                   '#66a61e', '#e6ab02', '#a6761d', '#666666']
        # first colour reserved for hull
        # penultimate colour reserved for off hull above cutoff
        # last colour reserved for OQMD
        self.colours = Dark2_8
        self.colours.append('#bc80bd')
        return

    def generic_voltage_curve(self):
        """ A more generic version of voltage curve.

        DEPCRECATED.

        """
        raise DeprecationWarning
        import matplotlib.pyplot as plt

        def get_voltage_profile_segment(structure_new, structure_old,
                                        chempot_Li=self.mu_enthalpy[0]):
            """ Return voltage between two structures. """
            if structure_old['num_intercalated'] == float('inf'):
                V = 0
            else:
                V = (-(structure_new['enthalpy_per_b'] - structure_old['enthalpy_per_b']) /
                     (structure_new['num_intercalated'] - structure_old['num_intercalated']) +
                     (chempot_Li))
            return V

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.style.use('fivethirtyeight')

        mu_enthalpy = get_array_from_cursor(self.match, 'enthalpy_per_atom')
        x = get_num_intercalated(self.hull_cursor)
        stable_enthalpy_per_b = get_array_from_cursor(self.hull_cursor, 'enthalpy_per_b')[np.argsort(x)]
        set_cursor_from_array(self.hull_cursor, x, 'num_intercalated')
        for i in range(len(x)):
            for j in range(len(x)):
                if(self.hull_cursor[i]['gravimetric_capacity'] > self.hull_cursor[j]['gravimetric_capacity'] and
                   self.hull_cursor[i]['num_intercalated'] < self.hull_cursor[j]['num_intercalated']+2):
                    V = get_voltage_profile_segment(self.hull_cursor[i], self.hull_cursor[j], mu_enthalpy[0])
                    ax.plot([self.hull_cursor[i]['gravimetric_capacity'],
                            self.hull_cursor[j]['gravimetric_capacity']],
                            [V, V], alpha=0.05, c='b', lw=4)
                    ax.scatter((self.hull_cursor[i]['gravimetric_capacity'] + self.hull_cursor[j]['gravimetric_capacity']) / 2,
                               V, alpha=1, c='b', s=10, marker='o', zorder=10000000)
        self.hull_cursor = filter_cursor(self.hull_cursor, 'hull_distance', 0, 0.0001)
        x = get_num_intercalated(self.hull_cursor)
        Q = get_binary_grav_capacities(x, get_molar_mass(self.elements[1]))
        stable_enthalpy_per_b = get_array_from_cursor(self.hull_cursor, 'enthalpy_per_b')[np.argsort(x)]
        Q = Q[np.argsort(x)]
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
        # make V, Q and x available for plotting
        self.voltages = V
        self.Q = Q
        for i in range(len(self.voltages)-1):
            ax.plot([self.Q[i-1], self.Q[i]], [self.voltages[i], self.voltages[i]],
                    c='k', zorder=99999, lw=6)
            ax.plot([self.Q[i], self.Q[i]], [self.voltages[i], self.voltages[i+1]],
                    c='k', zorder=99999, lw=6)
        ax.set_ylabel('Voltage (V)')
        ax.set_xlabel('Gravimetric capacity (mAh/g)')
        plt.show()
        return


class FakeHull:
    """ Implements a thin class to mimic a ConvexHull object
    that would otherwise be undefined for two points. """
    def __init__(self):
        """ Define the used hull properties. """
        self.vertices = [0, 1]
        self.simplices = []
