# coding: utf-8
""" This file implements convex hull functionality
from database queries.
"""

from __future__ import print_function
# matador modules
from .utils.print_utils import print_failure, print_notify, print_warning
from .utils.hull_utils import barycentric2cart, vertices2plane, vertices2line
from .utils.chem_utils import get_binary_grav_capacities, get_molar_mass, get_num_intercalated
from .utils.chem_utils import get_generic_grav_capacity
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


class QueryConvexHull():
    """ Implements a Convex Hull for formation energies
    from a fryan DBQuery.
    """
    def __init__(self, query, subcmd='hull', **kwargs):
        """ Accept query from fryan as argument. """
        self.args = kwargs
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd
        self.query = query
        self.cursor = list(query.cursor)
        self.K2eV = 8.61733e-5

        if self.args.get('hull_temp') is not None:
            self.hull_cutoff = float(self.args['hull_temp']*self.K2eV)
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
            if self.args.get('debug'):
                self.generic_voltage_curve()
                # self.voltage_curve()
                # self.metastable_voltage_profile()
            else:
                self.voltage_curve()
                if self.args.get('subplot'):
                    self.subplot_voltage_hull()
                else:
                    self.plot_voltage_curve()
                self.plot_2d_hull()
        elif self.args.get('volume'):
            self.volume_curve()

        if self.args['subcmd'] == 'hull' and not self.args.get('no_plot'):
            if self.args.get('bokeh'):
                self.plot_2d_hull_bokeh()
            else:
                if self.args.get('debug'):
                    self.plot_3d_ternary_hull()
                if self.ternary:
                    self.plot_ternary_hull()
                else:
                    self.plot_2d_hull()
            if self.args.get('volume'):
                self.plot_volume_curve()

    def get_chempots(self):
        """ Search for chemical potentials that match
        the structures in the query cursor,
        and add them to the cursor.
        """
        query = self.query
        self.mu_enthalpy = np.zeros((2))
        self.mu_volume = np.zeros((2))
        query_dict = dict()
        if self.chem_pots is not None:
            self.fake_chempots()
        else:
            print(60*'─')
            if not self.non_binary:
                elements = self.elements
            else:
                elements = self.chempot_search
            self.match = len(elements)*[None]
            # scan for suitable chem pots in database
            for ind, elem in enumerate(elements):
                print('Scanning for suitable', elem, 'chemical potential...')
                query_dict['$and'] = list(query.calc_dict['$and'])
                if self.args.get('ignore_warnings') is None:
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
            self.cursor.insert(0, self.match[0])
            for match in self.match[1:]:
                self.cursor.append(match)
        return

    def fake_chempots(self):
        """ Spoof documents for command-line
        chemical potentials.
        """
        self.match = [dict(), dict()]
        for i, mu in enumerate(self.match):
            self.mu_enthalpy[i] = self.chem_pots[i]
            self.match[i]['enthalpy_per_atom'] = self.mu_enthalpy[i]
            self.match[i]['enthalpy'] = self.mu_enthalpy[i]
            self.match[i]['num_fu'] = 1
            self.match[i]['text_id'] = ['command', 'line']
            self.match[i]['stoichiometry'] = [[self.elements[i], 1]]
            self.match[i]['space_group'] = 'xxx'
            self.match[i]['hull_distance'] = 0.0
            self.match[i]['enthalpy_per_b'] = self.match[i]['enthalpy_per_atom']
            self.match[i]['num_a'] = 0
        self.match[0]['num_a'] = float('inf')
        notify = ('Using custom energies of ' + str(self.mu_enthalpy[0]) + ' eV/atom ' +
                  'and ' + str(self.mu_enthalpy[1]) + ' eV/atom as chemical potentials.')
        print(len(notify)*'─')
        print(notify)
        print(len(notify)*'─')

    def get_hull_distances(self, structures):
        """ Returns array of hull distances. """
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

    def get_text_info(self, cursor=None, hull=False, html=False):
        """ Grab textual info for plot labels. """
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

    def hull_2d(self, dis=False):
        """ Create a convex hull for two elements. """
        query = self.query
        self.elements = query.args.get('composition')
        self.non_binary = False
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
            if not self.args.get('intersection'):
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
        if len(self.structure_slice) == 2:
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

        self.hull_cursor = [self.cursor[idx] for idx in np.where(self.hull_dist <= self.hull_cutoff + 1e-12)[0]]
        self.structures = structures
        try:
            self.info = self.get_text_info(html=self.args.get('bokeh'))
            self.hull_info = self.get_text_info(cursor=self.cursor, hull=True, html=self.args.get('bokeh'))
        except:
            print_exc()
            pass

    def voltage_curve(self):
        """ Take convex hull and calculate voltages. """
        print('Generating voltage curve...')
        mu_enthalpy = get_array_from_cursor(self.match, 'enthalpy_per_atom')
        x = get_num_intercalated(self.hull_cursor)
        # sort for voltage calculation
        # Q = get_binary_grav_capacities(x, get_molar_mass(self.elements[1]))
        Q = get_array_from_cursor(self.hull_cursor, 'gravimetric_capacity')
        Q = Q[np.argsort(x)]
        stable_enthalpy_per_b = get_array_from_cursor(self.hull_cursor, 'enthalpy_per_b')[np.argsort(x)]
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
        self.x = x
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

    def plot_2d_hull(self, dis=False):
        """ Plot calculated hull. """
        import matplotlib.pyplot as plt
        import matplotlib.colors as colours
        if self.args.get('pdf') or self.args.get('png') or self.args.get('svg'):
            fig = plt.figure(facecolor=None, figsize=(5, 4))
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
            # for ind, doc in enumerate(tie_line):
                # ax.annotate(self.hull_info[ind],
                            # xy=(tie_line[ind, 0], tie_line[ind, 1]),
                            # xytext=(tie_line[ind, 0], tie_line[ind, 1] - 0.08),
                            # textcoords='data',
                            # ha='center',
                            # zorder=99999)
            lw = self.scale * 0 if self.mpl_new_ver else 1
            # points for off hull structures
            if self.hull_cutoff == 0:
                # if no specified hull cutoff, ignore labels and colour
                # by distance from hull
                cmap_full = plt.cm.get_cmap('Dark2')
                cmap = colours.LinearSegmentedColormap.from_list(
                    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_full.name, a=0, b=1),
                    cmap_full(np.linspace(0.15, 0.4, 100)))
                # scatter = ax.scatter(self.structures[:, 0], self.structures[:, 1],
                                     # s=self.scale*40, lw=lw, alpha=0.9,
                                     # c=self.hull_dist,
                                     # edgecolor='k', zorder=300, cmap=cmap)
                scatter = ax.scatter(self.structures[np.argsort(self.hull_dist), 0][::-1],
                                     self.structures[np.argsort(self.hull_dist), -1][::-1],
                                     s=self.scale*40, lw=lw, alpha=1, c=np.sort(self.hull_dist)[::-1],
                                     edgecolor='k', zorder=10000, cmap=cmap, norm=colours.LogNorm(0.02, 2))
                cbar = plt.colorbar(scatter, aspect=30, pad=0.02, ticks=[0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                cbar.ax.set_yticklabels([0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                cbar.set_label('Distance from hull (eV)')
            if self.hull_cutoff != 0:
                # if specified hull cutoff, label and colour those below
                c = self.colours[1]
                for ind in range(len(self.structures)):
                    if self.hull_dist[ind] <= self.hull_cutoff or self.hull_cutoff == 0:
                        scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1],
                                       s=self.scale*40, lw=lw, alpha=0.9, c=c, edgecolor='k',
                                       label=self.info[ind], zorder=300))
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
                    scatter.append(ax.scatter(self.hull_cursor[ind]['concentration'][0], self.hull_cursor[ind]['formation_enthalpy_per_atom'],
                                   s=self.scale*40, lw=1.5, alpha=1, c=c, edgecolor='k',
                                   zorder=1000))
                else:
                    scatter.append(ax.scatter(self.hull_cursor[ind]['concentration'], self.hull_cursor[ind]['formation_enthalpy_per_atom'],
                                   s=self.scale*40, lw=1.5, alpha=1, c=c, edgecolor='k',
                                   zorder=1000))
                ax.plot([0, 1], [0, 0], lw=2, c=self.colours[0], zorder=900)
            for ind in range(len(self.structures)):
                scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1],
                               s=self.scale*40, lw=lw, alpha=0.9, c=c, edgecolor='k',
                               zorder=300))

        if len(one_minus_x_elem) == 1:
            ax.set_title(x_elem[0] + '$_\mathrm{x}$' + one_minus_x_elem[0] + '$_\mathrm{1-x}$')
        if self.non_binary:
            ax.set_title(self.chempot_search[0] + '$_\mathrm{x}$(' + self.chempot_search[1] + ')$_\mathrm{1-x}$')
        plt.locator_params(nbins=3)
        ax.set_xlabel('$\mathrm{x}$')
        ax.grid(False)
        ax.set_xlim(-0.05, 1.05)
        ax.set_xticks([0, 0.33, 0.5, 0.66, 1])
        ax.set_xticklabels(ax.get_xticks())
        ax.set_ylabel('E$_\mathrm{F}$ (eV/atom)')
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull.pdf',
                        dpi=400, bbox_inches='tight')
        elif self.args.get('svg'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull.svg',
                        dpi=200, bbox_inches='tight')
        elif self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull.png',
                        dpi=200, bbox_inches='tight')
        else:
            plt.show()

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

    def plot_ternary_hull(self):
        """ Plot calculated ternary hull as a 2D projection. """
        import ternary
        import matplotlib.pyplot as plt
        import matplotlib.colors as colours
        print('Plotting ternary hull...')
        if self.args.get('capmap') or self.args.get('efmap'):
            scale = 100
        elif self.args.get('sampmap'):
            scale = 20
        else:
            scale = 1
        fontsize = 18
        fig, ax = ternary.figure(scale=scale)
        fig.set_size_inches(10, 7.5)

        ax.boundary(linewidth=2.0, zorder=99)
        ax.gridlines(color='black', multiple=scale*0.1, linewidth=0.5)

        ax.clear_matplotlib_ticks()
        if scale == 1:
            ax.ticks(axis='lbr', linewidth=1, multiple=scale*0.1)
        else:
            ax.ticks(axis='lbr', linewidth=1, multiple=scale*0.1,
                     ticks=[str(round(num, 1)) for num in np.linspace(0.0, 1.0, 11)])

        ax.set_title(''.join(self.elements), fontsize=fontsize)
        ax.left_axis_label(self.elements[2], fontsize=fontsize)
        ax.right_axis_label(self.elements[1], fontsize=fontsize)
        ax.bottom_axis_label(self.elements[0], fontsize=fontsize)

        concs = np.zeros((len(self.structures), 3))

        concs[:, :-1] = self.structures[:, :-1]
        for i in range(len(concs)):
            # set third triangular coordinate
            concs[i, -1] = 1 - concs[i, 0] - concs[i, 1]

        stable = np.asarray([concs[ind] for ind in self.hull.vertices])

        # sort by hull distances so things are plotting the right order
        concs = concs[np.argsort(self.hull_dist)]
        hull_dist = np.sort(self.hull_dist)

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
                    ax.plot([scale*phase, [scale, 0, 0]], c='r', alpha=0.2, lw=6, zorder=100000000)

        cmap_full = plt.cm.get_cmap('Dark2')
        cmap = colours.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_full.name, a=0, b=1),
            # cmap_full(np.logspace(-0.8239, -0.3979, 100)))
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
        ax.scatter(scale*concs, colormap=cmap, colorbar=True,
                   c=hull_dist, vmax=max_cut, vmin=min_cut, zorder=1000, s=40, alpha=0)
        ax.scatter(scale*stable, marker='o', color=colours_hull[0], edgecolors='black', zorder=9999999,
                   s=150, lw=1.5)
        for i in range(len(concs)):
            ax.scatter(scale*concs[i].reshape(1, 3),
                       color=colours_hull[colours_list[i]],
                       marker='s',
                       zorder=10000-colours_list[i],
                       alpha=1,
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
        if self.args.get('png'):
            plt.savefig('ternary.png', dpi=400)
        elif self.args.get('pdf'):
            plt.savefig('ternary.pdf', dpi=400)
        ax.show()

    def plot_voltage_curve(self):
        """ Plot calculated voltage curve. """
        import matplotlib.pyplot as plt
        if self.args.get('pdf') or self.args.get('png'):
            fig = plt.figure(facecolor=None, figsize=(3, 2.7))
        else:
            fig = plt.figure(facecolor=None)
        axQ = fig.add_subplot(111)
        # axQ = ax.twiny()
        if self.args.get('expt') is not None:
            try:
                expt_data = np.loadtxt(self.args.get('expt'), delimiter=',')
                axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', ls='--')
            except:
                print_exc()
                pass
        for i in range(len(self.voltages)-1):
            axQ.plot([self.Q[i-1], self.Q[i]], [self.voltages[i], self.voltages[i]],
                     lw=2, c=self.colours[0])
            axQ.plot([self.Q[i], self.Q[i]], [self.voltages[i], self.voltages[i+1]],
                     lw=2, c=self.colours[0])
        for i in range(len(self.x)):
            if self.x[i] < 1e9:
                string_stoich = ''
                if abs(np.ceil(self.x[i])-self.x[i]) > 1e-8:
                    string_stoich = str(round(self.x[i], 1))
                else:
                    string_stoich = str(int(np.ceil(self.x[i])))
                if string_stoich is '1':
                    string_stoich = ''
                if string_stoich is '0':
                    string_stoich = ''
                else:
                    string_stoich = self.elements[0] + '$_{' + string_stoich + '}$' + self.elements[1]
                # axQ.annotate(string_stoich,
                             # xy=(self.Q[i], self.voltages[i]+0.001),
                             # textcoords='data',
                             # ha='center',
                             # zorder=9999)
        axQ.set_ylabel('Voltage (V)')
        axQ.set_xlabel('Gravimetric cap. (mAh/g)')
        start, end = axQ.get_ylim()
        axQ.set_ylim(start-0.01, end+0.01)
        axQ.grid('off')
        plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_voltage.pdf',
                        dpi=300)
        elif self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_voltage.png',
                        dpi=300, bbox_inches='tight')
        else:
            plt.show()

    def plot_volume_curve(self):
        """ Plot calculate volume curve. """
        import matplotlib.pyplot as plt
        if self.args.get('pdf') or self.args.get('png'):
            fig = plt.figure(facecolor=None, figsize=(2.7, 2.7))
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
            else:
                s = 30
                zorder = 900
                markeredgewidth = 0.5
                c = 'grey'
            ax.scatter(self.x[i], self.vol_per_y[i]/bulk_vol, marker='o', s=s, edgecolor='k', lw=markeredgewidth,
                       c=c, zorder=zorder)
        hull_comps, hull_vols = np.asarray(hull_comps), np.asarray(hull_vols)
        ax.plot(hull_comps, hull_vols/bulk_vol, marker='o', lw=4,
                c=self.colours[0], zorder=100)
        ax.set_xlabel('$\mathrm{u}$ in $\mathrm{'+self.elements[0]+'_u'+self.elements[1]+'}$')
        ax.set_ylabel('Volume ratio with bulk')
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(range(0, int(end)+1, 1))
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(range(1, int(end)+1, 1))
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.set_xticklabels(ax.get_xticks())
        ax.grid('off')
        plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_volume.pdf',
                        dpi=300)
        elif self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_volume.png',
                        dpi=300, bbox_inches='tight')
        else:
            plt.show()

    def subplot_voltage_hull(self, dis=False):
        """ Plot calculated hull with inset voltage curve. """
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

    def set_plot_param(self):
        """ Set some plotting options global to
        voltage and hull plots.
        """
        import matplotlib.pyplot as plt
        try:
            plt.style.use('bmh')
        except:
            pass
        if self.args.get('pdf') or self.args.get('png'):
            try:
                plt.style.use('article')
            except:
                pass
        self.scale = 1
        try:
            c = plt.cm.viridis(np.linspace(0, 1, 100))
            del c
            self.mpl_new_ver = True
        except:
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
        """ A more generic version of voltage curve. """

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
