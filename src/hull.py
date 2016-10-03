# coding: utf-8
""" This file implements convex hull functionality
from database queries.
"""

from __future__ import print_function
from scipy.spatial import ConvexHull
from traceback import print_exc
from bson.son import SON
from bisect import bisect_left
from print_utils import print_failure, print_notify, print_warning
from chem_utils import get_capacities, get_molar_mass
import pymongo as pm
import re
import numpy as np
from mpldatacursor import datacursor
import matplotlib.pyplot as plt
import matplotlib.colors as colours


class QueryConvexHull():
    """ Implements a Convex Hull for formation energies
    from a fryan DBQuery.
    """
    def __init__(self, query, **kwargs):
        """ Accept query from fryan as argument. """
        self.query = query
        self.cursor = list(query.cursor)
        self.args = kwargs
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

        self.query.display_results(self.hull_cursor, hull=True)

        if self.args['subcmd'] == 'voltage':
            if self.args.get('debug'):
                self.set_plot_param()
                self.voltage_curve(self.stable_enthalpy_per_b, self.stable_comp, self.mu_enthalpy)
                self.metastable_voltage_profile()
            else:
                self.voltage_curve(self.stable_enthalpy_per_b, self.stable_comp, self.mu_enthalpy)
                self.set_plot_param()
                if self.args.get('subplot'):
                    self.subplot_voltage_hull()
                else:
                    self.plot_voltage_curve()
                self.plot_2d_hull()
        elif self.args.get('volume'):
            self.volume_curve(self.stable_comp, self.stable_vol)

        if self.args['subcmd'] == 'hull' and not self.args['no_plot']:
            if self.args.get('bokeh'):
                self.plot_2d_hull_bokeh()
            else:
                self.set_plot_param()
                self.plot_2d_hull()
            if self.args.get('volume'):
                self.plot_volume_curve()

    def get_chempots(self):
        """ Search for chemical potentials that match
        the structures in the query cursor.
        """
        query = self.query
        self.mu_enthalpy = np.zeros((2))
        self.mu_volume = np.zeros((2))
        query_dict = dict()
        if self.chem_pots is not None:
            self.fake_chempots()
        else:
            print(60*'─')
            self.match = len(self.elements)*[None]
            # scan for suitable chem pots in database
            for ind, elem in enumerate(self.elements):
                print('Scanning for suitable', elem, 'chemical potential...')
                query_dict['$and'] = list(query.calc_dict['$and'])
                query_dict['$and'].append(query.query_quality())
                query_dict['$and'].append(query.query_composition(custom_elem=[elem]))
                # if oqmd, only query composition, not parameters
                if query.args.get('tags') is not None:
                    query_dict['$and'].append(query.query_tags())
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

    def get_atoms_per_fu(self, doc):
        """ Calculate the number of atoms per formula unit. """
        atoms_per_fu = 0
        for j in range(len(doc['stoichiometry'])):
            atoms_per_fu += doc['stoichiometry'][j][1]
        return atoms_per_fu

    def get_formation_energy(self, doc):
        """ From given chemical potentials, calculate the simplest
        formation energy of the desired document.
        """
        formation = doc['enthalpy_per_atom']
        for mu in self.match:
            for j in range(len(doc['stoichiometry'])):
                if mu['stoichiometry'][0][0] == doc['stoichiometry'][j][0]:
                    formation -= (mu['enthalpy_per_atom'] * doc['stoichiometry'][j][1] /
                                  self.get_atoms_per_fu(doc))
        return formation

    def get_concentration(self, doc):
        """ Returns x for A_x B_{1-x}. """
        stoich = 0.0
        for elem in doc['stoichiometry']:
            if self.x_elem[0] in elem[0]:
                stoich = elem[1]/float(self.get_atoms_per_fu(doc))
        return stoich

    def get_array_from_cursor(self, cursor, key):
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
        return array

    def set_cursor_from_array(self, array, key):
        """ Updates the key-value pair for documents in
        internal cursor from a numpy array.
        """
        assert(len(array) == len(self.cursor))
        for ind, doc in enumerate(self.cursor):
            self.cursor[ind][key] = array[ind]
        return

    def get_hull_distances(self, structures):
        """ Returns array of hull distances. """
        hull_dist = np.zeros((len(structures)))
        hull_comp = self.structure_slice[self.hull.vertices, 0]
        hull_energy = self.structure_slice[self.hull.vertices, 1]
        hull_comp = np.asarray(hull_comp)
        hull_energy = hull_energy[np.argsort(hull_comp)]
        hull_comp = hull_comp[np.argsort(hull_comp)]
        for ind in range(len(structures)):
            # get the index of the next stoich on the hull from the current structure
            i = bisect_left(hull_comp, structures[ind, 0])
            energy_pair = (hull_energy[i-1], hull_energy[i])
            comp_pair = (hull_comp[i-1], hull_comp[i])
            # calculate equation of line between the two
            gradient = (energy_pair[1] - energy_pair[0]) / (comp_pair[1] - comp_pair[0])
            intercept = ((energy_pair[1] + energy_pair[0]) -
                         gradient * (comp_pair[1] + comp_pair[0])) / 2
            # calculate hull_dist
            hull_dist[ind] = structures[ind, 1] - (gradient * structures[ind, 0] + intercept)
        return hull_dist, hull_energy, hull_comp

    def hull_2d(self, dis=False):
        """ Create a convex hull for two elements. """
        query = self.query
        self.elements = query.args.get('composition')
        print(query.args.get('composition'))
        self.elements = [elem for elem in re.split(r'([A-Z][a-z]*)', self.elements[0]) if elem]
        if len(self.elements) > 3:
            print('Cannot create binary hull for more or less than 2 elements (yet!).')
            return
        self.get_chempots()
        print('Constructing hull...')
        # num_structures = len(self.cursor)
        # source_ind = np.zeros((num_structures))
        info = []
        # self.source_list = []
        # define hull by order in command-line arguments
        self.x_elem = [self.elements[0]]
        self.one_minus_x_elem = list(self.elements[1:])
        one_minus_x_elem = self.one_minus_x_elem
        # grab relevant information from query results; also make function?
        for ind, doc in enumerate(self.cursor):
            # num_fu = doc['num_fu']
            # calculate number of atoms of type B per formula unit
            nums_b = len(one_minus_x_elem)*[0]
            for elem in doc['stoichiometry']:
                for chem_pot_ind, chem_pot in enumerate(one_minus_x_elem):
                    if elem[0] == chem_pot:
                        nums_b[chem_pot_ind] += elem[1]
            # num_b = sum(nums_b)
            # get enthalpy and volume per unit B
            # self.cursor[ind]['enthalpy_per_b'] = doc['enthalpy'] / (num_b*num_fu)
            # self.cursor[ind]['cell_volume_per_b'] = doc['cell_volume'] / (num_b*num_fu)
            # source_dir = ''.join(doc['source'][0].split('/')[:-1])
            # self.cursor[ind]['source_idx'] = self.source_list.index(source_dir) + 1
            # if source_dir not in self.source_list:
                # self.source_list.append(source_dir)
            self.cursor[ind]['formation_enthalpy_per_atom'] = self.get_formation_energy(doc)
            self.cursor[ind]['concentration'] = self.get_concentration(doc)
        # put chem pots in same array as formation for easy hulls
        structures = np.vstack((self.get_array_from_cursor(self.cursor, 'concentration'),
                                self.get_array_from_cursor(self.cursor, 'formation_enthalpy_per_atom'))).T

        # create hull with SciPy routine, including only points with E_F < 0
        self.structure_slice = structures[np.where(structures[:, 1] <= 0 + 1e-9)]
        self.hull = ConvexHull(self.structure_slice)
        hull_dist, hull_energy, hull_comp = self.get_hull_distances(structures)
        self.set_cursor_from_array(hull_dist, 'hull_distance')

        # create hull_cursor to pass to other modules
        hull_cursor = [self.cursor[idx] for idx in np.where(hull_dist <= self.hull_cutoff + 1e-12)[0]]

        def get_text_info(self, cursor):
            # grab info for datacursor
            info = []
            docs = self.match[1:]
            ind = 0
            stoich_string = ''
            if len(docs) > 1:
                stoich_string += '['
            for ind, doc in enumerate(docs):
                for elem in doc['stoichiometry']:
                    stoich_string += elem[0]
                    stoich_string += '$_{' + str(elem[1]) + '}$' if elem[1] != 1 else ''
                if ind != len(docs) - 1:
                    stoich_string += ', '
            if len(docs) > 1:
                stoich_string += ']'
            info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string,
                                                                    doc['text_id'][0] + ' ' +
                                                                    doc['text_id'][1],
                                                                    doc['space_group'],
                                                                    hull_dist[ind]))
            for ind, doc in enumerate(self.cursor):
                stoich_string = ''
                for elem in doc['stoichiometry']:
                    stoich_string += elem[0]
                    stoich_string += '$_{' + str(elem[1]) + '}$' if elem[1] != 1 else ''
                info.append("{0:^10}\n{1:^24}\n{2:^5s}\n{3:2f} eV".format(stoich_string,
                                                                          doc['text_id'][0] + ' ' +
                                                                          doc['text_id'][1],
                                                                          doc['space_group'],
                                                                          hull_dist[ind+1]))
            doc = self.match[0]
            ind = len(hull_dist)-1
            stoich_string = ''
            for elem in doc['stoichiometry']:
                stoich_string += elem[0]
                stoich_string += '$_{' + str(elem[1]) + '}$' if elem[1] != 1 else ''
            info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string,
                                                                    doc['text_id'][0] + ' ' +
                                                                    doc['text_id'][1],
                                                                    doc['space_group'],
                                                                    hull_dist[ind]))

        self.structures = structures
        self.info = info
        # self.source_ind = source_ind
        self.hull_cursor = hull_cursor
        self.hull_dist = hull_dist
        self.hull_comp = hull_comp
        self.hull_energy = hull_energy

    def voltage_curve(self, stable_enthalpy_per_b, stable_comp, mu_enthalpy):
        """ Take convex hull and calculate voltages. """
        print('Generating voltage curve...')
        stable_num = []
        V = []
        x = []
        for i in range(len(stable_comp)):
            if 1-stable_comp[i] == 0:
                stable_num.append(1e5)
            else:
                stable_num.append(stable_comp[i]/(1-stable_comp[i]))
        # V.append(0)
        for i in range(len(stable_num)-1, 0, -1):
            V.append(-(stable_enthalpy_per_b[i] - stable_enthalpy_per_b[i-1]) /
                      (stable_num[i] - stable_num[i-1]) +
                      (mu_enthalpy[0]))
            x.append(stable_num[i])
        V.append(V[-1])
        x.append(0)
        self.voltages = np.asarray(V)
        self.x = np.asarray(x)
        # gravimetric capacity
        self.Q = get_capacities(x, get_molar_mass(self.elements[1]))
        return

    def volume_curve(self, stable_comp, stable_vol):
        """ Take stable compositions and volume and calculate
        volume expansion per "B" in AB binary.
        """
        stable_comp = stable_comp[:-1]
        stable_vol = stable_vol[:-1]
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
        if self.args.get('pdf') or self.args.get('png') or self.args.get('svg'):
            fig = plt.figure(facecolor=None, figsize=(5, 4))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
        scatter = []
        hull_scatter = []
        x_elem = [self.elements[0]]
        one_minus_x_elem = list(self.elements[1:])
        plt.draw()
        # star structures on hull
        for ind in range(len(self.hull.vertices)):
            if self.structure_slice[self.hull.vertices[ind], 1] <= 0:
                hull_scatter.append(ax.scatter(self.structure_slice[self.hull.vertices[ind], 0],
                                               self.structure_slice[self.hull.vertices[ind], 1],
                                               c=self.colours[1],
                                               marker='o', zorder=99999, edgecolor='k',
                                               s=self.scale*40, lw=1.5, alpha=1))
                                               # label=self.info[self.hull.vertices[ind]]))
                # ax.annotate(self.info[self.hull.vertices[ind]].split('\n')[0],
                            # xy=(self.structure_slice[self.hull.vertices[ind], 0],
                                # self.structure_slice[self.hull.vertices[ind], 1]),
                            # textcoords='data',
                            # ha='center',
                            # zorder=99999,
                            # xytext=(self.structure_slice[self.hull.vertices[ind], 0],
                                    # self.structure_slice[self.hull.vertices[ind], 1]-0.05))
        lw = self.scale * 0 if self.mpl_new_ver else 1
        # points for off hull structures
        if self.hull_cutoff == 0:
            # if no specified hull cutoff, ignore labels and colour
            # by distance from hull
            cmap_full = plt.cm.get_cmap('Dark2')
            cmap = colours.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_full.name, a=0, b=1),
                cmap_full(np.linspace(0.15, 0.4, 100)))
            scatter = ax.scatter(self.structures[:, 0], self.structures[:, 1],
                                 s=self.scale*40, lw=lw, alpha=0.9, c=self.hull_dist,
                                 edgecolor='k', zorder=300, cmap=cmap)
                                 # edgecolor='k', zorder=300, cmap=cmap)
            scatter = ax.scatter(self.structures[np.argsort(self.hull_dist), 0][::-1],
                                 self.structures[np.argsort(self.hull_dist), 1][::-1],
                                 s=self.scale*40, lw=lw, alpha=1, c=np.sort(self.hull_dist)[::-1],
                                 edgecolor='k', zorder=10000, cmap=cmap, norm=colours.LogNorm(0.02, 2))
            cbar = plt.colorbar(scatter, aspect=30, pad=0.02, ticks=[0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
            cbar.ax.set_yticklabels([0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
            cbar.set_label('Distance from hull (eV)')
        if self.hull_cutoff != 0:
            # if specified hull cutoff, label and colour those below
            for ind in range(len(self.structures)):
                if self.hull_dist[ind] <= self.hull_cutoff or self.hull_cutoff == 0:
                    c = self.colours[1]
                    scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1],
                                   s=self.scale*40, lw=lw, alpha=0.9, c=c, edgecolor='k'))
                                   # label=self.info[ind], zorder=300))
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
        if not dis and self.hull_cutoff != 0:
            datacursor(scatter[:], formatter='{label}'.format, draggable=False,
                       bbox=dict(fc='white'),
                       arrowprops=dict(arrowstyle='simple', alpha=1))
        ax.set_ylim(-0.1 if np.min(self.structure_slice[self.hull.vertices, 1]) > 0
                    else np.min(self.structure_slice[self.hull.vertices, 1])-0.15,
                    0.1 if np.max(self.structure_slice[self.hull.vertices, 1]) > 1
                    else np.max(self.structure_slice[self.hull.vertices, 1])+0.1)
        if len(one_minus_x_elem) == 1:
            ax.set_title(x_elem[0] + '$_\mathrm{x}$' + one_minus_x_elem[0] + '$_\mathrm{1-x}$')
        else:
            ax.set_title(x_elem[0] + '$_\mathrm{x}$(' + ''.join([elem for elem in one_minus_x_elem]) + ')$_\mathrm{1-x}$')
        plt.locator_params(nbins=3)
        ax.set_xlabel('$\mathrm{x}$')
        ax.grid(False)
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
        from bokeh.plotting import figure, show, output_file
        from bokeh.models import ColumnDataSource, HoverTool
        # x_elem = self.elements[0]
        # one_minus_x_elem = self.elements[1]
        # prepare data for bokeh
        tie_line_data = dict()
        # tie_line_data['info'] = list()
        tie_line_data['composition'] = list()
        tie_line_data['energy'] = list()
        for ind in range(len(self.hull.vertices)):
            if self.structure_slice[self.hull.vertices[ind], 1] <= 0:
                tie_line_data['composition'].append(self.structure_slice[self.hull.vertices[ind], 0])
                tie_line_data['energy'].append(self.structure_slice[self.hull.vertices[ind], 1])
                # tie_line_data['info'].append(self.info[self.hull.vertices[ind]])
        tie_line_data['energy'] = np.asarray(tie_line_data['energy'])
        tie_line_data['composition'] = np.asarray(tie_line_data['composition'])
        tie_line_data['energy'] = tie_line_data['energy'][np.argsort(tie_line_data['composition'])]
        # tie_line_data['info'] = [tie_line_data['info'][i]
                                 # for i in list(np.argsort(tie_line_data['composition']))]
        tie_line_data['composition'] = np.sort(tie_line_data['composition'])
        hull_data = dict()
        # hull_data['info'] = list()
        # points for off hull structures
        hull_data['composition'] = self.structures[:, 0]
        hull_data['energy'] = self.structures[:, 1]
        hull_data['hull_dist'] = self.hull_dist
        # hull_data['info'] = self.info
        cmap_limits = [0, 0.5]
        colormap = plt.cm.get_cmap('Dark2')
        cmap_input = np.interp(hull_data['hull_dist'], cmap_limits, [0.15, 0.4], left=0.15, right=0.4)
        colours = colormap(cmap_input, 1, True)
        bokeh_colours = ["#%02x%02x%02x" % (r, g, b) for r, g, b in colours[:, 0:3]]
        fixed_colours = colormap([0.0, 0.15], 1, True)
        tie_line_colour, on_hull_colour = ["#%02x%02x%02x" % (r, g, b) for r, g, b in fixed_colours[:, 0:3]]

        tie_line_source = ColumnDataSource(data=tie_line_data)
        hull_source = ColumnDataSource(data=hull_data)

        # hover = HoverTool(tooltips="""
                          # <div>
                              # <div>
                                  # <span style="font-size: 12px;">@info</span>
                              # </div>
                          # </div>
                          # """)

        tools = ['pan', 'wheel_zoom']
        # tools.append(hover)
        fig = figure(tools=tools)

        fig.xaxis.axis_label = 'x'
        fig.yaxis.axis_label = 'Formation energy (eV/atom)'

        fig.xaxis.axis_label_text_font_size = '20pt'
        fig.yaxis.axis_label_text_font_size = '20pt'

        fig.line('composition', 'energy',
                 source=tie_line_source,
                 line_width=4,
                 line_color=tie_line_colour)
        fig.scatter('composition', 'energy',
                    source=hull_source,
                    alpha=1,
                    size=10,
                    fill_color=bokeh_colours,
                    line_color=None)
        fig.square('composition', 'energy',
                   source=tie_line_source,
                   line_color='black',
                   color=on_hull_colour,
                   line_width=2,
                   alpha=1,
                   size=10)

        fig.plot_width = 800
        fig.plot_height = 800
        output_file('test.html', title='test.py example')
        show(fig)

    def plot_voltage_curve(self):
        """ Plot calculated voltage curve. """
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
        for i in range(2, len(self.voltages)):
            axQ.plot([self.Q[i], self.Q[i]], [self.voltages[i], self.voltages[i-1]],
                     lw=2, c=self.colours[0])
            axQ.plot([self.Q[i-1], self.Q[i]], [self.voltages[i-1], self.voltages[i-1]],
                     lw=2, c=self.colours[0])
            if abs(np.ceil(self.x[i-1])-self.x[i-1]) > 1e-8:
                string_stoich = str(round(self.x[i-1], 1))
            else:
                string_stoich = str(int(np.ceil(self.x[i-1])))
                if string_stoich is '1':
                    string_stoich = ''
        axQ.set_ylabel('Voltage (V)')
        axQ.set_xlabel('Gravimetric cap. (mAh/g)')
        start, end = axQ.get_ylim()
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
        if self.args.get('pdf') or self.args.get('png'):
            fig = plt.figure(facecolor=None, figsize=(2.7, 2.7))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
        hull_vols = []
        hull_comps = []
        for i in range(len(self.vol_per_y)):
            if self.stable_hull_dist[i] <= 0:
                hull_vols.append(self.vol_per_y[i])
                hull_comps.append(self.x[i])
                ax.scatter(self.x[i], self.vol_per_y[i]/self.vol_per_y[0], marker='o', s=40, edgecolor='k', lw=1.5,
                           c=self.colours[1], zorder=1000)
            else:
                ax.scatter(self.x[i], self.vol_per_y[i]/self.vol_per_y[0], marker='o', s=30, edgecolor='k', lw=0,
                           c=self.colours[1], zorder=900, alpha=1)
        hull_comps, hull_vols = np.asarray(hull_comps), np.asarray(hull_vols)
        ax.plot(hull_comps, hull_vols/self.vol_per_y[0], marker='o', lw=4,
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
        if not dis:
            datacursor(scatter[:], formatter='{label}'.format, draggable=False,
                       bbox=dict(fc='white'),
                       arrowprops=dict(arrowstyle='simple', alpha=1))
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
        from palettable.colorbrewer.qualitative import Dark2_8
        from palettable.colorbrewer.qualitative import Set3_10
        # if len(self.source_list) < 6:
            # self.colours = Dark2_8.hex_colors[1:len(self.source_list)+1]
        # else:
        self.colours = Dark2_8.hex_colors[1:]
        self.colours.extend(Dark2_8.hex_colors[1:])
        # first colour reserved for hull
        self.colours.insert(0, Dark2_8.hex_colors[0])
        # penultimate colour reserved for off hull above cutoff
        self.colours.append(Dark2_8.hex_colors[-1])
        # last colour reserved for OQMD
        self.colours.append(Set3_10.hex_colors[-1])
        return

    def metastable_voltage_profile(self):
        """ Construct a smeared voltage profile from metastable hull,
        weighting either with Boltzmann, PDF overlap of nothing.
        """
        print('Generating metastable voltage profile...')
        structures = self.hull_cursor[:-1]
        guest_atoms = []
        for i in range(len(structures)):
            stoich = self.get_formation_energy(structures[i])
            try:
                structures['num_a'] = stoich / (1 - stoich)
            except:
                structures['num_a'] = float('inf')
        for i in range(len(structures)):
            guest_atoms.append(structures[i]['num_a'])
            structures[i]['capacity'] = get_capacities(guest_atoms[-1], get_molar_mass(self.elements[1]))
        # guest_atoms = np.asarray(guest_atoms)
        # max_guest = np.max(guest_atoms)
        # grab reference cathode
        reference_cathode = self.hull_cursor[-1]

        # set up running average voltage profiles
        num_divisions = 100
        running_average_profile = np.zeros((num_divisions))
        advanced_average_profile = np.zeros_like(running_average_profile)
        diff_average_profile = np.ones_like(running_average_profile)
        # replace with max cap
        capacity_space = np.linspace(0, 2000, num_divisions)
        all_voltage_profiles = []
        # set convergance of voltage profile tolerance
        tolerance = 1e-9 * num_divisions

        def boltzmann(structures):
            idx = np.random.randint(0, len(structures))
            return structures[idx]

        def pdf_overlap(structures):
            idx = np.random.randint(0, len(structures))
            return structures[idx]

        def pick_metastable_structure(structures, weighting=None, last=None):
            """ Pick an acceptable metastable structure. """
            count = 0
            while True:
                count += 1
                if count > 1000:
                    exit('Something went wrong...')
                test = structures[np.random.randint(0, len(structures))]
                if test['num_a'] < last['num_a']:
                    print(count)
                    return test
            else:
                return weighting(structures)

        def get_voltage_profile_segment(structure_new, structure_old):
            """ Return voltage between two structures. """
            V = (-(structure_new['enthalpy_per_b'] - structure_old['enthalpy_per_b']) /
                 (structure_new['num_a'] - structure_old['num_a']) +
                 (self.mu_enthalpy[0]))
            if structure_old['num_a'] == float('inf'):
                V = 0
            return V

        def interp_steps(xspace, x, y):
            """ Interpolate (x,y) across xspace with step functions. """
            yspace = np.zeros_like(xspace)
            for x_val, y_val in zip(x, y):
                yspace[np.where(xspace <= x_val)] = y_val
            return yspace

        weighting = None
        converged = False
        convergence_window = 3
        last_converged = False
        num_contribs = 0
        num_converged = 0
        num_trials = 100000
        while not converged:
            if np.abs(diff_average_profile).sum() < tolerance:
                if num_converged > convergence_window:
                    converged = True
                    break
                if not last_converged:
                    last_converged = True
                    num_converged = 1
                if last_converged:
                    num_converged += 1
            if num_contribs > num_trials:
                break
            # print('Sampling path', num_contribs)
            delithiated = False
            voltage_profile = []
            prev = reference_cathode
            path_structures = list(structures)
            while not delithiated:
                next = pick_metastable_structure(path_structures, weighting=weighting, last=prev)
                idx = path_structures.index(next)
                path_structures = path_structures[:idx]
                if prev['num_a'] != float('inf'):
                    voltage_profile.append([next['capacity'], get_voltage_profile_segment(next, prev)])
                if next['num_a'] == 0:
                    delithiated = True
                else:
                    prev = next
            if delithiated and len(voltage_profile) != 0:
                # piecewise linear interpolation of V and add to running average
                voltage_profile = np.asarray(voltage_profile)
                # print(voltage_profile)
                advanced_average_profile += interp_steps(capacity_space, voltage_profile[:, 0], voltage_profile[:, 1])
                num_contribs += 1
                diff_average_profile = advanced_average_profile/num_contribs - running_average_profile
                # print(advanced_average_profile)
                all_voltage_profiles.append(voltage_profile)
                running_average_profile = advanced_average_profile/num_contribs

        print('Convergence achieved with', num_contribs, 'paths.')
        print(np.abs(diff_average_profile).sum(), tolerance)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # for idx, vp in enumerate(all_voltage_profiles):
            # if idx % len(all_voltage_profiles) == 0:
                # ax.scatter(vp[:, 0], vp[:, 1], s=50)
            # print(vp)
        for i in range(2, len(self.voltages)):
            ax.plot([self.Q[i], self.Q[i]], [self.voltages[i], self.voltages[i-1]],
                    lw=2, c=self.colours[0])
            ax.plot([self.Q[i-1], self.Q[i]], [self.voltages[i-1], self.voltages[i-1]],
                    lw=2, c=self.colours[0])
        ax.scatter(capacity_space, running_average_profile, c='r', marker='*', s=50)
        plt.show()
