# coding: utf-8
""" This file implements convex hull functionality
from database queries.
"""
from __future__ import print_function
from scipy.spatial import ConvexHull
from bson.son import SON
from bisect import bisect_left
from print_utils import print_failure
import pymongo as pm
import re
import numpy as np
from mpldatacursor import datacursor
import matplotlib.pyplot as plt
import matplotlib.colors as colours


class QueryConvexHull():
    """
    Implements a Convex Hull for formation energies
    from a fryan DBQuery.
    """
    def __init__(self, query, *args):
        """ Accept query from fryan as argument. """
        self.query = query
        self.cursor = list(query.cursor)
        self.args = args[0]
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
        self.binary_hull()
        if self.args['subcmd'] == 'voltage':
            print('Generating voltage curve...')
            self.voltage_curve(self.stable_enthalpy_per_b, self.stable_comp, self.mu_enthalpy)
            self.set_plot_param()
            if self.args.get('subplot'):
                self.subplot_voltage_hull()
            else:
                self.plot_voltage_curve()
            self.set_plot_param()
            self.plot_hull()
        elif self.args.get('volume'):
            self.volume_curve(self.stable_comp, self.stable_vol)
        if self.args['subcmd'] == 'hull' and not self.args['no_plot']:
            if self.args.get('bokeh'):
                self.plot_hull_bokeh()
            else:
                self.set_plot_param()
                self.plot_hull()
            if self.args.get('volume'):
                self.plot_volume_curve()

    def get_chempots(self):
        """ Search for chemical potentials that match
        the structures in the query cursor.
        """
        query = self.query
        self.mu_enthalpy = np.zeros((2))
        self.mu_volume = np.zeros((2))
        self.match = [None, None]
        query_dict = dict()
        print(60*'─')
        if self.chem_pots is not None:
            # read chem pots from command line
            self.mu_enthalpy[0] = self.chem_pots[0]
            self.match[0] = dict()
            self.match[0]['enthalpy_per_atom'] = self.mu_enthalpy[0]
            self.match[0]['text_id'] = ['command', 'line']
            self.match[0]['stoichiometry'] = [[self.elements[0], 1]]
            self.match[0]['space_group'] = 'xxx'
            self.mu_enthalpy[1] = self.chem_pots[1]
            self.match[1] = dict()
            self.match[1]['enthalpy_per_atom'] = self.mu_enthalpy[1]
            self.match[1]['text_id'] = ['command', 'line']
            self.match[1]['text_id'] = ['command', 'line']
            self.match[1]['stoichiometry'] = [[self.elements[1], 1]]
            self.match[1]['space_group'] = 'xxx'
            print('Using custom energies of', self.mu_enthalpy[0], 'eV/atom '
                  'and', self.mu_enthalpy[1], 'eV/atom as chemical potentials.')
            print(60*'─')
        else:
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
                for doc_ind, doc in enumerate(mu_cursor):
                    if doc_ind == 0:
                        self.match[ind] = doc
                        break
                if self.match[ind] is not None:
                    self.mu_enthalpy[ind] = float(self.match[ind]['enthalpy_per_atom'])
                    self.mu_volume[ind] = float(self.match[ind]['cell_volume'] /
                                                self.match[ind]['num_atoms'])
                    print('Using', ''.join([self.match[ind]['text_id'][0], ' ',
                          self.match[ind]['text_id'][1]]), 'as chem pot for', elem)
                    print(60*'─')
                else:
                    print_failure('No possible chem pots found for ' + elem + '.')
                    exit()
        return

    def binary_hull(self, dis=False):
        """ Create a convex hull for two elements. """
        query = self.query
        self.include_oqmd = query.args.get('self.include_oqmd')
        self.elements = query.args.get('composition')
        self.elements = [elem for elem in re.split(r'([A-Z][a-z]*)', self.elements[0]) if elem]
        if len(self.elements) != 2:
            print('Cannot create binary hull for more or less than 2 elements (yet!).')
            return
        self.get_chempots()
        print('Constructing hull...')
        num_structures = len(self.cursor)
        formation = np.zeros((num_structures))
        stoich = np.zeros((num_structures))
        enthalpy = np.zeros((num_structures))
        volume = np.zeros((num_structures))
        disorder = np.zeros((num_structures))
        source_ind = np.zeros((num_structures+2), dtype=int)
        hull_dist = np.zeros((num_structures+2))
        info = []
        self.source_list = []
        if dis:
            from disorder import disorder_hull
        # define hull by order in command-line arguments
        x_elem = self.elements[0]
        one_minus_x_elem = self.elements[1]
        # grab relevant information from query results; also make function?
        for ind, doc in enumerate(self.cursor):
            atoms_per_fu = doc['stoichiometry'][0][1] + doc['stoichiometry'][1][1]
            num_fu = doc['num_fu']
            # calculate number of atoms of type B per formula unit
            if doc['stoichiometry'][0][0] == one_minus_x_elem:
                num_b = doc['stoichiometry'][0][1]
            elif doc['stoichiometry'][1][0] == one_minus_x_elem:
                num_b = doc['stoichiometry'][1][1]
            else:
                print_failure('Something went wrong!')
                exit()
            # get enthalpy and volume per unit B
            enthalpy[ind] = doc['enthalpy'] / (num_b*num_fu)
            formation[ind] = doc['enthalpy_per_atom']
            volume[ind] = doc['cell_volume'] / (num_b*num_fu)
            source_dir = ''.join(doc['source'][0].split('/')[:-1])
            if source_dir in self.source_list:
                source_ind[ind] = self.source_list.index(source_dir) + 1
            else:
                self.source_list.append(source_dir)
                source_ind[ind] = self.source_list.index(source_dir) + 1
            for mu in self.match:
                for j in range(len(doc['stoichiometry'])):
                    if mu['stoichiometry'][0][0] == doc['stoichiometry'][j][0]:
                        formation[ind] -= (mu['enthalpy_per_atom'] * doc['stoichiometry'][j][1] /
                                           atoms_per_fu)
            for elem in doc['stoichiometry']:
                if x_elem in elem[0]:
                    stoich[ind] = elem[1]/float(atoms_per_fu)
            if dis:
                disorder[ind], warren = disorder_hull(doc)
        # put chem pots in same array as formation for easy hulls
        formation = np.append([0.0], formation)
        formation = np.append(formation, [0.0])
        enthalpy = np.append(self.mu_enthalpy[1], enthalpy)
        enthalpy = np.append(enthalpy, self.mu_enthalpy[0])
        volume = np.append(self.mu_volume[1], volume)
        volume = np.append(volume, self.mu_volume[0])
        ind = len(formation)-3
        stoich = np.append([0.0], stoich)
        stoich = np.append(stoich, [1.0])
        structures = np.vstack((stoich, formation)).T

        # create hull with SciPy routine
        self.hull = ConvexHull(structures)

        hull_energy = []
        hull_comp = []
        hull_enthalpy = []
        hull_volume = []
        hull_cursor = []
        for ind in range(len(self.hull.vertices)):
            if structures[self.hull.vertices[ind], 1] <= 0:
                hull_energy.append(structures[self.hull.vertices[ind], 1])
                hull_enthalpy.append(enthalpy[self.hull.vertices[ind]])
                hull_comp.append(structures[self.hull.vertices[ind], 0])
                hull_volume.append(volume[ind])
        # calculate distance to hull of all structures
        hull_energy = np.asarray(hull_energy)
        hull_enthalpy = np.asarray(hull_enthalpy)
        hull_comp = np.asarray(hull_comp)
        hull_volume = np.asarray(hull_volume)
        hull_energy = hull_energy[np.argsort(hull_comp)]
        hull_enthalpy = hull_enthalpy[np.argsort(hull_comp)]
        hull_volume = hull_volume[np.argsort(hull_comp)]
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
        # if below cutoff, include in arg to voltage curve
        stable_energy = list(hull_energy)
        stable_enthalpy_per_b = list(hull_enthalpy)
        stable_comp = list(hull_comp)
        stable_vol = list(hull_volume)
        for ind in range(len(structures)):
            if hull_dist[ind] <= self.hull_cutoff:
                # recolour if under cutoff
                source_ind[ind] = 0
                # get lowest enthalpy at particular comp
                if structures[ind, 0] not in stable_comp:
                    stable_energy.append(structures[ind, 1])
                    stable_enthalpy_per_b.append(enthalpy[ind])
                    stable_comp.append(structures[ind, 0])
                    stable_vol.append(volume[ind])
        # create hull_cursor to pass to other modules
        # skip last and first as they are chem pots
        self.match[0]['enthalpy_per_atom'] = 0.0
        self.match[1]['enthalpy_per_atom'] = 0.0
        hull_cursor.append(self.match[0])
        for ind in range(1, len(hull_dist)-1):
            if hull_dist[ind] <= self.hull_cutoff+1e-12:
                self.cursor[ind-1]['enthalpy_per_atom'] = hull_dist[ind]
                # take ind-1 to ignore first chem pot
                hull_cursor.append(self.cursor[ind-1])
        hull_cursor.append(self.match[1])
        # grab info for datacursor
        info = []
        doc = self.match[0]
        ind = 0
        stoich_string = str(doc['stoichiometry'][0][0])
        info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string,
                                                                doc['text_id'][0] + ' ' +
                                                                doc['text_id'][1],
                                                                doc['space_group'],
                                                                hull_dist[ind]))
        for ind, doc in enumerate(self.cursor):
            stoich_string = str(doc['stoichiometry'][0][0])
            stoich_string += '$_{' + str(doc['stoichiometry'][0][1]) + '}$' \
                if doc['stoichiometry'][0][1] != 1 else ''
            stoich_string += str(doc['stoichiometry'][1][0])
            stoich_string += '$_{' + str(doc['stoichiometry'][1][1]) + '}$' \
                if doc['stoichiometry'][1][1] != 1 else ''
            info.append("{0:^10}\n{1:^24}\n{2:^5s}\n{3:2f} eV".format(stoich_string,
                                                                      doc['text_id'][0] + ' ' +
                                                                      doc['text_id'][1],
                                                                      doc['space_group'],
                                                                      hull_dist[ind+1]))
        doc = self.match[1]
        ind = len(hull_dist)-1
        stoich_string = str(doc['stoichiometry'][0][0])
        info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string,
                                                                doc['text_id'][0] + ' ' +
                                                                doc['text_id'][1],
                                                                doc['space_group'],
                                                                hull_dist[ind]))

        stable_energy = np.asarray(stable_energy)
        stable_comp = np.asarray(stable_comp)
        stable_enthalpy_per_b = np.asarray(stable_enthalpy_per_b)
        stable_vol = np.asarray(stable_vol)
        stable_energy = stable_energy[np.argsort(stable_comp)]
        stable_enthalpy_per_b = stable_enthalpy_per_b[np.argsort(stable_comp)]
        stable_vol = stable_vol[np.argsort(stable_comp)]
        stable_comp = stable_comp[np.argsort(stable_comp)]

        self.structures = structures
        self.info = info
        self.source_ind = source_ind
        self.source_ind = source_ind
        self.hull_cursor = hull_cursor
        self.hull_dist = hull_dist
        self.hull_comp = hull_comp
        self.hull_energy = hull_energy
        self.stable_enthalpy_per_b = stable_enthalpy_per_b
        self.stable_vol = stable_vol
        self.stable_comp = stable_comp

    def voltage_curve(self, stable_enthalpy_per_b, stable_comp, mu_enthalpy):
        """ Take convex hull and calculate voltages. """
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
        self.voltages = V
        self.x = x
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
            if 1-stable_comp[i] == 0:
                x.append(1e5)
            else:
                x.append(stable_comp[i]/(1-stable_comp[i]))
            v.append(stable_vol[i])
        self.x = x
        self.vol_per_y = v
        return

    def plot_hull(self, dis=False):
        """ Plot calculated hull. """
        if self.args.get('pdf'):
            fig = plt.figure(facecolor=None, figsize=(7, 4))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
        scatter = []
        hull_scatter = []
        x_elem = self.elements[0]
        one_minus_x_elem = self.elements[1]
        plt.draw()
        # star structures on hull
        for ind in range(len(self.hull.vertices)):
            if self.structures[self.hull.vertices[ind], 1] <= 0:
                hull_scatter.append(ax.scatter(self.structures[self.hull.vertices[ind], 0],
                                               self.structures[self.hull.vertices[ind], 1],
                                               c=self.colours[0],
                                               marker='*', zorder=99999, edgecolor='k',
                                               s=self.scale*150, lw=1, alpha=1,
                                               label=self.info[self.hull.vertices[ind]]))
                ax.annotate(self.info[self.hull.vertices[ind]].split('\n')[0],
                            xy=(self.structures[self.hull.vertices[ind], 0],
                                self.structures[self.hull.vertices[ind], 1]),
                            textcoords='data',
                            ha='center',
                            xytext=(self.structures[self.hull.vertices[ind], 0],
                                    self.structures[self.hull.vertices[ind], 1]-0.1))
        lw = self.scale * 0 if self.mpl_new_ver else 1
        # points for off hull structures
        if self.hull_cutoff == 0:
            # if no specified hull cutoff, ignore labels and colour
            # by distance from hull
            cmap_full = plt.cm.get_cmap('GnBu_r')
            cmap = colours.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap_full.name, a=0, b=0.9),
                cmap_full(np.linspace(0, 0.9, 100)))
            scatter = ax.scatter(self.structures[:, 0], self.structures[:, 1],
                                 s=self.scale*30, lw=lw, alpha=0.9, c=self.hull_dist,
                                 edgecolor='k', zorder=300, cmap=cmap)
            cbar = plt.colorbar(scatter, aspect=30, pad=0.02)
            cbar.set_label('Distance from hull (eV)')
        if self.hull_cutoff != 0:
            # if specified hull cutoff, label and colour those below
            for ind in range(len(self.structures)):
                if self.hull_dist[ind] <= self.hull_cutoff or self.hull_cutoff == 0:
                    c = self.colours[1]
                    scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1],
                                   s=self.scale*30, lw=lw, alpha=0.9, c=c, edgecolor='k',
                                   label=self.info[ind], zorder=300))
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
        ax.set_ylim(-0.1 if np.min(self.structures[self.hull.vertices, 1]) > 0
                    else np.min(self.structures[self.hull.vertices, 1])-0.15,
                    0.5 if np.max(self.structures[self.hull.vertices, 1]) > 1
                    else np.max(self.structures[self.hull.vertices, 1])+0.1)
        ax.set_title(x_elem+'$_\mathrm{x}$'+one_minus_x_elem+'$_\mathrm{1-x}$')
        plt.locator_params(nbins=3)
        ax.set_xlabel('$\mathrm{x}$')
        ax.grid(False)
        ax.set_xticks([0, 0.33, 0.5, 0.66, 1])
        ax.set_ylabel('E$_\mathrm{F}$ (eV/atom)')
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull.pdf',
                        dpi=200, bbox_inches='tight')
        else:
            plt.show()

    def plot_hull_bokeh(self):
        """ Plot interactive hull with Bokeh. """
        from bokeh.plotting import figure, show, output_file
        from bokeh.models import ColumnDataSource, HoverTool
        # x_elem = self.elements[0]
        # one_minus_x_elem = self.elements[1]
        # prepare data for bokeh
        tie_line_data = dict()
        tie_line_data['info'] = list()
        tie_line_data['composition'] = list()
        tie_line_data['energy'] = list()
        for ind in range(len(self.hull.vertices)):
            if self.structures[self.hull.vertices[ind], 1] <= 0:
                tie_line_data['composition'].append(self.structures[self.hull.vertices[ind], 0])
                tie_line_data['energy'].append(self.structures[self.hull.vertices[ind], 1])
                tie_line_data['info'].append(self.info[self.hull.vertices[ind]])
        tie_line_data['energy'] = np.asarray(tie_line_data['energy'])
        tie_line_data['composition'] = np.asarray(tie_line_data['composition'])
        tie_line_data['energy'] = tie_line_data['energy'][np.argsort(tie_line_data['composition'])]
        tie_line_data['info'] = [tie_line_data['info'][i]
                                 for i in list(np.argsort(tie_line_data['composition']))]
        tie_line_data['composition'] = np.sort(tie_line_data['composition'])
        tie_line_data['colour'] = len(tie_line_data['energy']) * ['blue']
        hull_data = dict()
        hull_data['energy'] = list()
        hull_data['composition'] = list()
        hull_data['info'] = list()
        # points for off hull structures
        hull_data['composition'] = self.structures[:, 0]
        hull_data['energy'] = self.structures[:, 1]
        hull_data['info'] = self.info
        hull_data['colour'] = len(hull_data['energy']) * ['red']

        tie_line_source = ColumnDataSource(data=tie_line_data)
        hull_source = ColumnDataSource(data=hull_data)

        hover = HoverTool(tooltips="""
                          <div>
                              <div>
                                  <span style="font-size: 12px;">@info</span>
                              </div>
                          </div>
                          """)

        tools = ['pan', 'wheel_zoom']
        tools.append(hover)
        fig = figure(tools=tools)

        fig.line('composition', 'energy',
                 source=tie_line_source,
                 line_color='blue')
        fig.circle('composition', 'energy',
                   source=hull_source,
                   alpha=1,
                   size=5,
                   color='colour')
        fig.square('composition', 'energy',
                   source=tie_line_source,
                   color='colour',
                   alpha=1,
                   size=10)

        fig.plot_width = 800
        fig.plot_height = 800
        output_file('test.html', title='test.py example')
        show(fig)

    def plot_voltage_curve(self):
        """ Plot calculated voltage curve. """
        if self.args.get('pdf'):
            fig = plt.figure(facecolor=None, figsize=(3, 1.5))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
        for i in range(2, len(self.voltages)):
            # ax.scatter(self.x[i-1], self.voltages[i-1],
                       # marker='*', s=100, edgecolor='k', c=self.colours[0], zorder=1000)
            ax.plot([self.x[i], self.x[i]], [self.voltages[i], self.voltages[i-1]],
                    lw=2, c=self.colours[0])
            ax.plot([self.x[i-1], self.x[i]], [self.voltages[i-1], self.voltages[i-1]],
                    lw=2, c=self.colours[0])
        ax.set_ylabel('Voltage (V)')
        plt.locator_params(nbins=4)
        ax.set_xlim(0, np.max(np.asarray(self.x[1:]))+1)
        ax.set_ylim(np.min(np.asarray(self.voltages[2:]))-0.1,
                    np.max(np.asarray(self.voltages[2:]))+0.1)
        ax.set_title('$\mathrm{'+self.elements[0]+'_x'+self.elements[1]+'}$')
        ax.set_xlabel('$x$')
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_voltage.pdf',
                        dpi=300, bbox_inches='tight')
        else:
            plt.show()

    def plot_volume_curve(self):
        """ Plot calculate volume curve. """
        if self.args.get('pdf'):
            fig = plt.figure(facecolor=None, figsize=(3, 1.5))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
        ax.scatter(self.x, self.vol_per_y, marker='*', s=100, edgecolor='k',
                   c=self.colours[0], zorder=1000)
        ax2 = ax.twinx()
        ax2.scatter(self.x, self.vol_per_y/self.vol_per_y[0], alpha=0)
        ax.plot(self.x, self.vol_per_y, marker='*', lw=2,
                c=self.colours[0], zorder=100)
        ax.set_ylabel('Volume per ' + self.elements[1] + ' atom ($\AA^3$)')
        ax2.set_ylabel('Volume ratio with to bulk ' + self.elements[1])
        ax.set_title('$\mathrm{'+self.elements[0]+'_x'+self.elements[1]+'}$')
        ax.set_xlabel('$x$')
        ax.set_ylim(0.9*self.vol_per_y[0])
        ax2.set_ylim(0.9)
        ax2.grid('off')
        if self.args.get('pdf'):
            plt.savefig(self.elements[0]+self.elements[1]+'_voltage.pdf',
                        dpi=300, bbox_inches='tight')
        else:
            plt.show()

    def subplot_voltage_hull(self, dis=False):
        """ Plot calculated hull with inset voltage curve. """
        if self.args.get('pdf'):
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
            if self.structures[self.hull.vertices[ind], 1] <= 0:
                hull_scatter.append(ax.scatter(self.structures[self.hull.vertices[ind], 0],
                                               self.structures[self.hull.vertices[ind], 1],
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
        ax.set_ylim(-0.1 if np.min(self.structures[self.hull.vertices, 1]) > 0
                    else np.min(self.structures[self.hull.vertices, 1])-0.1,
                    0.5 if np.max(self.structures[self.hull.vertices, 1]) > 1
                    else np.max(self.structures[self.hull.vertices, 1])+0.1)
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
        if self.args.get('pdf'):
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
        if len(self.source_list) < 6:
            self.colours = Dark2_8.hex_colors[1:len(self.source_list)+1]
        else:
            self.colours = Dark2_8.hex_colors[1:]
            self.colours.extend(Dark2_8.hex_colors[1:])
        # first colour reserved for hull
        self.colours.insert(0, Dark2_8.hex_colors[0])
        # penultimate colour reserved for off hull above cutoff
        self.colours.append(Dark2_8.hex_colors[-1])
        # last colour reserved for OQMD
        self.colours.append(Set3_10.hex_colors[-1])
        return
