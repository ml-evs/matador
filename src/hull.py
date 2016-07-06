#!/usr/bin/python
# coding: utf-8
""" This file implements convex hull functionality
from database queries.
"""
from __future__ import print_function
from scipy.spatial import ConvexHull
from bson.son import SON
from bisect import bisect_left
from print_utils import print_failure, print_success, print_warning
import pymongo as pm
import re
import numpy as np
from mpldatacursor import datacursor
import matplotlib.pyplot as plt


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
        if self.args.get('hull_cutoff') is not None:
            self.hull_cutoff = float(self.args['hull_cutoff'])
        else:
            self.hull_cutoff = 0.0
        if self.args.get('self.include_oqmd'):
            self.include_oqmd = True
        self.binary_hull()
        if self.args['subcmd'] == 'voltage':
            print('Generating voltage curve...')
            self.voltage_curve(self.stable_enthalpy, self.stable_comp, self.mu_enthalpy)
            self.set_plot_param()
            if self.args.get('subplot'):
                self.subplot_voltage_hull()
            else:
                self.plot_voltage_curve()
        elif self.args['subcmd'] == 'hull':
            self.set_plot_param()
            self.plot_hull()

    def binary_hull(self, dis=False):
        """ Create a convex hull for two elements. """
        query = self.query
        self.include_oqmd = query.args.get('self.include_oqmd')
        self.elements = query.args.get('composition')
        self.elements = [elem for elem in re.split(r'([A-Z][a-z]*)', self.elements[0]) if elem]
        if len(self.elements) != 2:
            print('Cannot create binary hull for more or less than 2 elements (yet!).')
            return
        # try to get decent chemical potentials:
        # this relies on all of the first composition query
        # having the same parameters; need to think about this
        # should probably refactor this into a function
        mu_enthalpy = np.zeros((2))
        match = [None, None]
        query_dict = dict()
        print(60*'─')
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
                    match[ind] = doc
                    break
            if match[ind] is not None:
                mu_enthalpy[ind] = float(match[ind]['enthalpy_per_atom'])
                print('Using', ''.join([match[ind]['text_id'][0], ' ',
                      match[ind]['text_id'][1]]), 'as chem pot for', elem)
                print(60*'─')
            else:
                print_failure('No possible chem pots found for ' + elem + '.')
                exit()
        # include OQMD structures if desired, first find chem pots
        if self.include_oqmd:
            oqmd_mu_enthalpy = np.zeros((2))
            oqmd_match = [None, None]
            oqmd_query_dict = dict()
            for ind, elem in enumerate(self.elements):
                print('Scanning for suitable', elem, 'OQMD chemical potential...')
                oqmd_query_dict = query.query_composition(custom_elem=[elem])
                oqmd_mu_cursor = query.oqmd_repo.find(SON(oqmd_query_dict))
                oqmd_mu_cursor.sort('enthalpy_per_atom', pm.ASCENDING)
                for doc_ind, doc in enumerate(oqmd_mu_cursor):
                    if doc_ind == 0:
                        oqmd_match[ind] = doc
                        break
                if oqmd_match[ind] is not None:
                    oqmd_mu_enthalpy[ind] = float(oqmd_match[ind]['enthalpy_per_atom'])
                    print_success('Using ' + ''.join([oqmd_match[ind]['text_id'][0] + ' ' +
                          match[ind]['text_id'][1]]) + ' as OQMD chem pot for ' + elem)
                    print(60*'─')
                else:
                    print_failure('No possible chem pots found for ' + elem + '.')
                    exit()
        print('Constructing hull...')
        num_structures = len(self.cursor)
        formation = np.zeros((num_structures))
        stoich = np.zeros((num_structures))
        enthalpy = np.zeros((num_structures))
        disorder = np.zeros((num_structures))
        source_ind = np.zeros((num_structures+2), dtype=int)
        hull_dist = np.zeros((num_structures+2))
        info = []
        self.source_list = []
        if self.include_oqmd:
            oqmd_num_structures = query.oqmd_cursor.count()
            oqmd_formation = np.zeros((oqmd_num_structures))
            oqmd_stoich = np.zeros((oqmd_num_structures))
            oqmd_enthalpy = np.zeros((oqmd_num_structures))
            oqmd_info = []
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
            # get enthalpy per unit B
            enthalpy[ind] = doc['enthalpy'] / (num_b*num_fu)
            formation[ind] = doc['enthalpy_per_atom']
            source_dir = ''.join(doc['source'][0].split('/')[:-1])
            if source_dir in self.source_list:
                source_ind[ind] = self.source_list.index(source_dir) + 1
            else:
                self.source_list.append(source_dir)
                source_ind[ind] = self.source_list.index(source_dir) + 1
            for mu in match:
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
        enthalpy = np.append(mu_enthalpy[1], enthalpy)
        enthalpy = np.append(enthalpy, mu_enthalpy[0])
        ind = len(formation)-3
        stoich = np.append([0.0], stoich)
        stoich = np.append(stoich, [1.0])
        structures = np.vstack((stoich, formation)).T
        if self.include_oqmd:
            for ind, doc in enumerate(query.oqmd_cursor):
                oqmd_formation[ind] = doc['enthalpy_per_atom']
                atoms_per_fu = doc['stoichiometry'][0][1] + doc['stoichiometry'][1][1]
                num_fu = (doc['enthalpy']/doc['enthalpy_per_atom']) / float(atoms_per_fu)
                if doc['stoichiometry'][0][0] == one_minus_x_elem:
                    num_b = doc['stoichiometry'][0][1]
                elif doc['stoichiometry'][1][0] == one_minus_x_elem:
                    num_b = doc['stoichiometry'][1][1]
                else:
                    print_failure('Something went wrong!')
                    exit()
                oqmd_enthalpy[ind] = doc['enthalpy'] / (num_b * num_fu)
                oqmd_formation[ind] = doc['enthalpy_per_atom']
                for mu in oqmd_match:
                    for j in range(len(doc['stoichiometry'])):
                        if mu['stoichiometry'][0][0] == doc['stoichiometry'][j][0]:
                            oqmd_formation[ind] -= (mu['enthalpy_per_atom'] *
                                                    doc['stoichiometry'][j][1] /
                                                    atoms_per_fu)
                for elem in doc['stoichiometry']:
                    stoich_string = (str(doc['stoichiometry'][0][0]) +
                                     str(doc['stoichiometry'][0][1]) +
                                     str(doc['stoichiometry'][1][0]) +
                                     str(doc['stoichiometry'][1][1]))
                    if x_elem in elem[0]:
                        oqmd_stoich[ind] = (elem[1])/float(atoms_per_fu)
                oqmd_info.append("{0:^10}\n{1:^24}\n{2:^5s}\n{3:2f} eV".format(
                    stoich_string, 'OQMD' + ' ' + doc['text_id'][0] + ' ' + doc['text_id'][1],
                    doc['space_group'], formation[ind]))
            oqmd_stoich = np.append([0.0], oqmd_stoich)
            oqmd_stoich = np.append(oqmd_stoich, [1.0])
            oqmd_formation = np.append([0.0], oqmd_formation)
            oqmd_formation = np.append(oqmd_formation, [0.0])
            oqmd_enthalpy = np.append(oqmd_mu_enthalpy[1], oqmd_enthalpy)
            oqmd_enthalpy = np.append(oqmd_enthalpy, oqmd_mu_enthalpy[0])
            oqmd_structures = np.vstack((oqmd_stoich, oqmd_formation)).T
            ind = len(oqmd_formation)-3
            for doc in match:
                stoich_string = str(doc['stoichiometry'][0][0]) 
                oqmd_info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(
                    stoich_string, 'OQMD' + ' ' + doc['text_id'][0] + ' ' + doc['text_id'][1],
                    doc['space_group'], oqmd_formation[ind]))
                ind += 1
        # create hull with SciPy routine
        self.hull = ConvexHull(structures)
        if self.include_oqmd:
            oqmd_hull = ConvexHull(oqmd_structures)

        hull_energy = []
        hull_comp = []
        hull_enthalpy = []
        hull_cursor = []
        for ind in range(len(self.hull.vertices)):
            if structures[self.hull.vertices[ind], 1] <= 0:
                hull_energy.append(structures[self.hull.vertices[ind], 1])
                hull_enthalpy.append(enthalpy[self.hull.vertices[ind]])
                hull_comp.append(structures[self.hull.vertices[ind], 0])
        # calculate distance to hull of all structures
        hull_energy = np.asarray(hull_energy)
        hull_enthalpy = np.asarray(hull_enthalpy)
        hull_comp = np.asarray(hull_comp)
        hull_energy = hull_energy[np.argsort(hull_comp)]
        hull_enthalpy = hull_enthalpy[np.argsort(hull_comp)]
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
        stable_enthalpy = list(hull_enthalpy)
        stable_comp = list(hull_comp)
        for ind in range(len(structures)):
            if hull_dist[ind] <= self.hull_cutoff:
                # recolour if under cutoff
                source_ind[ind] = 0
                # get lowest enthalpy at particular comp
                if structures[ind, 0] not in stable_comp:
                    stable_energy.append(structures[ind, 1])
                    stable_enthalpy.append(enthalpy[ind])
                    stable_comp.append(structures[ind, 0])
        # create hull_cursor to pass to other modules
        # skip last and first as they are chem pots
        hull_cursor.append(match[0])
        for ind in range(1, len(hull_dist)-1):
            if hull_dist[ind] <= self.hull_cutoff:
                # take ind-1 to ignore first chem pot
                hull_cursor.append(self.cursor[ind-1])
        hull_cursor.append(match[1])
        if self.include_oqmd:
            oqmd_stable_comp = []
            oqmd_stable_energy = []
            oqmd_stable_enthalpy = []
            for ind in range(len(oqmd_hull.vertices)):
                if oqmd_structures[oqmd_hull.vertices[ind], 1] <= 0:
                    oqmd_stable_energy.append(oqmd_structures[oqmd_hull.vertices[ind], 1])
                    oqmd_stable_enthalpy.append(oqmd_enthalpy[oqmd_hull.vertices[ind]])
                    oqmd_stable_comp.append(oqmd_structures[oqmd_hull.vertices[ind], 0])
            oqmd_stable_comp = np.asarray(oqmd_stable_comp)
            oqmd_stable_energy = np.asarray(oqmd_stable_energy)
            oqmd_stable_enthalpy = np.asarray(oqmd_stable_enthalpy)
            oqmd_stable_energy = oqmd_stable_energy[np.argsort(oqmd_stable_comp)]
            oqmd_stable_enthalpy = oqmd_stable_enthalpy[np.argsort(oqmd_stable_comp)]
            oqmd_stable_comp = oqmd_stable_comp[np.argsort(oqmd_stable_comp)]
        # grab info for datacursor
        info = []
        doc = match[0]
        ind = 0
        stoich_string = str(doc['stoichiometry'][0][0]) 
        info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string,
                                                                doc['text_id'][0] + ' ' +
                                                                doc['text_id'][1],
                                                                doc['space_group'],
                                                                hull_dist[ind]))
        for ind, doc in enumerate(self.cursor):
            stoich_string = str(doc['stoichiometry'][0][0])
            stoich_string += str(doc['stoichiometry'][0][1]) if doc['stoichiometry'][0][1] != 1 else ''
            stoich_string += str(doc['stoichiometry'][1][0])
            stoich_string += str(doc['stoichiometry'][1][1]) if doc['stoichiometry'][1][1] != 1 else ''
            info.append("{0:^10}\n{1:^24}\n{2:^5s}\n{3:2f} eV".format(stoich_string,
                                                                      doc['text_id'][0] + ' ' +
                                                                      doc['text_id'][1],
                                                                      doc['space_group'],
                                                                      hull_dist[ind+1]))
        doc = match[1]
        ind = len(hull_dist)-1
        stoich_string = str(doc['stoichiometry'][0][0]) 
        info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string,
                                                                doc['text_id'][0] + ' ' +
                                                                doc['text_id'][1],
                                                                doc['space_group'],
                                                                hull_dist[ind]))

        stable_energy = np.asarray(stable_energy)
        stable_comp = np.asarray(stable_comp)
        stable_enthalpy = np.asarray(stable_enthalpy)
        stable_energy = stable_energy[np.argsort(stable_comp)]
        stable_enthalpy = stable_enthalpy[np.argsort(stable_comp)]
        stable_comp = stable_comp[np.argsort(stable_comp)]

        self.structures = structures 
        self.info = info
        self.source_ind = source_ind
        self.source_ind = source_ind
        self.hull_cursor = hull_cursor
        self.hull_dist = hull_dist
        self.hull_comp = hull_comp
        self.hull_energy = hull_energy
        self.stable_enthalpy = stable_enthalpy
        self.stable_comp = stable_comp
        self.mu_enthalpy = mu_enthalpy

    def voltage_curve(self, stable_enthalpy, stable_comp, mu_enthalpy):
        """ Take convex hull and calculate voltages. """
        stable_num = []
        stable_comp = stable_comp
        stable_enthalpy = stable_enthalpy
        V = []
        x = []
        for i in range(len(stable_comp)):
            if 1-stable_comp[i] == 0:
                stable_num.append(1e5)
            else:
                stable_num.append(stable_comp[i]/(1-stable_comp[i]))
        # V.append(0)
        for i in range(len(stable_num)-1, 0, -1):
            V.append(-(stable_enthalpy[i] - stable_enthalpy[i-1]) /
                      (stable_num[i] - stable_num[i-1]) +
                      (mu_enthalpy[0]))
            x.append(stable_num[i])
        V.append(V[-1])
        x.append(0)
        self.voltages = V
        self.x = x
        return 

    def plot_hull(self, dis=False):
        """ Plot calculated hull. """
        if self.args.get('png'):
            fig = plt.figure(facecolor=None, figsize=(3,1.5))
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
                                    self.structures[self.hull.vertices[ind], 1]-0.07))
                                                                
        lw = self.scale * 0.05 if self.mpl_new_ver else 1
        # points for off hull structures
        for ind in range(len(self.structures)):
            if self.hull_dist[ind] <= self.hull_cutoff or self.hull_cutoff == 0:
                c = self.colours[self.source_ind[ind]] if self.hull_cutoff == 0 else self.colours[1]
                scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1], s=self.scale*30, lw=lw,
                               alpha=0.9, c=c, edgecolor='k', label=self.info[ind], zorder=100))
            # if dis and warren:
                # ax.plot([self.structures[ind, 0]-disorder[ind]/10, self.structures[ind, 0]],
                        # [self.structures[ind, 1], self.structures[ind, 1]],
                        # c='g', alpha=0.5, lw=0.5)
            # if dis and not warren:
                # ax.plot([self.structures[ind, 0]-disorder[ind]/10, self.structures[ind, 0] + disorder[ind]],
                        # [self.structures[ind, 1], self.structures[ind, 1]],
                        # c='#28B453', alpha=0.5, lw=0.5)
        if self.hull_cutoff != 0:
            c = self.colours[self.source_ind[ind]] if self.hull_cutoff == 0 else self.colours[1]
            ax.scatter(self.structures[1:-1, 0], self.structures[1:-1, 1], s=self.scale*30, lw=lw,
                       alpha=0.3, c=self.colours[-2],
                       edgecolor='k', zorder=10)
        if self.include_oqmd:
            for ind in range(len(oqmd_hull.vertices)):
                if oqmd_structures[oqmd_hull.vertices[ind], 1] <= 0:
                    hull_scatter.append(ax.scatter(oqmd_structures[oqmd_hull.vertices[ind], 0],
                                                   oqmd_structures[oqmd_hull.vertices[ind], 1],
                                                   c=self.colours[-1], marker='*', zorder=10000,
                                                   edgecolor='k',
                                                   s=self.scale*150, lw=1, alpha=1,
                                                   label=oqmd_info[oqmd_hull.vertices[ind]]))
            for ind in range(len(oqmd_stoich)):
                scatter.append(ax.scatter(oqmd_stoich[ind], oqmd_formation[ind], s=self.scale*20, lw=1,
                               alpha=1, c=self.colours[-1], edgecolor='k', marker='D',
                               label=oqmd_info[ind],
                               zorder=200))
        # tie lines
        for ind in range(len(self.hull_comp)-1):
            ax.plot([self.hull_comp[ind], self.hull_comp[ind+1]],
                    [self.hull_energy[ind], self.hull_energy[ind+1]],
                    c=self.colours[0], lw=2, alpha=1, zorder=1000, label='')
            if self.hull_cutoff > 0:
                ax.plot([self.hull_comp[ind], self.hull_comp[ind+1]],
                        [self.hull_energy[ind]+self.hull_cutoff, self.hull_energy[ind+1]+self.hull_cutoff],
                        '--', c=self.colours[1], lw=1, alpha=0.5, zorder=1000, label='')
        if self.include_oqmd:
            for ind in range(len(oqmd_stable_comp)-1):
                ax.plot([oqmd_stable_comp[ind], oqmd_stable_comp[ind+1]],
                        [oqmd_stable_energy[ind], oqmd_stable_energy[ind+1]],
                        c=self.colours[-1], lw=2, alpha=1, zorder=900, label='')
        ax.set_xlim(-0.05, 1.05)
        # data cursor
        if not dis:
            datacursor(scatter[:], formatter='{label}'.format, draggable=False,
                       bbox=dict(fc='white'),
                       arrowprops=dict(arrowstyle='simple', alpha=1))
            # datacursor(hull_scatter[:], formatter='{label}'.format, draggable=False,
                       # bbox=dict(fc='white'),
                       # arrowprops=dict(arrowstyle='simple', alpha=1))
        ax.set_ylim(-0.1 if np.min(self.structures[self.hull.vertices, 1]) > 0
                    else np.min(self.structures[self.hull.vertices, 1])-0.15,
                    0.5 if np.max(self.structures[self.hull.vertices, 1]) > 1
                    else np.max(self.structures[self.hull.vertices, 1])+0.1)
        ax.set_title('$\mathrm{'+x_elem+'_x'+one_minus_x_elem+'_{1-x}}$')
        plt.locator_params(nbins=3)
        ax.set_xlabel('$x$')
        ax.set_xticks([0, 0.33, 0.5, 0.66, 1])
        ax.set_ylabel('$E_\mathrm{F}$ (eV/atom)')
        if self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull.png', dpi=300, bbox_inches='tight')
        else:
            plt.show()
    
    def plot_voltage_curve(self):
        """ Plot calculated voltage curve. """
        if self.args['png']:
            fig = plt.figure(facecolor=None, figsize=(3,1.5))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
        for i in range(2, len(self.voltages)):
            ax.scatter(self.x[i-1], self.voltages[i-1], marker='*', s=100, edgecolor='k', c=self.colours[0], zorder=1000)
            ax.plot([self.x[i], self.x[i]], [self.voltages[i], self.voltages[i-1]], lw=2, c=self.colours[0])
            ax.plot([self.x[i-1], self.x[i]], [self.voltages[i-1], self.voltages[i-1]], lw=2, c=self.colours[0])
        ax.set_ylabel('Voltage (V)')
        plt.locator_params(nbins=4)
        ax.set_xlim(0, np.max(np.asarray(self.x[1:]))+1)
        ax.set_ylim(np.min(np.asarray(self.voltages[2:]))-0.1,np.max(np.asarray(self.voltages[2:]))+0.1)
        ax.set_title('$\mathrm{'+self.elements[0]+'_x'+self.elements[1]+'}$')
        ax.set_xlabel('$x$')
        if self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_voltage.png', dpi=300, bbox_inches='tight')
        else:
            plt.show()

    def subplot_voltage_hull(self, dis=False):
        """ Plot calculated hull with inset voltage curve. """
        if self.args['png']:
            fig = plt.figure(facecolor=None, figsize=(4.5,1.5))
        else:
            fig = plt.figure(facecolor=None)
        ax = plt.subplot2grid((1,3), (0,0), colspan=2)
        ax2 = plt.subplot2grid((1,3), (0,2))
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
                c = self.colours[self.source_ind[ind]] if self.hull_cutoff == 0 else self.colours[1]
                scatter.append(ax.scatter(self.structures[ind, 0], self.structures[ind, 1], s=self.scale*30, lw=lw,
                               alpha=0.9, c=c, edgecolor='k', label=self.info[ind], zorder=100))
            # if dis and warren:
                # ax.plot([self.structures[ind, 0]-disorder[ind]/10, self.structures[ind, 0]],
                        # [self.structures[ind, 1], self.structures[ind, 1]],
                        # c='g', alpha=0.5, lw=0.5)
            # if dis and not warren:
                # ax.plot([self.structures[ind, 0]-disorder[ind]/10, self.structures[ind, 0] + disorder[ind]],
                        # [self.structures[ind, 1], self.structures[ind, 1]],
                        # c='#28B453', alpha=0.5, lw=0.5)
        if self.hull_cutoff != 0:
            c = self.colours[self.source_ind[ind]] if self.hull_cutoff == 0 else self.colours[1]
            ax.scatter(self.structures[1:-1, 0], self.structures[1:-1, 1], s=self.scale*30, lw=lw,
                       alpha=0.3, c=self.colours[-2],
                       edgecolor='k', zorder=10)
        if self.include_oqmd:
            for ind in range(len(oqmd_hull.vertices)):
                if oqmd_structures[oqmd_hull.vertices[ind], 1] <= 0:
                    hull_scatter.append(ax.scatter(oqmd_structures[oqmd_hull.vertices[ind], 0],
                                                   oqmd_structures[oqmd_hull.vertices[ind], 1],
                                                   c=self.colours[-1], marker='*', zorder=10000,
                                                   edgecolor='k',
                                                   s=self.scale*150, lw=1, alpha=1,
                                                   label=oqmd_info[oqmd_hull.vertices[ind]]))
            for ind in range(len(oqmd_stoich)):
                scatter.append(ax.scatter(oqmd_stoich[ind], oqmd_formation[ind], s=self.scale*20, lw=1,
                               alpha=1, c=self.colours[-1], edgecolor='k', marker='D',
                               label=oqmd_info[ind],
                               zorder=200))
        # tie lines
        for ind in range(len(self.hull_comp)-1):
            ax.plot([self.hull_comp[ind], self.hull_comp[ind+1]],
                    [self.hull_energy[ind], self.hull_energy[ind+1]],
                    c=self.colours[0], lw=2, alpha=1, zorder=1000, label='')
            if self.hull_cutoff > 0:
                ax.plot([self.hull_comp[ind], self.hull_comp[ind+1]],
                        [self.hull_energy[ind]+self.hull_cutoff, self.hull_energy[ind+1]+self.hull_cutoff],
                        '--', c=self.colours[1], lw=1, alpha=0.5, zorder=1000, label='')
        if self.include_oqmd:
            for ind in range(len(oqmd_stable_comp)-1):
                ax.plot([oqmd_stable_comp[ind], oqmd_stable_comp[ind+1]],
                        [oqmd_stable_energy[ind], oqmd_stable_energy[ind+1]],
                        c=self.colours[-1], lw=2, alpha=1, zorder=900, label='')
        ax.set_xlim(-0.05, 1.05)
        # data cursor
        if not dis:
            datacursor(scatter[:], formatter='{label}'.format, draggable=False,
                       bbox=dict(fc='white'),
                       arrowprops=dict(arrowstyle='simple', alpha=1))
            # datacursor(hull_scatter[:], formatter='{label}'.format, draggable=False,
                       # bbox=dict(fc='white'),
                       # arrowprops=dict(arrowstyle='simple', alpha=1))
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
            ax2.scatter(self.x[i-1], self.voltages[i-1], marker='*', s=100, edgecolor='k', c=self.colours[0], zorder=1000)
            ax2.plot([self.x[i], self.x[i]], [self.voltages[i], self.voltages[i-1]], lw=2, c=self.colours[0])
            ax2.plot([self.x[i-1], self.x[i]], [self.voltages[i-1], self.voltages[i-1]], lw=2, c=self.colours[0])
        ax2.set_ylabel('Voltage (V)')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_xlim(0, np.max(np.asarray(self.x[1:]))+1)
        ax2.set_ylim(np.min(np.asarray(self.voltages[1:]))-0.1,np.max(np.asarray(self.voltages[1:]))+0.1)
        ax2.set_xlabel('$n_\mathrm{Li}$', labelpad=-3)
        if self.args.get('png'):
            plt.savefig(self.elements[0]+self.elements[1]+'_hull_voltage.png', dpi=300, bbox_inches='tight')
        else:
            plt.show()

    def set_plot_param(self):
        """ Set some plotting options global to
        voltage and hull plots.
        """
        try:
            plt.style.use('bmh')
        except:
            pass
        if self.args.get('png'):
            try:
                plt.style.use('poster')
            except:
                pass
        self.scale = 1
        try:
            c = plt.cm.viridis(np.linspace(0, 1, 100))
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
