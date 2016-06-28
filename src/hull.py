#!/usr/bin/python
# coding: utf-8
""" This file implements convex hull functionality
from database queries.
"""
from __future__ import print_function
from scipy.spatial import ConvexHull
from mpldatacursor import datacursor
from bson.son import SON
from bisect import bisect_left
import pymongo as pm
import matplotlib.pyplot as plt
import re
import numpy as np

try:
    plt.style.use('bmh')
    plt.style.use('presentation')
except:
    pass


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
        self.binary_hull()

    def binary_hull(self, dis=False):
        """ Create a convex hull for two elements. """
        query = self.query
        include_oqmd = query.args.get('include_oqmd')
        elements = query.args.get('composition')
        elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements[0]) if elem]
        if len(elements) != 2:
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
        for ind, elem in enumerate(elements):
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
                print('No possible chem pots found for', elem, '.')
                return
        # include OQMD structures if desired, first find chem pots
        if include_oqmd:
            oqmd_mu_enthalpy = np.zeros((2))
            oqmd_match = [None, None]
            oqmd_query_dict = dict()
            for ind, elem in enumerate(elements):
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
                    print('Using', ''.join([oqmd_match[ind]['text_id'][0], ' ',
                          match[ind]['text_id'][1]]), 'as OQMD chem pot for', elem)
                    print(60*'─')
                else:
                    print('No possible chem pots found for', elem, '.')
                    return
        print('Constructing hull...')
        num_structures = len(self.cursor)
        formation = np.zeros((num_structures))
        stoich = np.zeros((num_structures))
        enthalpy = np.zeros((num_structures))
        disorder = np.zeros((num_structures))
        source_ind = np.zeros((num_structures+2), dtype=int)
        hull_dist = np.zeros((num_structures+2))
        info = []
        source_list = []
        if include_oqmd:
            oqmd_num_structures = query.oqmd_cursor.count()
            oqmd_formation = np.zeros((oqmd_num_structures))
            oqmd_stoich = np.zeros((oqmd_num_structures))
            oqmd_enthalpy = np.zeros((oqmd_num_structures))
            oqmd_info = []
        if dis:
            from disorder import disorder_hull
        # define hull by order in command-line arguments
        x_elem = elements[0]
        one_minus_x_elem = elements[1]
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
                exit('Something went wrong!')
            # get enthalpy per unit B
            enthalpy[ind] = doc['enthalpy'] / (num_b*num_fu)
            formation[ind] = doc['enthalpy_per_atom']
            source_dir = ''.join(doc['source'][0].split('/')[:-1])
            if source_dir in source_list:
                source_ind[ind] = source_list.index(source_dir) + 1
            else:
                source_list.append(source_dir)
                source_ind[ind] = source_list.index(source_dir) + 1
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
        if include_oqmd:
            for ind, doc in enumerate(query.oqmd_cursor):
                oqmd_formation[ind] = doc['enthalpy_per_atom']
                atoms_per_fu = doc['stoichiometry'][0][1] + doc['stoichiometry'][1][1]
                num_fu = (doc['enthalpy']/doc['enthalpy_per_atom']) / float(atoms_per_fu)
                if doc['stoichiometry'][0][0] == one_minus_x_elem:
                    num_b = doc['stoichiometry'][0][1]
                elif doc['stoichiometry'][1][0] == one_minus_x_elem:
                    num_b = doc['stoichiometry'][1][1]
                else:
                    exit('Something went wrong!')
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
                stoich_string = str(doc['stoichiometry'][0][0]) + str(doc['stoichiometry'][0][1])
                oqmd_info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(
                    stoich_string, 'OQMD' + ' ' + doc['text_id'][0] + ' ' + doc['text_id'][1],
                    doc['space_group'], oqmd_formation[ind]))
                ind += 1
        # create hull with SciPy routine
        hull = ConvexHull(structures)
        if include_oqmd:
            oqmd_hull = ConvexHull(oqmd_structures)

        hull_energy = []
        hull_comp = []
        hull_enthalpy = []
        hull_cursor = []
        for ind in range(len(hull.vertices)):
            if structures[hull.vertices[ind], 1] <= 0:
                hull_energy.append(structures[hull.vertices[ind], 1])
                hull_enthalpy.append(enthalpy[hull.vertices[ind]])
                hull_comp.append(structures[hull.vertices[ind], 0])
        # calculate distance to hull of all structures
        hull_energy = np.asarray(hull_energy)
        hull_enthalpy = np.asarray(hull_energy)
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
        if include_oqmd:
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
        stoich_string = str(doc['stoichiometry'][0][0]) + str(doc['stoichiometry'][0][1])
        info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string,
                                                                doc['text_id'][0] + ' ' +
                                                                doc['text_id'][1],
                                                                doc['space_group'],
                                                                hull_dist[ind]))
        for ind, doc in enumerate(self.cursor):
            stoich_string = (str(doc['stoichiometry'][0][0]) +  str(doc['stoichiometry'][0][1]) +
                             str(doc['stoichiometry'][1][0]) + str(doc['stoichiometry'][1][1]))
            info.append("{0:^10}\n{1:^24}\n{2:^5s}\n{3:2f} eV".format(stoich_string,
                                                                      doc['text_id'][0] + ' ' +
                                                                      doc['text_id'][1],
                                                                      doc['space_group'],
                                                                      hull_dist[ind+1]))
        doc = match[1]
        ind = len(hull_dist)-1
        stoich_string = str(doc['stoichiometry'][0][0]) + str(doc['stoichiometry'][0][1])
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

        # PLOTTING ONLY
        fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
        scale = 1
        try:
            c = plt.cm.viridis(np.linspace(0, 1, 100))
            mpl_new_ver = True
        except:
            mpl_new_ver = False
        from palettable.colorbrewer.qualitative import Dark2_8
        from palettable.colorbrewer.qualitative import Set3_10
        if len(source_list) < 6:
            colours = Dark2_8.hex_colors[1:len(source_list)+1]
        # first colour reserved for hull
        colours.insert(0, Dark2_8.hex_colors[0])
        # penultimate colour reserved for off hull above cutoff
        colours.append(Dark2_8.hex_colors[-1])
        # last colour reserved for OQMD
        colours.append(Set3_10.hex_colors[-1])
        scatter = []
        hull_scatter = []
        plt.draw()
        # structures on hull
        for ind in range(len(hull.vertices)):
            if structures[hull.vertices[ind], 1] <= 0:
                hull_scatter.append(ax.scatter(structures[hull.vertices[ind], 0],
                                               structures[hull.vertices[ind], 1],
                                               c=colours[source_ind[hull.vertices[ind]]],
                                               marker='*', zorder=99999, edgecolor='k',
                                               s=scale*150, lw=1, alpha=1,
                                               label=info[hull.vertices[ind]]))
        lw = scale * 0.05 if mpl_new_ver else 1
        # off hull structures
        for ind in range(len(structures)):
            if hull_dist[ind] <= self.hull_cutoff or self.hull_cutoff == 0:
                c = colours[source_ind[ind]] if self.hull_cutoff == 0 else colours[1]
                # print(c)
                scatter.append(ax.scatter(structures[ind, 0], structures[ind, 1], s=scale*30, lw=lw,
                               alpha=0.9, c=c, edgecolor='#a64902', label=info[ind], zorder=100))
            # if dis and warren:
                # ax.plot([structures[ind, 0]-disorder[ind]/10, structures[ind, 0]],
                        # [structures[ind, 1], structures[ind, 1]],
                        # c='g', alpha=0.5, lw=0.5)
            # if dis and not warren:
                # ax.plot([structures[ind, 0]-disorder[ind]/10, structures[ind, 0] + disorder[ind]],
                        # [structures[ind, 1], structures[ind, 1]],
                        # c='#28B453', alpha=0.5, lw=0.5)
        if self.hull_cutoff != 0:
            c = colours[source_ind[ind]] if self.hull_cutoff == 0 else colours[1]
            ax.scatter(structures[1:-1, 0], structures[1:-1, 1], s=scale*30, lw=lw,
                       alpha=0.3, c=colours[-2],
                       edgecolor='k', zorder=10)
        if include_oqmd:
            for ind in range(len(oqmd_hull.vertices)):
                if oqmd_structures[oqmd_hull.vertices[ind], 1] <= 0:
                    hull_scatter.append(ax.scatter(oqmd_structures[oqmd_hull.vertices[ind], 0],
                                                   oqmd_structures[oqmd_hull.vertices[ind], 1],
                                                   c=colours[-1], marker='*', zorder=10000,
                                                   edgecolor='k',
                                                   s=scale*150, lw=1, alpha=1,
                                                   label=oqmd_info[oqmd_hull.vertices[ind]]))
            for ind in range(len(oqmd_stoich)):
                scatter.append(ax.scatter(oqmd_stoich[ind], oqmd_formation[ind], s=scale*20, lw=1,
                               alpha=1, c=colours[-1], edgecolor='k', marker='D',
                               label=oqmd_info[ind],
                               zorder=200))
        # tie lines
        for ind in range(len(hull_comp)-1):
            ax.plot([hull_comp[ind], hull_comp[ind+1]],
                    [hull_energy[ind], hull_energy[ind+1]],
                    c=colours[0], lw=2, alpha=1, zorder=1000, label='')
            if self.hull_cutoff > 0:
                ax.plot([hull_comp[ind], hull_comp[ind+1]],
                        [hull_energy[ind]+self.hull_cutoff, hull_energy[ind+1]+self.hull_cutoff],
                        '--', c=colours[1], lw=1, alpha=0.5, zorder=1000, label='')
        if include_oqmd:
            for ind in range(len(oqmd_stable_comp)-1):
                ax.plot([oqmd_stable_comp[ind], oqmd_stable_comp[ind+1]],
                        [oqmd_stable_energy[ind], oqmd_stable_energy[ind+1]],
                        c=colours[-1], lw=2, alpha=1, zorder=900, label='')
        ax.set_xlim(-0.05, 1.05)
        # data cursor
        if not dis:
            datacursor(scatter[:], formatter='{label}'.format, draggable=False,
                       bbox=dict(fc='white'),
                       arrowprops=dict(arrowstyle='simple', alpha=1))
            # datacursor(hull_scatter[:], formatter='{label}'.format, draggable=False,
                       # bbox=dict(fc='white'),
                       # arrowprops=dict(arrowstyle='simple', alpha=1))
        ax.set_ylim(-0.1 if np.min(structures[hull.vertices, 1]) > 0
                    else np.min(structures[hull.vertices, 1])-0.1,
                    0.5 if np.max(structures[hull.vertices, 1]) > 1
                    else np.max(structures[hull.vertices, 1])+0.1)
        ax.set_title('$\mathrm{'+str(x_elem)+'_x'+str(one_minus_x_elem)+'_{1-x}}$')
        ax.set_xlabel('$x$')
        ax.set_ylabel('formation enthalpy per atom (eV)')
        if self.args['subcmd'] == 'voltage':
            print('Generating voltage curve...')
            self.voltage_curve(stable_enthalpy, stable_comp, mu_enthalpy, elements)
        plt.show()
        # plt.savefig('LiAs_hull.png', dpi=500)
        self.hull_cursor = hull_cursor

    def voltage_curve(self, stable_enthalpy, stable_comp, mu_enthalpy, elements):
        """ Take convex hull and plot voltage curves. """
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
            V.append(-(stable_enthalpy[i] - stable_enthalpy[i-1]) /
                      (stable_num[i] - stable_num[i-1]) +
                      (mu_enthalpy[0]))
            x.append(stable_num[i])
        V.append(V[-1])
        x.append(0)
        V[0] = 0
        fig = plt.figure(facecolor='w', figsize=(3,2))
        ax = fig.add_subplot(111)
        colour = []
        try:
            colour.append(ax._get_lines.prop_cycler.next()['color'])
        except:
            colour.append('blue')
        for i in range(0, len(V)):
            ax.scatter(x[i], V[i],
                       marker='*', c=colour[0], zorder=1000, edgecolor='k', s=200, lw=1)
        for i in range(1, len(V)):
            ax.scatter(x[i], V[i-1],
                       marker='*', c=colour[0], zorder=1000, edgecolor='k', s=200, lw=1)
            ax.plot([x[i], x[i]], [V[i], V[i-1]], lw=2, c=colour[0])
            ax.plot([x[i-1], x[i]], [V[i-1], V[i-1]], lw=2, c=colour[0])
        ax.set_ylabel('V')
        ax.set_xlim(0, np.max(np.asarray(x[1:]))+1)
        # ax.set_ylim(0)
        ax.set_title('$\mathrm{'+elements[0]+'_x'+elements[1]+'}$')
        ax.set_xlabel('$x$')
        plt.show()
