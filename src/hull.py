#!/usr/bin/python
# coding: utf-8
from __future__ import print_function
from scipy.spatial import ConvexHull
from mpldatacursor import datacursor
from bson.son import SON
import pymongo as pm
import matplotlib.pyplot as plt
import re
import numpy as np

try:
    plt.style.use('bmh')
except:
    pass


class FryanConvexHull():
    """
    Implements a Convex Hull for formation energies
    from a fryan DBQuery.
    """
    def __init__(self, query):
        """ Accept query from fryan as argument. """
        self.query = query
        self.binary_hull()

    def binary_hull(self, dis=False):
        """ Create a convex hull for two elements. """
        query = self.query
        # local_cursor = query.cursor.clone()
        elements = query.args.get('composition')
        elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements[0]) if elem]
        if len(elements) != 2:
            print('Cannot create binary hull for more or less than 2 elements (yet!).')
            return
        # try to get decent chemical potentials:
        # this relies on all of the first composition query
        # having the same parameters; need to think about this
        mu_enthalpy = np.zeros((2))
        match = [None, None]
        query_dict = []
        for ind, elem in enumerate(elements):
            print('Scanning for suitable', elem, 'chemical potential...')
            self.query.args['composition'] = [elem]
            query_dict.append(dict())
            query_dict[-1]['$and'] = list(self.query.calc_dict['$and'])
            query_dict[-1]['$and'].append(self.query.query_composition(custom_elem=[elem]))
            query_dict[-1]['$and'].append(self.query.query_tags())
            mu_cursor = self.query.repo.find(SON(query_dict[-1])).sort('enthalpy_per_atom',
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
        print('Constructing hull...')

        formation = np.zeros((query.cursor.count()))
        stoich = np.zeros((query.cursor.count()))
        disorder = np.zeros((query.cursor.count()))
        enthalpy = np.zeros((query.cursor.count()))
        info = []
        x_elem = elements[0]
        one_minus_x_elem = elements[1]
        for ind, doc in enumerate(query.cursor):
            atoms_per_fu = doc['stoichiometry'][0][1] + doc['stoichiometry'][1][1]
            num_fu = (doc['enthalpy']/doc['enthalpy_per_atom']) / float(atoms_per_fu)
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
            for mu in match:
                for j in range(len(doc['stoichiometry'])):
                    if mu['stoichiometry'][0][0] == doc['stoichiometry'][j][0]:
                        formation[ind] -= (mu['enthalpy_per_atom'] * doc['stoichiometry'][j][1] /
                                           atoms_per_fu)
            for elem in doc['stoichiometry']:
                stoich_string = (str(doc['stoichiometry'][0][0]) +
                                 str(doc['stoichiometry'][0][1]) +
                                 str(doc['stoichiometry'][1][0]) +
                                 str(doc['stoichiometry'][1][1]))
                if x_elem in elem[0]:
                    stoich[ind] = elem[1]/float(atoms_per_fu)
            info.append("{0:^10}\n{1:^24}\n{2:^5s}\n{3:2f} eV".format(stoich_string,
                                                                      doc['text_id'][0] + ' ' +
                                                                      doc['text_id'][1],
                                                                      doc['space_group'],
                                                                      formation[ind]))
            if dis:
                disorder[ind], warren = self.disorder_hull(doc)
        formation = np.append([0.0], formation)
        formation = np.append(formation, [0.0])
        enthalpy = np.append(mu_enthalpy[1], enthalpy)
        enthalpy = np.append(enthalpy, mu_enthalpy[0])
        ind = len(formation)-3
        for doc in match:
            stoich_string = str(doc['stoichiometry'][0][0]) + str(doc['stoichiometry'][0][1])
            info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string,
                                                                    doc['text_id'][0] + ' ' +
                                                                    doc['text_id'][1],
                                                                    doc['space_group'],
                                                                    formation[ind]))
            ind += 1
        stoich = np.append([0.0], stoich)
        stoich = np.append(stoich, [1.0])
        structures = np.vstack((stoich, formation)).T
        hull = ConvexHull(structures)
        try:
            colours = plt.cm.plasma(np.linspace(0, 1, 100))
        except:
            colours = plt.cm.winter(np.linspace(0, 1, 100))
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        plt.draw()
        # plot all structures
        scatter = []
        for ind in range(len(structures)-2):
            scatter.append(ax.scatter(structures[ind, 0], structures[ind, 1], s=35, lw=1, alpha=1,
                                      c=colours[int(100*structures[ind, 0])], edgecolor='k',
                                      label=info[ind], zorder=100))
            if dis and warren:
                ax.plot([structures[ind, 0]-disorder[ind], structures[ind, 0]],
                        [structures[ind, 1], structures[ind, 1]],
                        c='g', alpha=0.5, lw=0.5)
            if dis and not warren:
                ax.plot([structures[ind, 0]-disorder[ind], structures[ind, 0] + disorder[ind]],
                        [structures[ind, 1], structures[ind, 1]],
                        c='m', alpha=0.5, lw=0.5)
        stable_energy = []
        stable_comp = []
        stable_enthalpy = []
        hull_docs = []
        for ind in range(len(hull.vertices)):
            if structures[hull.vertices[ind], 1] <= 0:
                stable_energy.append(structures[hull.vertices[ind], 1])
                stable_enthalpy.append(enthalpy[hull.vertices[ind]])
                stable_comp.append(structures[hull.vertices[ind], 0])
                scatter.append(ax.scatter(structures[hull.vertices[ind], 0],
                                          structures[hull.vertices[ind], 1],
                                          c='r', marker='*', zorder=1000, edgecolor='k',
                                          s=250, lw=1, alpha=1, label=info[hull.vertices[ind]]))
        # skip last and first as they are chem pots
        for ind in range(1, len(hull.vertices)-1):
            query.cursor.rewind()
            hull_docs.append(query.cursor[int(np.sort(hull.vertices)[ind])])
        stable_energy = np.asarray(stable_energy)
        stable_comp = np.asarray(stable_comp)
        stable_enthalpy = np.asarray(stable_enthalpy)
        stable_energy = stable_energy[np.argsort(stable_comp)]
        stable_enthalpy = stable_enthalpy[np.argsort(stable_comp)]
        stable_comp = stable_comp[np.argsort(stable_comp)]
        for ind in range(len(stable_comp)-1):
                ax.plot([stable_comp[ind], stable_comp[ind+1]],
                        [stable_energy[ind], stable_energy[ind+1]],
                        'k--', lw=2, alpha=1, zorder=1, label='')
        ax.set_xlim(-0.05, 1.05)
        if not dis:
            datacursor(scatter[:], formatter='{label}'.format, draggable=False,
                       bbox=dict(fc='yellow'),
                       arrowprops=dict(arrowstyle='simple', alpha=1))
        ax.set_ylim(-0.1 if np.min(structures[hull.vertices, 1]) > 0
                    else np.min(structures[hull.vertices, 1])-0.1,
                    0.5 if np.max(structures[hull.vertices, 1]) > 1
                    else np.max(structures[hull.vertices, 1])+0.1)
        ax.set_title('$\mathrm{'+str(x_elem)+'_x'+str(one_minus_x_elem)+'_{1-x}}$')
        ax.set_xlabel('$x$')
        ax.set_ylabel('formation enthalpy per atom (eV)')
        if query.args.get('voltage'):
            print('Generating voltage curve...')
            self.voltage_curve(stable_enthalpy, stable_comp, mu_enthalpy, elements)
        plt.show()
        self.hull_docs = hull_docs
        return hull_docs

    def voltage_curve(self, stable_enthalpy, stable_comp, mu_enthalpy, elements):
        """ Take convex hull and plot voltage curves. """
        stable_num = []
        V = []
        x = []
        for i in range(len(stable_comp)):
            stable_num.append(stable_comp[i]/(1-stable_comp[i]))
        V.append(0)
        x.append(1e5)
        for i in range(len(stable_num)-2, 0, -1):
            V.append(-(stable_enthalpy[i] - stable_enthalpy[i-1]) /
                      (stable_num[i] - stable_num[i-1]) +
                      (mu_enthalpy[0]))
            x.append(stable_num[i])
        V.append(V[-1])
        x.append(0)
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        plt.style.use('bmh')
        colour = []
        try:
            colour.append(ax._get_lines.prop_cycler.next()['color'])
        except:
            colour.append('blue')
        # print(zip(x, V))
        for i in range(1, len(V)):
            ax.scatter(x[i], V[i],
                       marker='*', c=colour[0], zorder=1000, edgecolor='k', s=200, lw=1)
        for i in range(2, len(V)):
            ax.scatter(x[i], V[i-1],
                       marker='*', c=colour[0], zorder=1000, edgecolor='k', s=200, lw=1)
            ax.plot([x[i], x[i]], [V[i], V[i-1]], lw=2, c=colour[0])
            ax.plot([x[i-1], x[i]], [V[i-1], V[i-1]], lw=2, c=colour[0])
        ax.set_ylabel('V')
        ax.set_xlim(0)
        ax.set_title('$\mathrm{'+elements[0]+'_x'+elements[1]+'}$')
        ax.set_xlabel('$x$')
        plt.show()
