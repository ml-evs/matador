#!/usr/bin/python
# coding: utf-8
from __future__ import print_function
# import related crysdb functionality
# from fryan import DBQuery
# import external libraries
from scipy.spatial import ConvexHull
from mpldatacursor import datacursor
from bson.son import SON
import pymongo as pm
import matplotlib.pyplot as plt
try:
    plt.style.use('bmh')
except:
    pass
import re
import numpy as np

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
        elements = query.args.get('composition')
        elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements[0]) if elem]
        if len(elements) != 2:
            print('Cannot create binary hull for more or less than 2 elements.')
            return
        cursor = self.query.cursor.clone()
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
            mu_cursor = self.query.repo.find(SON(query_dict[-1])).sort('enthalpy_per_atom', pm.ASCENDING)
            for doc_ind, doc in enumerate(mu_cursor):
                if doc_ind == 0:
                    match[ind] = doc
                    break
            if match[ind] != None:
                mu_enthalpy[ind] = float(match[ind]['enthalpy_per_atom'])
                print('Using', ''.join([match[ind]['text_id'][0], ' ', match[ind]['text_id'][1]]), 'as chem pot for', elem)
                print(60*'─')
            else:
                print('No possible chem pots found for', elem, '.')
                return
        print('Plotting hull...')
        formation = np.zeros((query.cursor.count()))
        stoich = np.zeros((query.cursor.count()))
        disorder = np.zeros((query.cursor.count()))
        info = []
        one_minus_x_elem = ''
        for ind, doc in enumerate(query.cursor):
            atoms_per_fu = doc['stoichiometry'][0][1] + doc['stoichiometry'][1][1]  
            # this is probably better done by spatula; then can plot hull for given chem pot
            formation[ind] = doc['enthalpy_per_atom']
            for mu in match:
                for j in range(len(doc['stoichiometry'])):
                    if mu['stoichiometry'][0][0] == doc['stoichiometry'][j][0]:
                        formation[ind] -= mu['enthalpy_per_atom']*doc['stoichiometry'][j][1] / atoms_per_fu
            stoich[ind] = doc['stoichiometry'][1][1]/float(atoms_per_fu)
            if one_minus_x_elem != doc['stoichiometry'][1][0] and len(one_minus_x_elem) != 0:
                print('A problem has occurred...')
            one_minus_x_elem = doc['stoichiometry'][1][0]
            x_elem = doc['stoichiometry'][0][0]
            info.append("{0:^24}\n{1:5s}\n{2:2f} eV\n{3:^10}\n{4:^24}".format(doc['text_id'][0]+' '+doc['text_id'][1], doc['space_group'], formation[ind], doc['stoichiometry'], doc['source'][0].split('/')[-1]))
            if dis:
                disorder[ind], warren = self.disorder_hull(doc)
        formation = np.append(formation, [0.0, 0.0])
        ind = len(formation)-3
        for doc in match:
            info.append("{0:^24}\n{1:5s}\n{2:2f} eV\n{3:^10}\n{4:^24}".format(doc['text_id'][0]+' '+doc['text_id'][1], doc['space_group'], formation[ind], doc['stoichiometry'], doc['source'][0].split('/')[-1]))
            ind += 1 
        stoich = np.append(stoich, [0.0, 1.0])
        structures = np.vstack((stoich, formation)).T
        hull = ConvexHull(structures)
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        plt.draw()
        try:
            colours = plt.cm.plasma(np.linspace(0, 1, 100))
        except:
            colours = plt.cm.winter(np.linspace(0, 1, 100))
        for ind in range(len(structures)-2):
            ax.scatter(structures[ind,0], structures[ind,1], s=50, lw=1, alpha=1, c=colours[int(100*structures[ind,0])], edgecolor='k', label=info[ind], zorder=100)
            # if dis and warren:
                # ax.plot([structures[ind,0]-disorder[ind], structures[ind,0]], [structures[ind,1], structures[ind,1]],
                        # c='g', alpha=0.5, lw=0.5)
            # if dis and not warren:
                # ax.plot([structures[ind,0]-disorder[ind], structures[ind,0]+disorder[ind]], [structures[ind,1], structures[ind,1]],
                        # c='m', alpha=0.5, lw=0.5)
        stable_energy = []
        stable_comp = []
        for ind in range(len(hull.vertices)):
            if structures[hull.vertices[ind], 1] <= 0:
                stable_energy.append(structures[hull.vertices[ind], 1])
                stable_comp.append(structures[hull.vertices[ind], 0])
                ax.scatter(structures[hull.vertices[ind], 0], structures[hull.vertices[ind], 1], 
                           c='r', marker='*', zorder=1000, edgecolor='k', s=250, lw=1, alpha=1, label=info[ind])
        for ind in range(len(hull.vertices)-1):
            if structures[hull.vertices[ind+1], 1] <= 0 and structures[hull.vertices[ind], 1] <= 0:
                ax.plot([structures[hull.vertices[ind], 0], structures[hull.vertices[ind+1], 0]], 
                        [structures[hull.vertices[ind], 1], structures[hull.vertices[ind+1], 1]],
                        'k--', lw=2, alpha=1, zorder=1)
        stable_energy =  np.asarray(stable_energy)
        stable_comp =  np.asarray(stable_comp)
        ax.set_xlim(-0.05, 1.05)
        if not dis:
            datacursor(formatter='{label}'.format, draggable=False, bbox=dict(fc='white'),
                    arrowprops=dict(arrowstyle='simple', alpha=0.8))
        ax.set_ylim(-0.1 if np.min(structures[hull.vertices,1]) > 0 else np.min(structures[hull.vertices,1])-0.1,
                    0.5 if np.max(structures[hull.vertices,1]) > 1 else np.max(structures[hull.vertices,1])+0.1)
        ax.set_title('$\mathrm{'+str(x_elem)+'_x'+str(one_minus_x_elem)+'_{1-x}}$')
        ax.set_xlabel('$x$')
        ax.set_ylabel('formation enthalpy per atom (eV)')
        if query.args.get('voltage'):
            print('Generating voltage curve...')
            print(mu_enthalpy)
            self.voltage_curve(stable_energy, stable_comp, mu_enthalpy)# info[stable_list])
        plt.show()
        return structures, disorder, hull, mu, info, fig
    
    def voltage_curve(self, stable_energy, stable_comp, mu_enthalpy):
        """ Take convex hull and plot voltage curves. """
        voltages = np.zeros((len(stable_comp)))
        for i in range(len(voltages)-1):
            voltages[i] = stable_energy[i+1] - stable_energy[i]
            voltages[i] /= stable_comp[i+1] - stable_comp[i]
            voltages[i] -= mu_enthalpy[1]
            i += 1
        voltages[-1] = voltages[-2]
        fig = plt.figure(facecolor='w')
        ax  = fig.add_subplot(111)
        # datacursor(formatter='{label}'.format, draggable=False, bbox=dict(fc='white'),
            # arrowprops=dict(arrowstyle='simple', alpha=0.8))
        for i in range(len(voltages)-1):
            # ax.plot([1-stable_comp[i], 1-stable_comp[i]], [voltages[i], voltages[i]])
            # ax.plot([1-stable_comp[i], 1-stable_comp[i+1]], [voltages[i], voltages[i+1]])
            ax.scatter(1-stable_comp[i], voltages[i], s=100)#,# label=info[i])
            ax.scatter(1-stable_comp[i], voltages[i+1], s=100)#,# label=info[i])
        plt.show()

def disorder_hull(doc):
    """ Broaden structures on phase diagram by 
    a measure of local stoichiometry.
    """
    num_atoms = doc['num_atoms']
    lat_cart  = doc['lattice_cart']
    disps = np.zeros((num_atoms, num_atoms-1))
    atoms = np.empty((num_atoms, num_atoms-1), dtype=str)
    for i in range(num_atoms):
        jindex = 0
        for j in range(num_atoms):
            temp_disp = np.zeros((3))
            real_disp = np.zeros((3))
            if i != j:
                atoms[i, jindex] = doc['atom_types'][j]
                for k in range(3):
                    temp_disp[k] = (doc['positions_frac'][j][k] - doc['positions_frac'][i][k])
                    if temp_disp[k] > 0.5:
                        temp_disp[k] -= 1
                    elif temp_disp[k] < -0.5:
                        temp_disp[k] += 1
                for k in range(3):
                    for q in range(3):
                        real_disp[q] += temp_disp[k]*lat_cart[k][q]
                for k in range(3):
                    disps[i,jindex] += real_disp[k]**2
                jindex += 1
    disps = np.sqrt(disps)
    
    def warren_cowley(atoms, disps):
        nn_atoms = []
        for i in range(len(atoms)):
            nn_atoms.append(atoms[i][np.where(disps[i] < 3)])
        count = np.zeros((2), dtype=float)
        for i in range(len(nn_atoms)):
            same_elem = doc['atom_types'][i][0]
            for j in range(len(nn_atoms[i])):
                if nn_atoms[i][j] == same_elem:
                    count[0] += 1.0
                    count[1] += 1.0
                else:
                    count[0] -= 1.0
                    count[1] += 1.0
        return count[0] / (4*(count[1]))
    
    def bond_disorder(atoms, disps):
        nn_atoms = []
        for i in range(len(atoms)):
            nn_atoms.append(atoms[i][np.where(disps[i] < 3)])
        count = np.zeros((2), dtype=float)
        for i in range(len(nn_atoms)):
            same_elem = doc['atom_types'][i][0]
            for j in range(len(nn_atoms[i])):
                if nn_atoms[i][j] == same_elem:
                    count[0] += 1.0
                else:
                    count[1] += 1.0
        
        return count[0] / (4*(count[1]+count[0]))
    
    warren = False
    if warren:
        return warren_cowley(atoms, disps), warren
    else:
        return bond_disorder(atoms, disps), warren
