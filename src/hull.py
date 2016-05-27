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
        x_elem = elements[0]
        one_minus_x_elem = elements[1]
        for ind, doc in enumerate(query.cursor):
            atoms_per_fu = doc['stoichiometry'][0][1] + doc['stoichiometry'][1][1]  
            formation[ind] = doc['enthalpy_per_atom']
            for mu in match:
                for j in range(len(doc['stoichiometry'])):
                    if mu['stoichiometry'][0][0] == doc['stoichiometry'][j][0]:
                        formation[ind] -= mu['enthalpy_per_atom']*doc['stoichiometry'][j][1] / atoms_per_fu
            for elem in doc['stoichiometry']:
                stoich_string = str(doc['stoichiometry'][0][0]) + str(doc['stoichiometry'][0][1]) + str(doc['stoichiometry'][1][0]) + str(doc['stoichiometry'][1][1])
                if x_elem in elem[0]:
                    stoich[ind] = elem[1]/float(atoms_per_fu)
            info.append("{0:^10}\n{1:^24}\n{2:^5s}\n{3:2f} eV".format(stoich_string, doc['text_id'][0]+' '+doc['text_id'][1], doc['space_group'], formation[ind]))
            if dis:
                disorder[ind], warren = self.disorder_hull(doc)
        formation = np.append(formation, [0.0, 0.0])
        ind = len(formation)-3
        for doc in match:
            stoich_string = str(doc['stoichiometry'][0][0]) + str(doc['stoichiometry'][0][1]) 
            info.append("{0:^10}\n{1:24}\n{2:5s}\n{3:2f} eV".format(stoich_string, doc['text_id'][0]+' '+doc['text_id'][1], doc['space_group'], formation[ind]))
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
        scatter = []
        for ind in range(len(structures)-2):
            scatter.append(ax.scatter(structures[ind,0], structures[ind,1], s=35, lw=1, alpha=1, c=colours[int(100*structures[ind,0])], edgecolor='k', label=info[ind], zorder=100))
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
                scatter.append(ax.scatter(structures[hull.vertices[ind], 0], structures[hull.vertices[ind], 1], 
                           c='r', marker='*', zorder=1000, edgecolor='k', s=250, lw=1, alpha=1, label=info[hull.vertices[ind]]))
        stable_energy =  np.asarray(stable_energy)
        stable_comp =  np.asarray(stable_comp)
        stable_energy = stable_energy[np.argsort(stable_comp)]
        stable_comp =  stable_comp[np.argsort(stable_comp)]
        for ind in range(len(stable_comp)-1):
                ax.plot([stable_comp[ind], stable_comp[ind+1]], 
                        [stable_energy[ind], stable_energy[ind+1]],
                        'k--', lw=2, alpha=1, zorder=1, label='')
        ax.set_xlim(-0.05, 1.05)
        if not dis:
            datacursor(scatter[:], formatter='{label}'.format, draggable=False, bbox=dict(fc='yellow'),
                    arrowprops=dict(arrowstyle='simple', alpha=1))
        ax.set_ylim(-0.1 if np.min(structures[hull.vertices,1]) > 0 else np.min(structures[hull.vertices,1])-0.1,
                    0.5 if np.max(structures[hull.vertices,1]) > 1 else np.max(structures[hull.vertices,1])+0.1)
        ax.set_title('$\mathrm{'+str(x_elem)+'_x'+str(one_minus_x_elem)+'_{1-x}}$')
        ax.set_xlabel('$x$')
        ax.set_ylabel('formation enthalpy per atom (eV)')
        if query.args.get('voltage'):
            print('Generating voltage curve...')
            self.voltage_curve(stable_energy, stable_comp, mu_enthalpy, elements)# info[stable_list])
        plt.show()
        return structures, disorder, hull, mu, info, fig
    
    def voltage_curve(self, stable_energy, stable_comp, mu_enthalpy, elements):
        """ Take convex hull and plot voltage curves. """
        voltages = np.zeros((len(stable_comp))-1)
        for i in range(1,len(voltages)):
            voltages[i] = -(stable_energy[i] - stable_energy[i-1])
            voltages[i] /= stable_comp[i] - stable_comp[i-1]
            voltages[i] += mu_enthalpy[0]
        # voltages /= (188.01/2)
        # print(voltages)
        voltages[0] = np.NaN
        fig = plt.figure(facecolor='w')
        ax  = fig.add_subplot(111)
        # datacursor(formatter='{label}'.format, draggable=False, bbox=dict(fc='white'),
            # arrowprops=dict(arrowstyle='simple', alpha=0.8))
        plt.style.use('bmh')
        colour = []
        colour.append(ax._get_lines.prop_cycler.next()['color'])
        # ax.plot([stable_comp[0], stable_comp[0]], [voltages[0], voltages[1]], ls='-.', c=colour[0], lw=2)
        # ax.plot([stable_comp[-2], stable_comp[-2]], [voltages[-1], 0], ls='-.', c=colour[0], lw=2)
        for i in range(1,len(voltages)-1):
            ax.plot([stable_comp[i], stable_comp[i]], [voltages[i], voltages[i+1]], c=colour[0], lw=2)
        for i in range(1, len(voltages)):
            ax.plot([stable_comp[i-1], stable_comp[i]], [voltages[i], voltages[i]], c=colour[0], lw=2)
            ax.scatter(stable_comp[i], voltages[i],
                       c=colour[0], marker='*', zorder=1000, edgecolor='k', s=250, lw=1, alpha=1)
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylabel('V')
        ax.set_title('$\mathrm{'+elements[0]+'_x'+elements[1]+'_{1-x}}$')
        ax.set_xlabel('$x$')
        ax.set_yticks([])
        # ax.set_ylim(0,np.max(voltages[1:])*1.5)
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
