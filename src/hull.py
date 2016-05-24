#!/usr/bin/python
# coding: utf-8
from __future__ import print_function
# import related crysdb functionality
# from fryan import DBQuery
# import external libraries
from scipy.spatial import ConvexHull
from mpldatacursor import datacursor
import matplotlib.pyplot as plt
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
            print('Cannot create binary hull for more than 2 elements.')
            return
        # try to get decent chemical potentials:
        # this relies on all of the first composition query
        # having the same parameters; need to think about this
        mu = np.array([0.0, 0.0])
        match = [None, None]
        for ind, elem in enumerate(elements):
            print('Scanning for suitable', elem, 'chemical potential...')
            query.args['composition'] = [elem]
            mu_cursor = query.query_composition()
            if mu_cursor.count() > 0:
                for doc in mu_cursor:
                    try:
                        # temporary fix to pspot from dir issue
                        if type(doc['species_pot']) != dict or elem not in doc['species_pot'] or not '.usp' in doc['species_pot'][elem]:
                            continue
                        else:
                            print('\n', doc['species_pot'][elem], 'vs', query.cursor[0]['species_pot'][elem], end='')
                            if doc['species_pot'][elem] == query.cursor[0]['species_pot'][elem]:
                                print(' ✓')
                                print('\n\t', doc['external_pressure'][0][0], 'GPa vs', query.cursor[0]['external_pressure'][0][0], 'GPa', end='')
                                if doc['external_pressure'][0] == query.cursor[0]['external_pressure'][0]:
                                    print(' ✓')
                                    print('\n\t\t', doc['xc_functional'], 'vs', query.cursor[0]['xc_functional'], end='')
                                    if doc['xc_functional'] == query.cursor[0]['xc_functional']:
                                        print(' ✓')
                                        print('\n\t\t\t', doc['cut_off_energy'], 'eV vs', query.cursor[0]['cut_off_energy'], 'eV', end='')
                                        if doc['cut_off_energy'] >= query.cursor[0]['cut_off_energy']:
                                            print(' ✓')
                                            print(60*'─')
                                            match[ind] = doc
                                            print('Match found!')
                                            break
                    except Exception as oops:
                        print(oops)
                        continue
                if match[ind] != None:
                    mu[ind] = float(match[ind]['enthalpy_per_atom'])
                    print('Using', ''.join([match[ind]['text_id'][0], ' ', match[ind]['text_id'][1]]), 'as chem pot for', elem)
                    print(60*'─')
                else:
                    print('No possible chem pots found for', elem, '.')
                    return
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
        points = np.vstack((stoich, formation)).T
        hull = ConvexHull(points)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for ind in range(len(points)-2):
            ax.scatter(points[ind,0], points[ind,1], s=50, lw=1, alpha=0.6, label=info[ind], zorder=100)
            # if dis and warren:
                # ax.plot([points[ind,0]-disorder[ind], points[ind,0]], [points[ind,1], points[ind,1]],
                        # c='g', alpha=0.5, lw=0.5)
            # if dis and not warren:
                # ax.plot([points[ind,0]-disorder[ind], points[ind,0]+disorder[ind]], [points[ind,1], points[ind,1]],
                        # c='m', alpha=0.5, lw=0.5)
        for ind in range(len(hull.vertices)):
            if points[hull.vertices[ind], 1] <= 0:
                ax.scatter(points[hull.vertices[ind], 0], points[hull.vertices[ind], 1], 
                           c='r', marker='*', zorder=1000, s=250, lw=1, alpha=1, label=info[ind])
        for ind in range(len(hull.vertices)-1):
            if points[hull.vertices[ind+1], 1] <= 0:
                ax.plot([points[hull.vertices[ind], 0], points[hull.vertices[ind+1], 0]], 
                        [points[hull.vertices[ind], 1], points[hull.vertices[ind+1], 1]],
                        'k--', lw=1, alpha=0.6, zorder=1)
        ax.set_xlim(-0.05, 1.05)
        if not dis:
            datacursor(formatter='{label}'.format, draggable=False)
        # for i in range(len(hull.vertices)):
            # if points[hull.vertices[i], 1] <= 0:
        ax.set_ylim(-0.1 if np.min(points[:,1]) < 0 else np.min(points[:,1])-0.1,
                    1 if np.max(points[:,1])>1 else np.max(points[:,1])+0.1)
        ax.set_title('$\mathrm{'+str(x_elem)+'_x'+str(one_minus_x_elem)+'_{1-x}}$')
        ax.set_xlabel('$x$')
        ax.set_ylabel('formation enthalpy per atom (eV)')
        plt.show()
        return points, disorder, hull, fig
    
def disorder_hull(doc):
    """ Broaden points on phase diagram by 
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
