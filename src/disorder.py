#!/usr/bin/python
import numpy as np


def disorder_hull(doc, warren=True):
    """ Broaden structures on phase diagram by
    a measure of local stoichiometry.
    """
    num_atoms = doc['num_atoms']
    if 'lattice_cart' not in doc:
        return np.NaN
    lat_cart = doc['lattice_cart']
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
                    disps[i, jindex] += real_disp[k]**2
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

    if warren:
        return warren_cowley(atoms, disps)
    else:
        return bond_disorder(atoms, disps)
