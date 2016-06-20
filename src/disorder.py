#!/usr/bin/python
import numpy as np
import string
from copy import deepcopy


def disorder_hull(doc, warren=True):
    """ Broaden structures on phase diagram by
    a measure of local stoichiometry.
    """
    num_atoms = doc['num_atoms']
    if 'lattice_cart' not in doc:
        return np.NaN, warren
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
        return warren_cowley(atoms, disps), warren
    else:
        return bond_disorder(atoms, disps), warren


def disorder_swaps(self, doc, pairs=1, template_param=None):
    """ Take a db document as input and perform atomic swaps. """
    for source in doc['source']:
        if '.castep' or '.res' in source:
            name = source.split('/')[-1].split('.')[0]
    name = name + '-' + str(pairs) + '-pair-swaps/' + name
    swapDoc = deepcopy(doc)
    swapAtoms = swapDoc['atom_types']
    for i in range(pairs):
        valid = False
        while not valid:
            swap = np.random.randint(0, len(swapAtoms) - 1, size=2)
            if swap[0] != swap[1] and swapAtoms[swap[0]] != swapAtoms[swap[1]]:
                    valid = True
        swapAtoms[swap[1]], swapAtoms[swap[0]] = swapAtoms[swap[0]], swapAtoms[swap[1]]
    swapPos = np.asarray(swapDoc['positions_frac'])
    for i in range(len(swapAtoms)):
        swapPos[i] += np.random.rand(3) * (0.1 / 7.9)
    swapDoc['positions_frac'] = swapPos
    hash = self.generate_hash(8)
    self.doc2cell(swapDoc, name + '-' + hash)
    self.doc2param(swapDoc, name + '-' + hash, template_param)


def generate_hash(self, hashLen=6):
    """ Quick hash generator, based on implementation in PyAIRSS by J. Wynn. """
    hashChars = [str(x) for x in range(0, 10)] + [x for x in string.ascii_lowercase]
    hash = ''
    for i in range(hashLen):
        hash += np.random.choice(hashChars)
    return hash
