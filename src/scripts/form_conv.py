#!/usr/bin/python
from __future__ import print_function
from scrapers.castep_scrapers import castep2dict
from os import walk
from sys import exit
from collections import defaultdict
import numpy as np

structure_files = defaultdict(list)
for root, dirs, files in walk('.', topdown=True, followlinks=True):
    for file in files:
        if file.endswith('.castep'):
            castep_dict, success = castep2dict(root + '/' + file)
            if not success:
                exit('Failure to read castep')
            else:
                structure_files[castep_dict['source'][0].split('/')[-1]].append(castep_dict)

cutoff_chempots_dict = defaultdict(dict)
for key in structure_files:
    for doc in structure_files[key]:
        if len(doc['stoichiometry']) == 1:
            doc['formation_enthalpy_per_atom'] = 0
            cutoff_chempots_dict[str(doc['cut_off_energy'])][doc['atom_types'][0]] = doc['enthalpy_per_atom']

cutoff_form = defaultdict(list)
stoich_list = dict()
for key in structure_files:
    for doc in structure_files[key]:
        if len(doc['stoichiometry']) == 2:
            doc['formation_enthalpy_per_atom'] = doc['enthalpy_per_atom']
            for atom in doc['atom_types']:
                doc['formation_enthalpy_per_atom'] -= cutoff_chempots_dict[str(doc['cut_off_energy'])][atom] / len(doc['atom_types'])
            cutoff_form[key].append([doc['cut_off_energy'], doc['formation_enthalpy_per_atom']])
            stoich_list[key] = (doc['stoichiometry'][0][0] + str(doc['stoichiometry'][0][1]) + 
                                doc['stoichiometry'][1][0] + str(doc['stoichiometry'][1][1]))
for key in cutoff_form:
    cutoff_form[key].sort()
    cutoff_form[key] = np.asarray(cutoff_form[key])

import matplotlib.pyplot as plt
plt.style.use('bmh')
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111)
for key in cutoff_form:
    ax.plot(1/cutoff_form[key][:,0], cutoff_form[key][:,1]-cutoff_form[key][-1,1], 'o--', alpha=1, label=stoich_list[key], lw=2)
# ax.plot(1/struct1[:,0], struct1[:,1]-struct1[-1,1], 'o--', alpha=1, label=struct1_stoich, lw=2)
# ax.plot(1/struct2[:,0], struct2[:,1]-struct2[-1,1], 'o--', alpha=0.8, label=struct2_stoich, lw=2)
plt.legend(loc=3)
ax.axhline(0, linestyle='--', c='k', alpha=0.4)
ax.set_xticks(1/cutoff_form[key][:,0])
ax.set_ylabel('$E(e_c) - E(' + str(cutoff_form[key][-1,0]) + ' \mathrm{eV})$')
ax.set_xlabel('$1/e_c  (\mathrm{eV}^{-1})$')
ax.set_xticklabels(['1/' + str(E) for E in cutoff_form[key][:,0]])
plt.show()
