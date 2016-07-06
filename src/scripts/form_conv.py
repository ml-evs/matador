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
            print(root + '/' + file)
            castep_dict, success = castep2dict(root + '/' + file)
            if not success:
                print('Failure to read castep')
            else:
                source = castep_dict['source'][0].split('/')[-1]
                source = source.replace('.castep', '')
                source = ''.join(source.split('_')[:-1])
                print(source)
                structure_files[source].append(castep_dict)

cutoff_chempots_dict = defaultdict(dict)
cutoff_chempots = defaultdict(list)
chempot_list = dict()
for key in structure_files:
    for doc in structure_files[key]:
        if len(doc['stoichiometry']) == 1:
            doc['formation_enthalpy_per_atom'] = 0
            cutoff_chempots_dict[str(doc['cut_off_energy'])][doc['atom_types'][0]] = doc['enthalpy_per_atom']
            cutoff_chempots[key].append([doc['cut_off_energy'], doc['enthalpy_per_atom']])
            chempot_list[key] = doc['stoichiometry'][0][0] 

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

for key in cutoff_chempots:
    cutoff_chempots[key].sort()
    cutoff_chempots[key] = np.asarray(cutoff_chempots[key])

import matplotlib.pyplot as plt
plt.style.use('bmh')
fig = plt.figure(facecolor='w', figsize=(10,10))
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
for key in cutoff_form:
    ax.plot(1/cutoff_form[key][:,0], cutoff_form[key][:,1]-cutoff_form[key][-1,1], 'o--', alpha=1, label=stoich_list[key], lw=2)
ax.legend(loc=2)
ax.axhline(0, linestyle='--', c='k', alpha=0.4)
ax.set_xticks(1/cutoff_form[key][:,0])
ax.set_ylabel('$E(e_c) - E(' + str(cutoff_form[key][-1,0]) + ' \mathrm{eV})$')
ax.set_xlabel('$1/e_c  (\mathrm{eV}^{-1})$')
ax.set_xticklabels(['1/' + str(E) for E in cutoff_form[key][:,0]])
for key in cutoff_chempots:
    ax2.plot(1/cutoff_chempots[key][:,0], cutoff_chempots[key][:,1]-cutoff_chempots[key][-1, 1], 'o--', alpha=1, label=chempot_list[key], lw=2)
ax2.legend(loc=2)
ax2.axhline(0, linestyle='--', c='k', alpha=0.4)
ax2.set_xticks(1/cutoff_chempots[key][:,0])
ax2.set_ylabel('$E(e_c) - E(' + str(cutoff_chempots[key][-1,0]) + ' \mathrm{eV})$')
ax2.set_xlabel('$1/e_c  (\mathrm{eV}^{-1})$')
ax2.set_xticklabels(['1/' + str(E) for E in cutoff_chempots[key][:,0]])
plt.show()
