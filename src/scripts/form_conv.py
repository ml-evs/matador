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
            castep_dict, success = castep2dict(root + '/' + file, db=False)
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
            doc['formation_0K_energy_per_atom'] = 0
            cutoff_chempots_dict[str(doc['cut_off_energy'])][doc['atom_types'][0]] = doc['0K_energy_per_atom']
            cutoff_chempots[key].append([doc['cut_off_energy'], doc['0K_energy_per_atom']])
            chempot_list[key] = doc['stoichiometry'][0][0] 

cutoff_form = defaultdict(list)
stoich_list = dict()
for key in structure_files:
    for doc in structure_files[key]:
        if len(doc['stoichiometry']) == 2:
            doc['formation_0K_energy_per_atom'] = doc['0K_energy_per_atom']
            for atom in doc['atom_types']:
                doc['formation_0K_energy_per_atom'] -= cutoff_chempots_dict[str(doc['cut_off_energy'])][atom] / len(doc['atom_types'])
            cutoff_form[key].append([doc['cut_off_energy'], doc['formation_0K_energy_per_atom']])
            stoich_list[key] = (doc['stoichiometry'][1][0])
            stoich_list[key] += str(doc['stoichiometry'][1][1]) if doc['stoichiometry'][1][1] != 1 else ''
            stoich_list[key] += (doc['stoichiometry'][0][0])
            stoich_list[key] += str(doc['stoichiometry'][0][1]) if doc['stoichiometry'][0][1] != 1 else ''
for key in cutoff_form:
    cutoff_form[key].sort()
    cutoff_form[key] = np.asarray(cutoff_form[key])

for key in cutoff_chempots:
    cutoff_chempots[key].sort()
    cutoff_chempots[key] = np.asarray(cutoff_chempots[key])

import matplotlib.pyplot as plt
plt.style.use('bmh')
plt.style.use('article')
fig = plt.figure(facecolor='w', figsize=(5,3))
ax = fig.add_subplot(111)
subax = plt.axes([.6, .5, .25, .3])
# ax2 = fig.add_subplot(122)
for key in cutoff_form:
    # ax.plot(1/cutoff_form[key][:,0], cutoff_form[key][:,1]-cutoff_form[key][-1,1], 'o--', alpha=1, label=stoich_list[key], lw=2)
    ax.plot(cutoff_form[key][:,0], np.abs(cutoff_form[key][:,1]-cutoff_form[key][-1, 1]), 'o-', alpha=1, label=stoich_list[key], lw=2, zorder=1000)
    subax.plot(cutoff_form[key][:,0], np.abs(cutoff_form[key][:,1]-cutoff_form[key][-1, 1]), 'o-', alpha=1, label=stoich_list[key], lw=2, zorder=1000)
ax.set_ylabel('$|\mathrm{E(G_{max}) - E(' + str(cutoff_form[key][-1,0]) + ' eV)}|$')
for key in cutoff_chempots:
    ax.plot(cutoff_chempots[key][:,0], np.abs(cutoff_chempots[key][:,1]-cutoff_chempots[key][-1, 1]), 'o-', alpha=1, label=chempot_list[key], lw=2)
    subax.plot(cutoff_chempots[key][:,0], np.abs(cutoff_chempots[key][:,1]-cutoff_chempots[key][-1, 1]), 'o-', alpha=1, label=chempot_list[key], lw=2)
ax.set_xlabel('$\mathrm{G_{max}} (\mathrm{eV}^{-1}$)')
subax.set_ylim(-0.0001, 0.001)
subax.set_xlim(500, 850)
subax.set_xticks([500, 600, 700, 800])
subax.set_yticks([0, 0.0005, 0.001])
subax.grid('off')
ax.set_xlim(200, 900)
legend = ax.legend(loc='lower center', fontsize=12, ncol=4, shadow=True, bbox_to_anchor=(0.5, 0.01))
# legend.get_frame().set_facecolor('w')
ax.set_ylim(-0.01, 0.03)
ax.grid('off')
# plt.show()
plt.savefig('/home/matthew/proj/thesis/img/lias_conv.pdf', bbox_inches='tight')
