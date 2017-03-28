#!/usr/bin/env python
from __future__ import print_function
from matador.scrapers.castep_scrapers import castep2dict
from os import walk, chdir
from os.path import isdir
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')


def get_files(path):
    chdir(path)
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
                    # print(source)
                    structure_files[source].append(castep_dict)
    chdir('..')
    return structure_files


def get_cutoffs(structure_files):

    cutoff_chempots_dict = defaultdict(dict)
    cutoff_chempots = defaultdict(list)
    chempot_list = dict()
    for key in structure_files:
        for doc in structure_files[key]:
            if len(doc['stoichiometry']) == 1:
                doc['formation_energy_per_atom'] = 0
                cutoff_chempots_dict[str(doc['cut_off_energy'])][doc['atom_types'][0]] = doc['enthalpy_per_atom']
                cutoff_chempots[key].append([doc['cut_off_energy'], doc['enthalpy_per_atom']])
                chempot_list[key] = doc['stoichiometry'][0][0]
    cutoff_form = defaultdict(list)
    stoich_list = dict()
    # print(cutoff_chempots_dict)
    for key in structure_files:
        for doc in structure_files[key]:
            # print(doc['cut_off_energy'])
            if len(doc['stoichiometry']) == 2:
                doc['formation_energy_per_atom'] = doc['enthalpy_per_atom']
                for atom in doc['atom_types']:
                    doc['formation_energy_per_atom'] -= cutoff_chempots_dict[str(doc['cut_off_energy'])][atom] / len(doc['atom_types'])
                cutoff_form[key].append([doc['cut_off_energy'], doc['formation_energy_per_atom']])
                stoich_list[key] = (doc['stoichiometry'][1][0])
                stoich_list[key] += ('$_\mathrm{' + str(doc['stoichiometry'][1][1]) + '}$') if doc['stoichiometry'][1][1] != 1 else ''
                stoich_list[key] += (doc['stoichiometry'][0][0])
                stoich_list[key] += ('$_\mathrm{' + str(doc['stoichiometry'][0][1]) + '}$') if doc['stoichiometry'][0][1] != 1 else ''
    for key in cutoff_form:
        cutoff_form[key].sort()
        cutoff_form[key] = np.asarray(cutoff_form[key])

    for key in cutoff_chempots:
        cutoff_chempots[key].sort()
        cutoff_chempots[key] = np.asarray(cutoff_chempots[key])

    print(cutoff_chempots)

    return cutoff_chempots, cutoff_form, stoich_list, chempot_list


def get_kpts(structure_files):

    kpt_chempots_dict = defaultdict(dict)
    kpt_chempots = defaultdict(list)
    chempot_list = dict()
    for key in structure_files:
        for doc in structure_files[key]:
            if len(doc['stoichiometry']) == 1:
                doc['formation_energy_per_atom'] = 0
                # grid = doc['kpoints_mp_grid']
                # num_k = sum(grid)
                num_k = round(doc['kpoints_mp_spacing'], 2)
                print(doc['atom_types'][0], num_k, doc['kpoints_mp_spacing'])
                kpt_chempots_dict[str(num_k)][doc['atom_types'][0]] = doc['enthalpy_per_atom']
                kpt_chempots[key].append([num_k, doc['enthalpy_per_atom']])
                chempot_list[key] = doc['stoichiometry'][0][0]

    # print(kpt_chempots_dict)
    # print(kpt_chempots)

    kpt_form = defaultdict(list)
    stoich_list = dict()
    for key in structure_files:
        for doc in structure_files[key]:
            print(doc['source'])
            grid = doc['kpoints_mp_grid']
            num_k = sum(grid)
            num_k = doc['kpoints_mp_spacing']
            if len(doc['stoichiometry']) == 2:
                doc['formation_energy_per_atom'] = doc['enthalpy_per_atom']
                print('!', round(doc['kpoints_mp_spacing'], 2), doc['kpoints_mp_spacing'])
                for atom in doc['atom_types']:
                    doc['formation_energy_per_atom'] -= kpt_chempots_dict[str(round(doc['kpoints_mp_spacing'], 2))][atom] / len(doc['atom_types'])
                kpt_form[key].append([num_k, doc['formation_energy_per_atom']])
                stoich_list[key] = (doc['stoichiometry'][1][0])
                stoich_list[key] += ('$_\mathrm{' + str(doc['stoichiometry'][1][1]) + '}$') if doc['stoichiometry'][1][1] != 1 else ''
                stoich_list[key] += (doc['stoichiometry'][0][0])
                stoich_list[key] += ('$_\mathrm{' + str(doc['stoichiometry'][0][1]) + '}$') if doc['stoichiometry'][0][1] != 1 else ''
    for key in kpt_form:
        kpt_form[key].sort()
        kpt_form[key] = np.asarray(kpt_form[key])

    for key in kpt_chempots:
        kpt_chempots[key].sort()
        kpt_chempots[key] = np.asarray(kpt_chempots[key])

    return kpt_chempots, kpt_form, stoich_list, chempot_list


def plot_both(cutoff_chempots, kpt_chempots,
              cutoff_form, kpt_form,
              cutoff_stoich_list, kpt_stoich_list,
              cutoff_chempot_list, kpt_chempot_list):
    import matplotlib.pyplot as plt
    plt.style.use('bmh')
    try:
        plt.style.use('article')
    except:
        pass
    fig = plt.figure(facecolor='w', figsize=(10, 6))
    ax = fig.add_subplot(121, axisbg='w')
    ax2 = fig.add_subplot(122, axisbg='w')
    for ind, key in enumerate(cutoff_form):
        ax.plot(-1/cutoff_form[key][:, 0], np.abs(cutoff_form[key][:, 1]-cutoff_form[key][0, 1])*1000,
                'o-', markersize=5, alpha=1, label=cutoff_stoich_list[key], lw=1, zorder=1000)
    # for key in cutoff_chempots:
        # ax_chempots.plot(-1/cutoff_chempots[key][:, 0], np.abs(cutoff_chempots[key][:, 1]-cutoff_chempots[key][-1, 1])*1000, 'o-', markersize=5, alpha=1, label=cutoff_chempot_list[key], lw=1)
    ax.set_ylabel('Relative energy difference (meV/atom)')
    ax.set_xlabel('plane wave cutoff (eV)')
    cutoffs = np.loadtxt('cutoff.conv')
    ax.set_xticks(-1/cutoffs)
    ax.set_xticklabels(cutoffs)
    ax.legend(loc='upper center', fontsize=10, ncol=4, shadow=True, bbox_to_anchor=(1.0, 1.25))
    # ax.set_ylim(-0.002e3, 0.03e3)
    # ax.set_yticklabels(ax.get_yticks())
    ax.grid('off')

    for key in kpt_form:
        ax2.plot(kpt_form[key][:, 0], np.abs(kpt_form[key][:, 1]-kpt_form[key][0, 1])*1000, 'o-', markersize=5, alpha=1, label=kpt_stoich_list[key], lw=1, zorder=1000)
    for key in kpt_chempots:
        ax2.plot(kpt_chempots[key][:, 0], np.abs(kpt_chempots[key][:, 1]-kpt_chempots[key][0, 1])*1000, 'o-', markersize=5, alpha=1, label=kpt_chempot_list[key], lw=1)
    kpts = list(reversed(np.loadtxt('kpt.conv').tolist()))
    ax2.set_xticks(kpts)
    ax2.set_xlabel('max k-point spacing (1/A)')
    ax2.grid('off')

    # subax = plt.axes([.3, .50, .16, .36], axisbg='w')
    # for key in cutoff_form:
        # subax.plot(1/cutoff_form[key][:, 0], np.abs(cutoff_form[key][:, 1]-cutoff_form[key][-1, 1])*1000,
                   # 'o-', markersize=5, alpha=1, label=cutoff_stoich_list[key], lw=1, zorder=1000)
    # for key in cutoff_chempots:
        # subax.plot(1/cutoff_chempots[key][:, 0], np.abs(cutoff_chempots[key][:, 1]-cutoff_chempots[key][-1, 1])*1000, 'o-', markersize=5, alpha=1, label=cutoff_chempot_list[key], lw=1)
    # subax.set_ylim(-0.0001e3, 0.0005e3)
    # subax.set_xlim(500, 800)
    # subax.set_xticks([500, 600, 700, 800])
    # subax.set_yticks([0, 0.00025e3, 0.0005e3])
    # subax.set_xticklabels(subax.get_xticks())
    # subax.set_yticklabels(subax.get_yticks())
    # subax.grid('off')
    # plt.show()
    # plt.tight_layout()
    plt.savefig('conv.png', bbox_inches='tight')


if __name__ == '__main__':
    if isdir('completed_cutoff'):
        print('Finding cutoffs...')
        cutoff_structure_files = get_files('completed_cutoff')
        cutoff_chempots, cutoff_form, cutoff_stoich_list, cutoff_chempot_list = get_cutoffs(cutoff_structure_files)
    if isdir('completed_kpts'):
        print('Finding kpts...')
        kpt_structure_files = get_files('completed_kpts')
        kpt_chempots, kpt_form, kpt_stoich_list, kpt_chempot_list = get_kpts(kpt_structure_files)
    plot_both(cutoff_chempots, kpt_chempots,
              cutoff_form, kpt_form,
              cutoff_stoich_list, kpt_stoich_list,
              cutoff_chempot_list, kpt_chempot_list)
