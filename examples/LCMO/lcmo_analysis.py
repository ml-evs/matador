from __future__ import print_function
import pymongo as pm
import matplotlib.pyplot as plt
import query as fryan
import numpy as np
from numpy.linalg import norm
from collections import defaultdict
# plt.style.use('bmh')
# plt.style.use('dark_background')
client = pm.MongoClient()
db = client.crystals


def lcmo_info_from_tag(tag, db=['scratch']):
    args = {'db': 'scratch', 'tags': [tag]}
    query = fryan.DBQuery(**args)
    pressure = []
    free_energy = []
    cell_volume = []
    lattice_abc = []
    neighbours = []
    octahedra = []
    octahedra_angles = []
    popn = []
    bonds = []
    charges = []
    spins = []
    net_spin = []
    abs_spin = []
    for doc in query.cursor.sort('pressure', pm.ASCENDING):
        if doc['pressure'] == 'xxx':
            continue
        else:
            pressure.append(doc['pressure'])
            free_energy.append(doc['free_energy_per_atom'])
            cell_volume.append(doc['cell_volume'] / doc['num_fu'])
            lattice_abc.append(doc['lattice_abc'])
            temp_neighbours, temp_oct, temp_ang = b_site_nn(doc)
            neighbours.append(temp_neighbours)
            octahedra.append(temp_oct)
            octahedra_angles.append(temp_ang)
            try:
                b_site_bond, b_site_popn = b_site_bonds(doc)
                popn.append(b_site_popn)
                bonds.append(b_site_bond)
            except:
                pass
            try:
                charges.append(mulliken_charges(doc))
                spins.append(mulliken_spin(doc))
            except:
                pass
            try:
                net_spin.append(doc['mulliken_net_spin'])
                abs_spin.append(doc['mulliken_abs_spin'])
            except:
                net_spin.append(0)
                abs_spin.append(0)
    pressure = np.asarray(pressure)
    free_energy = np.asarray(free_energy)
    cell_volume = np.asarray(cell_volume)
    lattice_abc = np.asarray(lattice_abc)

    unstrained_lattice = lattice_abc[np.argmin(np.abs(pressure))]
    tag_dict = dict()
    if len(pressure) > 2:
        strain = np.zeros((len(pressure), 3))
        for i in range(len(pressure)):
            for j in range(3):
                # strain along <100>, <010> and <001>
                strain[i, j] = (lattice_abc[i][0][j]-unstrained_lattice[0][j]) / \
                    unstrained_lattice[0][j]
        strain *= 100
        tag_dict['strain'] = strain
        tag_dict['poisson_ratio'] = np.zeros((len(pressure)))
        for i in range(len(pressure)):
            tag_dict['poisson_ratio'][i] = tag_dict['strain'][i][2] / \
                (tag_dict['strain'][i][2] - (tag_dict['strain'][i][0]+tag_dict['strain'][i][1]))
    tag_dict['neighbours'] = neighbours
    tag_dict['bonds'] = bonds
    tag_dict['popn'] = popn
    tag_dict['charges'] = charges
    tag_dict['spins'] = spins
    tag_dict['volume'] = cell_volume
    tag_dict['octahedra'] = octahedra
    tag_dict['octahedra_angles'] = octahedra_angles
    tag_dict['net_spin'] = net_spin
    tag_dict['abs_spin'] = abs_spin
    tag_dict['pressure'] = pressure
    tag_dict['free_energy'] = free_energy
    tag_dict['lattice_abc'] = lattice_abc
    return tag_dict

def b_site_nn(doc):
    positions_frac = np.asarray(doc['positions_frac'])
    atom_types = np.asarray(doc['atom_types'])
    lattice_cart = np.asarray(doc['lattice_cart'])
    # print(lattice_cart)
    num_images = 27
    supercell_frac = np.zeros((num_images*len(positions_frac), 3))
    num_atoms = len(positions_frac)
    image = 0
    while image < num_images:
        for a in range(-1,2,1):
            for b in range(-1,2,1):
                for c in range(-1,2,1):
                    supercell_frac[image*num_atoms:(image+1)*num_atoms][:,0] = positions_frac[:, 0] + a
                    supercell_frac[image*num_atoms:(image+1)*num_atoms][:,1] = positions_frac[:, 1] + b
                    supercell_frac[image*num_atoms:(image+1)*num_atoms][:,2] = positions_frac[:, 2] + c
                    image += 1
    positions_abs = np.zeros_like(positions_frac)
    for i in range(len(positions_frac)):
        for j in range(3):
            # print(lattice_cart[j])
            positions_abs[i] += lattice_cart[j] * positions_frac[i, j]
    supercell_abs = np.zeros_like(supercell_frac)
    for i in range(len(supercell_frac)):
        for j in range(3):
            supercell_abs[i] += lattice_cart[j] * supercell_frac[i, j]
    np.savetxt('test2.cell', positions_frac)
    supercell_types = num_images*[atom_types]
    supercell_types = np.asarray(supercell_types).flatten()
    dists = defaultdict(list)
    oct_dict = defaultdict(list)
    ang_dict = defaultdict(list)
    pairs = [['Mn', 'Co', 6], ['Mn', 'Mn', 6], ['Co', 'Co', 6], ['Mn', 'O', 2.2], ['Co', 'O', 2.2], ['La', 'O', 2.5], ['Ti', 'O', 2.5], ['Ti', 'Ti', 4]]
    vec_001 = np.array([0, 0, 1])
    vec_010 = np.array([0, 1, 0])
    vec_100 = np.array([1, 0, 0])
    for pair in pairs:
        elem = pair[0]
        elem2 = pair[1]
        scale = pair[2]
        for i in range(len(positions_abs)):
            if atom_types[i] == elem:
                for j in range(len(supercell_types)):
                    if supercell_types[j] == elem2:
                        temp_dist = 0
                        for k in range(3):
                            temp_dist += (positions_abs[i][k] - supercell_abs[j][k])**2
                        if temp_dist != 0 and np.sqrt(temp_dist) <= scale:
                            dists[pair[0]+pair[1]].append(np.sqrt(temp_dist))
                            if pair[1] == 'O' and pair[0] != 'O':
                                temp_disp = np.zeros((3))
                                for k in range(3):
                                    temp_disp[k] = supercell_abs[j][k] - positions_abs[i][k]
                                oct_dict[elem+'_'+str(i)].append(temp_disp)
                                ang_100 = (180/np.pi)*np.arccos(np.dot(vec_100, temp_disp) / norm(vec_100) / norm(temp_disp))
                                ang_010 = (180/np.pi)*np.arccos(np.dot(vec_010, temp_disp) / norm(vec_010) / norm(temp_disp))
                                ang_001 = (180/np.pi)*np.arccos(np.dot(vec_001, temp_disp) / norm(vec_001) / norm(temp_disp))
                                ang_dict[elem+'_'+str(i)].append([ang_100, ang_010, ang_001])
    return dists, oct_dict, ang_dict

def mulliken_charges(doc):
    charge_dict = defaultdict(list)
    for ind, charge in enumerate(doc['mulliken_charges']):
        charge_dict[doc['atom_types'][ind]].append(charge)
    return charge_dict

def mulliken_spin(doc):
    spin_dict = defaultdict(list)
    for ind, spin in enumerate(doc['mulliken_spins']):
        spin_dict[doc['atom_types'][ind]].append(spin)
    return spin_dict

def b_site_bonds(doc):
    bond_dict = defaultdict(list)
    popn_dict = defaultdict(list)
    positions_frac = np.asarray(doc['positions_frac'])
    atom_types = np.asarray(doc['atom_types'])
    lattice_cart = np.asarray(doc['lattice_cart'])
    bonds = np.asarray(doc['bonds'])
    pairs = [['Mn', 'Co'], ['Mn', 'Mn'], ['Co', 'Co'], ['Mn', 'O'], ['Co', 'O'], ['La', 'O'], ['O', 'O']]
    for pair in pairs:
        for bond in bonds:
            if ((atom_types[bond[0][0]] == pair[0] and atom_types[bond[0][1]] == pair[1]) or
               (atom_types[bond[0][1]] == pair[0] and atom_types[bond[0][0]] == pair[1])):
                bond_dict[pair[0]+pair[1]].append(bond[2])
                popn_dict[pair[0]+pair[1]].append(bond[1])
    return bond_dict, popn_dict

def plot_popn(tag_list, calc_dict, species=['Co','Mn','La','O'], bonds=['CoO', 'MnO', 'LaO']):
    """ Plot Mulliken populations of bonds and species. """
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    num_bonds = np.zeros((len(bonds)), dtype=int)
    num_charges = np.zeros((len(species)), dtype=int)
    markers = ['^', 'v']
    if len(tag_list) == 1:
        lines = ['-']
    else:
        lines = ['-', '--']
    # colours = ['purple', 'blue', 'green', 'red']
    colours = []
    for elem in range(10):
        colours.append(ax._get_lines.prop_cycler.next()['color'])
    tag_colours = ['red', 'blue']
    for bond_ind, bond in enumerate(bonds):
        num_bonds[bond_ind] = len(calc_dict[tag_list[0]]['bonds'][0][bond])
    for elem_ind, elem in enumerate(species):
        num_charges[elem_ind] = len(calc_dict[tag_list[0]]['charges'][0][elem])
    for ind, tag in enumerate(tag_list):
        num_pressures = len(calc_dict[tag]['pressure'])
        bond_popn = []
        charges = []
        for j in range(num_pressures):
            for elem_ind, elem in enumerate(bonds):
                bond_popn.append(np.zeros((num_bonds[elem_ind], num_pressures)))
                for k in range(num_bonds[elem_ind]):
                    bond_popn[elem_ind][k, j] = calc_dict[tag]['popn'][j][elem][k]
            for elem_ind, elem in enumerate(species):
                charges.append(np.zeros((num_charges[elem_ind], num_pressures)))
                for k in range(num_charges[elem_ind]):
                    charges[elem_ind][k, j] = calc_dict[tag]['charges'][j][elem][k]
        for elem_ind, elem in enumerate(bonds):
            for k in range(num_bonds[elem_ind]):
                if k == 0:
                    label = elem + ' $\\langle' + tag + '\\rangle$'
                else:
                    label = None
                ax.plot(calc_dict[tag]['pressure'], bond_popn[elem_ind][k, :], ls=lines[ind], c=colours[elem_ind], lw=0.5, alpha=0.5, marker=markers[ind])   
                ax.plot(np.NaN, np.NaN, label=label, ls=lines[ind], c=colours[elem_ind], lw=0, alpha=1, marker=markers[ind])   
            ax.legend(bbox_to_anchor=(-0.1,1))
            ax.set_xlabel('stress on cell (GPa)')
            ax.set_ylabel('bond population')
    for ind, tag in enumerate(tag_list):
        for elem_ind, elem in enumerate(species):
            for k in range(num_charges[elem_ind]):
                if k == 0:
                    label = elem + ' $\\langle' + tag + '\\rangle$'
                else:
                    label = None
                if len(tag_list) > 1:
                    colour = tag_colours[ind]
                else:
                    colour = colours[elem_ind]
                ax2.plot(calc_dict[tag]['pressure'], charges[elem_ind][k, :], c=colours[elem_ind], ls=lines[ind], lw=1, alpha=0.5, marker=markers[ind])
                ax2.plot(np.NaN, np.NaN, c=colours[elem_ind], ls=lines[ind], label=label, lw=1, alpha=1, marker=markers[ind])
                ax2.set_xlabel('stress on cell (GPa)')
                ax2.legend(loc=1, bbox_to_anchor=(1.34,1.05))
                ax2.set_ylabel('Mulliken charge on atom')


def plot_atom_dist(tag_list, calc_dict, label_list=None, elems=['MnMn'], shift=1, lw=4, xaxis=None, figsize=6):
    """ Plot atomic distances. """
    num_plots = len(elems)
    elem_list = [elem for elem in elems]
    fig = plt.figure(figsize=(figsize*num_plots, figsize))
    ax_list = []
    for i in range(num_plots):
        ax_list.append(fig.add_subplot(1,num_plots,i+1))
    if label_list is None:
        label_list = tag_list
    for i in range(num_plots):
        elems = elem_list[i]
        ax = ax_list[i]
        markers = 10*['^', 'v']
        if len(elems) == 1:
            lines = 10*['-']
        else:
            lines = ['-', '--', '-.', ':']
        tag_colours = []
        if '111' in tag_list[0]:
            col = ax._get_lines.prop_cycler.next()['color']
            col = ax._get_lines.prop_cycler.next()['color']
            col = ax._get_lines.prop_cycler.next()['color']
        for ind, tag in enumerate(tag_list):
            if ind == 2:
                col = ax._get_lines.prop_cycler.next()['color']
            tag_colours.append(ax._get_lines.prop_cycler.next()['color'])
        for ind, tag in enumerate(tag_list):
            label = label_list[ind]
            num_pressures = len(calc_dict[tag]['pressure'])
            num_bonds = len(calc_dict[tag]['neighbours'][0][elems])
            atom_dist = np.zeros((num_bonds, num_pressures))
            for j in range(num_pressures):
                for k in range(num_bonds):
                    atom_dist[k][j] = calc_dict[tag]['neighbours'][j][elems][k]
            # sort atom dists by lowest pressure dists
            for j in range(0,num_pressures-1):
                atom_dist[:,j] = atom_dist[:, j][np.argsort(atom_dist[:, -1])]
            atom_dist[:, -1] = np.sort(atom_dist[:, -1])
            if xaxis is None:
                plot_xaxis = calc_dict[tag]['pressure']
                ax.set_xlabel('stress on cell (GPa)')
            else:
                plot_xaxis = xaxis
                ax.set_xlabel('\% transverse strain')
            for k in range(len(atom_dist)):
                if k == 0:
                    label = label
                else:
                    label = None
                ax_list[i].plot(plot_xaxis, atom_dist[k, :], c=tag_colours[ind], ls = '-', lw=lw,alpha=0.2, marker=None, zorder=0, solid_capstyle='round')
                ax_list[i].scatter(plot_xaxis, atom_dist[k, :], c=tag_colours[ind] if len(tag_list)!=1 else colours[tag][k], lw=0.5, alpha=1, marker=markers[ind], zorder=1000, s=20, edgecolors='k')
                ax_list[i].plot(np.NaN, np.NaN, c=tag_colours[ind] if len(tag_list)!=1 else colours[tag][k], label=label, ls='-', lw=5, marker=markers[ind], alpha=0.5, solid_capstyle='round')
            if elems == 'MnO' or elems == 'CoO':
                ax_list[i].set_ylim(1.88, 2.04)
            if i == 0:
                ax_list[i].set_ylabel('neighbour distances (A)')
            else:
                ax_list[i].yaxis.tick_right()
            ax_list[i].set_ylim(ax_list[0].get_ylim())
            ax_list[i].grid('off')
            if i == num_plots-1:
                ax_list[i].legend(loc='lower center',fontsize=12)#, bbox_to_anchor=(1.38,1))

    plt.tight_layout()
    plt.savefig('octo' + tag + '.pdf', dpi=300)


def plot_free_energy(tag_list, calc_dict, label_list=None, xaxis=None, figsize=8):
    fig = plt.figure(figsize=(figsize, figsize))
    ax = fig.add_subplot(111)
    if label_list is None:
        label_list = tag_list
    colours = []
    min_energy = 0
    markers = 10*['^', 'v']
    for ind, tag in enumerate(tag_list):
        if ind == 2:
            col = ax._get_lines.prop_cycler.next()['color']
        colours.append(ax._get_lines.prop_cycler.next()['color'])
        for energy in calc_dict[tag]['free_energy']:
            if energy < min_energy:
                min_energy = energy
    for ind, tag in enumerate(tag_list):
        if xaxis is None:
            plot_xaxis = calc_dict[tag]['pressure']
            ax.set_xlabel('pressure on cell (GPa)')
        else:
            plot_xaxis = xaxis
            ax.set_xlabel('\% transverse strain')
            ax.set_xlim(-3, 3)
        strain_space = np.linspace(-5, 5, 100)
        ax.plot(strain_space,
                np.poly1d(np.polyfit(
                    plot_xaxis,
                    calc_dict[tag]['free_energy'] - min_energy, deg=2))(strain_space),
                lw=2, c=colours[ind], markersize=0)
        ax.scatter(plot_xaxis, calc_dict[tag]['free_energy']-min_energy,
                   c=colours[ind], s=30, lw=0.5, edgecolors='k', zorder=1000, marker=markers[ind])
        ax.plot(np.NaN, np.NaN, label=label_list[ind],
                lw=2, c=colours[ind])
    ax.set_ylabel('Relative free energy (eV/atom)')
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles, labels, bbox_to_anchor=(1.7, 0.8), shadow=True, fontsize=12, ncol=1)
    ax.grid('off')
    plt.tight_layout()
    plt.savefig('energy.pdf', dpi=200)
