#!/usr/bin/env python
from sys import argv
from matador.utils.cell_utils import get_special_kpoints_for_lattice


def plot_spectral(seed, bandstructure=True, dos=True):
    from matador.scrapers.castep_scrapers import bands2dict
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set(style='whitegrid', font_scale=1.4)
    sns.set_palette('Dark2')
    colours = sns.color_palette()
    valence = colours[0]
    conduction = colours[1]

    plot_window = (-10, 10)

    if bandstructure and not dos:
        fig, ax_bs = plt.subplots(figsize=(5, 5))
    elif dos and not bandstructure:
        fig, ax_dos = plt.subplots(figsize=(3, 7))
    elif bandstructure and dos:
        fig, ax_grid = plt.subplots(1, 2, figsize=(6.5, 5),
                                    gridspec_kw={'width_ratios': [5, 1],
                                                 'wspace': 0.05,
                                                 'left': 0.15})

    ax_bs = ax_grid[0]
    ax_dos = ax_grid[1]

    if bandstructure:
        bs, s = bands2dict(seed + '.bands')
        path = np.linspace(0, 1, int(bs['num_kpoints']))
        for ns in range(bs['num_spins']):
            for nb in range(bs['num_bands']):
                if np.mean(bs['eigenvalues_k_s'][ns][nb]) < 0:
                    colour = valence
                else:
                    colour = conduction
                ax_bs.plot(path, bs['eigenvalues_k_s'][ns][nb], c=colour, lw=1, marker=None)
        ax_bs.axhline(0, ls='--', lw=1, c='grey')
        ax_bs.set_ylim(plot_window)
        ax_bs.set_ylabel('$\epsilon_k$ (eV)')
        ax_bs.set_xlim(0, 1)
        ax_bs.yaxis.grid(False)
        lattice = 'hexagonal'
        special_points = get_special_kpoints_for_lattice(lattice)
        xticks = []
        xticklabels = []
        for ind, kpt in enumerate(bs['kpoint_path']):
            for label, point in special_points.items():
                if np.allclose(point, kpt):
                    if label == 'G':
                        label = '$\Gamma$'
                    else:
                        label = '${}$'.format(label)
                    xticklabels.append(label)
                    xticks.append(path[ind])
                    break

        ax_bs.set_xticks(xticks)
        ax_bs.set_xticklabels(xticklabels)

    if dos:
        dos_data = np.loadtxt(seed + '.adaptive.dat')
        energies = dos_data[:, 0]
        dos = dos_data[:, 1]
        ax_dos.set_ylim(plot_window)
        ax_dos.axvline(0, c='k')
        # ax_dos.set_xlabel('DOS)')
        ax_dos.set_yticks([])
        max_density = np.max(dos[np.where(energies > plot_window[0])])
        ax_dos.set_xticks([])
        ax_dos.set_xlim(0, max_density*1.2)
        ax_dos.axhline(0, c='grey', ls='--', lw=1)
        ax_dos.grid('off')
        ax_dos.plot(dos, energies, lw=1, c='k')
        # ax_dos.plot(dos[np.where(energies < 0)], energies[np.where(energies < 0)], lw=1.5, c=valence)
        # ax_dos.plot(dos[np.where(energies >= 0)], energies[np.where(energies >= 0)], lw=1.5, c='k')
        ax_dos.fill_betweenx(energies[np.where(energies < 0)], 0, dos[np.where(energies < 0)], facecolor=valence, alpha=0.5)
        ax_dos.fill_betweenx(energies[np.where(energies >= 0)], 0, dos[np.where(energies >= 0)], facecolor=conduction, alpha=0.5)

    plt.show()


if __name__ == '__main__':
    plot_spectral(argv[1])
