""" This file implements several useful plotting routines. """


def plot_spectral(seeds, plot_bandstructure=True, plot_dos=False, cell=False, gap=False, **kwargs):
    """ Plot bandstructure and optional DOS from <seed>.bands and
    <seed>.adaptive.dat file.

    Input:

        | seeds: list, list of filenames of bands files

    Args:

        | plot_bandstructure : bool, whether to plot bandstructure
        | plot_dos           : bool, whether to plot density of states
        | cell          : bool, whether to work out correct labels from structure in cell file

    """
    from matador.scrapers.castep_scrapers import bands2dict, cell2dict
    from matador.utils.cell_utils import doc2spg
    from seekpath import get_path
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set(style='whitegrid', font_scale=1.4)
    sns.set_palette('Dark2')
    colours = sns.color_palette()
    valence = colours[0]
    conduction = colours[1]
    crossing = colours[2]

    if not isinstance(seeds, list):
        seeds = [seeds]
    print(seeds)
    ls = []
    for i in range(len(seeds)):
        if i % 3 == 0:
            ls.append('-')
        elif i % 3 == 1:
            ls.append('--')
        elif i % 3 == 2:
            ls.append('-.')

    if kwargs.get('plot_window') is not None:
        plot_window = (-kwargs.get('plot_window'), kwargs.get('plot_window'))
    else:
        plot_window = (-5, 5)

    if plot_bandstructure and not plot_dos:
        fig, ax_bs = plt.subplots(figsize=(5, 5))
    elif plot_bandstructure and plot_dos:
        fig, ax_grid = plt.subplots(1, 2, figsize=(6.5, 5), sharey=True,
                                    gridspec_kw={'width_ratios': [5, 1],
                                                 'wspace': 0.05,
                                                 'left': 0.15})
        ax_bs = ax_grid[0]
        ax_dos = ax_grid[1]
    elif not plot_bandstructure and plot_dos:
        fig, ax_dos = plt.subplots(figsize=(7, 5))

    for seed_ind, seed in enumerate(seeds):
        seed = seed.replace('.bands', '')
        if plot_bandstructure:
            bs, s = bands2dict(seed + '.bands', summary=True)
            path = np.linspace(0, 1, int(bs['num_kpoints'])-len(bs['kpoint_branches'])+1)
            for branch_ind, branch in enumerate(bs['kpoint_branches']):
                for ns in range(bs['num_spins']):
                    for nb in range(bs['num_bands']):
                        if np.max(bs['eigenvalues_k_s'][ns][nb][branch]) < 0:
                            colour = valence
                        elif np.min(bs['eigenvalues_k_s'][ns][nb][branch]) > 0:
                            colour = conduction
                        elif np.min(bs['eigenvalues_k_s'][ns][nb][branch]) < 0 and np.max(bs['eigenvalues_k_s'][ns][nb][branch]) > 0:
                            colour = crossing
                        else:
                            colour = 'black'
                        ax_bs.plot(path[(np.asarray(branch)-branch_ind).tolist()], bs['eigenvalues_k_s'][ns][nb][branch], c=colour, lw=1, marker=None, ls=ls[seed_ind])
                if branch_ind != len(bs['kpoint_branches'])-1:
                    ax_bs.axvline(path[branch[-1]-branch_ind], ls='-.', lw=1, c='grey')
            ax_bs.axhline(0, ls='--', lw=1, c='grey')
            ax_bs.set_ylim(plot_window)
            ax_bs.set_ylabel('$\epsilon_k$ (eV)')
            ax_bs.set_xlim(0, 1)
            ax_bs.yaxis.grid(False)
            xticks = []
            xticklabels = []
            shear_planes = []
            labelled = []

            if cell:
                doc, success = cell2dict(seed + '.cell', db=False, verbosity=10, outcell=True, positions=True)
                if success:
                    spg_structure = doc2spg(doc)
                    if spg_structure is not False:
                        seekpath_results = get_path(spg_structure)
                    path_labels = seekpath_results['point_coords']

                for branch_ind, branch in enumerate(bs['kpoint_branches']):
                    for sub_ind, ind in enumerate(branch):
                        kpt = bs['kpoint_path'][ind]
                        for label, point in path_labels.items():
                            if np.allclose(point, kpt):
                                if ind-branch_ind not in labelled:
                                    label = label.replace('GAMMA', '\Gamma')
                                    label = label.replace('SIGMA', '\Sigma')
                                    if sub_ind == len(branch)-1:
                                        if branch_ind < len(bs['kpoint_branches'])-1:
                                            next_point = bs['kpoint_path'][bs['kpoint_branches'][branch_ind+1][0]]
                                            for new_label, new_point in path_labels.items():
                                                new_label = new_label.replace('GAMMA', '\Gamma')
                                                new_label = new_label.replace('SIGMA', '\Sigma')
                                                if np.allclose(new_point, next_point):
                                                    label = '{}|{}'.format(label, new_label)
                                                    labelled.append(ind-branch_ind)
                                                    shear_planes.append(ind-branch_ind)
                                    label = '${}$'.format(label)
                                    xticklabels.append(label)
                                    xticks.append(path[ind-branch_ind])
                                    break
            else:
                if np.allclose(next_point, [0, 0, 0]):
                    label = '$\Gamma$'
                    xticklabels.append(label)
                    xticks.append(path[ind-branch_ind])
            if gap and bs['band_gap'] > 0:
                vbm_pos = bs['band_gap_path_inds'][1]
                vbm = bs['valence_band_min']
                cbm_pos = bs['band_gap_path_inds'][0]
                cbm = bs['conduction_band_max']
                vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                ax_bs.plot([path[vbm_pos-vbm_offset], path[cbm_pos-cbm_offset]], [vbm, cbm], c='red')

            ax_bs.set_xticks(xticks)
            ax_bs.set_xticklabels(xticklabels)
        if plot_dos:
            dos_data = np.loadtxt(seed + '.adaptive.dat')
            energies = dos_data[:, 0]
            dos = dos_data[:, 1]
            max_density = np.max(dos[np.where(energies > plot_window[0])])
            if plot_bandstructure:
                ax_dos.set_xticks([0.6*max_density])
                ax_dos.set_xticklabels(['DOS'])
                ax_dos.axhline(0, c='grey', ls='--', lw=1)
                ax_dos.set_xlim(0, max_density*1.2)
                ax_dos.set_ylim(plot_window)
                ax_dos.axvline(0, c='k')
                ax_dos.set_yticks([])
                ax_dos.plot(dos, energies, lw=1, c='k', ls=ls[seed_ind])
                if seed_ind == 0:
                    ax_dos.fill_betweenx(energies[np.where(energies <= 0)], 0, dos[np.where(energies <= 0)], facecolor=valence, alpha=0.5)
                    ax_dos.fill_betweenx(energies[np.where(energies >= 0)], 0, dos[np.where(energies >= 0)], facecolor=conduction, alpha=0.5)
            else:
                ax_dos.set_ylabel('DOS')
                ax_dos.set_xlabel('Energy (eV)')
                ax_dos.axvline(0, c='grey', ls='--', lw=1)
                ax_dos.set_ylim(0, max_density*1.2)
                ax_dos.set_xlim(plot_window)
                ax_dos.axhline(0, c='k')
                ax_dos.plot(energies, dos, lw=1, c='k', ls=ls[seed_ind])
                if seed_ind == 0:
                    ax_dos.fill_between(energies[np.where(energies <= 0)], 0, dos[np.where(energies <= 0)], facecolor=valence, alpha=0.5)
                    ax_dos.fill_between(energies[np.where(energies >= 0)], 0, dos[np.where(energies >= 0)], facecolor=conduction, alpha=0.5)
            ax_dos.grid('off')

    if kwargs.get('pdf'):
        plt.savefig(seed.replace('.bands', '') + '_spectral.pdf', bbox_inches='tight')
    if kwargs.get('png'):
        plt.savefig(seed.replace('.bands', '') + '_spectral.png', bbox_inches='tight', dpi=300)

    plt.show()
