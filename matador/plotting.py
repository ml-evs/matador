""" This file implements several useful plotting routines. """


def plot_spectral(seeds,
                  plot_bandstructure=True, plot_dos=False, phonons=False,
                  cell=False, gap=False, colour_by_seed=False, verbosity=0, **kwargs):
    """ Plot bandstructure and optional DOS from <seed>.bands and
    <seed>.adaptive.dat file.

    Input:

        | seeds: list, list of filenames of bands files

    Args:

        | plot_bandstructure : bool, whether to plot bandstructure
        | plot_dos           : bool, whether to plot density of states
        | phonons            : bool, whether to plot phonon or electronic data
        | cell               : bool, whether to work out correct labels from structure in cell file
        | gap                : bool, draw on the band gap
        | colour_by_seed     : bool, plot with a separate colour per bandstructure

    """
    from matador.scrapers.castep_scrapers import bands2dict, cell2dict, phonon2dict
    from matador.utils.cell_utils import doc2spg
    from seekpath import get_path
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from os.path import isfile
    sns.set(style='whitegrid', font_scale=1.2)
    sns.set_style({
        'axes.facecolor': 'white', 'figure.facecolor': 'white',
        'font.sans-serif': ['Linux Biolinum O', 'Helvetica', 'Arial'],
        'axes.linewidth': 0.5,
        'axes.grid': False,
        'legend.frameon': False,
        'axes.axisbelow': True})
    sns.set_palette('Dark2')
    colours = sns.color_palette()
    valence = colours[0]
    conduction = colours[1]
    crossing = colours[2]

    if not isinstance(seeds, list):
        seeds = [seeds]

    if len(seeds) > 1 and colour_by_seed:
        seed_colours = colours
        ls = ['-']*len(seeds)
    else:
        ls = []
        for i in range(len(seeds)):
            if i % 3 == 0:
                ls.append('-')
            elif i % 3 == 1:
                ls.append('--')
            elif i % 3 == 2:
                ls.append('-.')
    if len(seeds) == 1:
        colour_by_seed = False

    if kwargs.get('plot_window') is not None:
        plot_window = (-kwargs.get('plot_window'), kwargs.get('plot_window'))
    else:
        plot_window = (-5, 5)

    if plot_bandstructure and not plot_dos:
        fig, ax_dispersion = plt.subplots(figsize=(5, 5))
    elif plot_bandstructure and plot_dos:
        fig, ax_grid = plt.subplots(1, 2, figsize=(6.5, 5), sharey=True,
                                    gridspec_kw={'width_ratios': [5, 1],
                                                 'wspace': 0.05,
                                                 'left': 0.15})
        ax_dispersion = ax_grid[0]
        ax_dos = ax_grid[1]
    elif not plot_bandstructure and plot_dos:
        fig, ax_dos = plt.subplots(figsize=(7, 5))

    for seed_ind, seed in enumerate(seeds):
        seed = seed.replace('.bands', '').replace('.phonon', '')
        if plot_bandstructure:
            if phonons:
                dispersion, s = phonon2dict(seed + '.phonon', verbosity=verbosity)
                branch_key = 'qpoint_branches'
                num_key = 'num_qpoints'
                path_key = 'qpoint_path'
                eig_key = 'eigenvalues_q'
                band_key = 'num_branches'
                dispersion['num_spins'] = 1
                spin_key = 'num_spins'
                plot_window = [np.min(dispersion[eig_key]), np.max(dispersion[eig_key])]
            else:
                dispersion, s = bands2dict(seed + '.bands', summary=True, gap=gap, verbosity=verbosity)
                branch_key = 'kpoint_branches'
                num_key = 'num_kpoints'
                path_key = 'kpoint_path'
                eig_key = 'eigenvalues_k_s'
                band_key = 'num_bands'
                spin_key = 'num_spins'
            path = [0]
            for branch in dispersion[branch_key]:
                for ind, kpt in enumerate(dispersion[path_key][branch]):
                    if ind != len(branch) - 1:
                        diff = np.sqrt(np.sum((kpt - dispersion[path_key][branch[ind+1]])**2))
                        path.append(path[-1] + diff)
            path = np.asarray(path)
            path /= np.max(path)
            assert len(path) == int(dispersion[num_key]) - len(dispersion[branch_key]) + 1
            if phonons:
                dispersion[eig_key] = modes2bands(dispersion[eig_key], path, dispersion[branch_key])
            for branch_ind, branch in enumerate(dispersion[branch_key]):
                for ns in range(dispersion[spin_key]):
                    for nb in range(dispersion[band_key]):
                        if not phonons:
                            if np.max(dispersion[eig_key][ns][nb][branch]) < 0:
                                if colour_by_seed:
                                    colour = seed_colours[seed_ind]
                                else:
                                    colour = valence
                            elif np.min(dispersion[eig_key][ns][nb][branch]) > 0:
                                if colour_by_seed:
                                    colour = seed_colours[seed_ind]
                                else:
                                    colour = conduction
                            elif np.min(dispersion[eig_key][ns][nb][branch]) < 0 and np.max(dispersion[eig_key][ns][nb][branch]) > 0:
                                if colour_by_seed:
                                    colour = seed_colours[seed_ind]
                                else:
                                    colour = crossing
                            else:
                                if colour_by_seed:
                                    colour = seed_colours[seed_ind]
                                else:
                                    colour = 'black'
                        else:
                            colour = 'black'
                        ax_dispersion.plot(path[(np.asarray(branch)-branch_ind).tolist()], dispersion[eig_key][ns][nb][branch], c=colour, lw=1, marker=None, ls=ls[seed_ind], alpha=1)
                        # ax_dispersion.scatter(path[(np.asarray(branch)-branch_ind).tolist()], dispersion[eig_key][ns][nb][branch], c=colour, s=5, marker='o')
                if branch_ind != len(dispersion[branch_key])-1:
                    ax_dispersion.axvline(path[branch[-1]-branch_ind], ls='-.', lw=1, c='grey')
            ax_dispersion.axhline(0, ls='--', lw=1, c='grey')
            ax_dispersion.set_ylim(plot_window)
            if phonons:
                ylabel = 'Wavenumber (cm$^{-1}$)'
            else:
                ylabel = '$\epsilon_k$ (eV)'
            ax_dispersion.set_ylabel(ylabel)
            ax_dispersion.set_xlim(0, 1)
            if phonons:
                dispersion['freq_unit'] = dispersion['freq_unit'].replace('-1', '$^{-1}$')
                ax_dispersion.axhline(np.min(dispersion['softest_mode_freq']), ls=ls[seed_ind], c='r',
                                      label='$\omega_\mathrm{{min}} = {:5.3f}$ {}'.format(np.min(dispersion['softest_mode_freq']), dispersion['freq_unit']))
                ax_dispersion.legend()
            xticks = []
            xticklabels = []
            shear_planes = []
            labelled = []

            if cell:
                doc, success = cell2dict(seed + '.cell', db=False, verbosity=verbosity, outcell=True, positions=True)
                if success:
                    spg_structure = doc2spg(doc)
                    if spg_structure is not False:
                        seekpath_results = get_path(spg_structure)
                    path_labels = seekpath_results['point_coords']

                for branch_ind, branch in enumerate(dispersion[branch_key]):
                    for sub_ind, ind in enumerate(branch):
                        kpt = dispersion[path_key][ind]
                        for label, point in path_labels.items():
                            if np.allclose(point, kpt):
                                if ind-branch_ind not in labelled:
                                    label = label.replace('GAMMA', '\Gamma')
                                    label = label.replace('SIGMA', '\Sigma')
                                    if sub_ind == len(branch)-1:
                                        if branch_ind < len(dispersion[branch_key])-1:
                                            next_point = dispersion[path_key][dispersion[branch_key][branch_ind+1][0]]
                                            for new_label, new_point in path_labels.items():
                                                new_label = new_label.replace('GAMMA', '\Gamma')
                                                new_label = new_label.replace('SIGMA', '\Sigma')
                                                if np.allclose(new_point, next_point):
                                                    label = '{}|{}'.format(label, new_label)
                                                    labelled.append(ind-branch_ind)
                                                    shear_planes.append(ind)
                                    label = '${}$'.format(label)
                                    xticklabels.append(label)
                                    xticks.append(path[ind-branch_ind])
                                    break
            # else:
                # for point in dispersion[branch_key]:
                    # if np.allclose(next_point, [0, 0, 0]):
                        # label = '$\Gamma$'
                        # xticklabels.append(label)
                        # xticks.append(path[ind-branch_ind])
            if not phonons and gap and dispersion['band_gap'] > 0:
                vbm_pos = dispersion['band_gap_path_inds'][1]
                vbm = dispersion['valence_band_min']
                cbm_pos = dispersion['band_gap_path_inds'][0]
                cbm = dispersion['conduction_band_max']
                vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                ax_dispersion.plot([path[vbm_pos-vbm_offset], path[cbm_pos-cbm_offset]], [vbm, cbm], ls=ls[seed_ind], c='blue', label='indirect gap {:3.3f} eV'.format(cbm-vbm))
                if cbm_pos != vbm_pos:
                    vbm_pos = dispersion['direct_gap_path_inds'][1]
                    vbm = dispersion['direct_valence_band_min']
                    cbm_pos = dispersion['direct_gap_path_inds'][0]
                    cbm = dispersion['direct_conduction_band_max']
                    vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                    cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                    ax_dispersion.plot([path[vbm_pos-vbm_offset], path[cbm_pos-cbm_offset]], [vbm, cbm], ls=ls[seed_ind], c='red', label='direct gap {:3.3f} eV'.format(cbm-vbm))
                    ax_dispersion.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True, ncol=2, handlelength=1)
                else:
                    ax_dispersion.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True, ncol=1, handlelength=1)
            ax_dispersion.set_xticks(xticks)
            ax_dispersion.set_xticklabels(xticklabels)
        if plot_dos:
            if not phonons:
                dos_data = np.loadtxt(seed + '.adaptive.dat')
                energies = dos_data[:, 0]
                dos = dos_data[:, 1]
                max_density = np.max(dos[np.where(energies > plot_window[0])])
            else:
                if not isfile(seed + '.phonon_dos'):
                    phonon_data, s = phonon2dict(seed + '.phonon')
                    plot_window = [np.min(phonon_data['eigenvalues_q'])-10, np.max(phonon_data['eigenvalues_q'])]
                    if s:
                        space_size = 1000
                        gaussian_width = 10
                        raw_weights = []
                        raw_eigenvalues = []
                        for qind, qpt in enumerate(phonon_data['eigenvalues_q']):
                            weight = phonon_data['qpoint_weights'][qind]
                            for eig in qpt:
                                raw_weights.append(weight)
                                raw_eigenvalues.append(eig)
                        hist, energies = np.histogram(raw_eigenvalues, weights=raw_weights, bins=space_size)
                        # shift bin edges to bin centres
                        energies -= energies[1] - energies[0]
                        energies = energies[:-1]
                        new_energies = np.reshape(energies, (1, len(energies))) - np.reshape(energies, (1, len(energies))).T
                        dos = np.sum(hist*np.exp(-(new_energies)**2 / gaussian_width), axis=1)
                        dos = np.divide(dos, np.sqrt(2 * np.pi * gaussian_width**2))
                        max_density = np.max(dos)
                        phonon_data['freq_unit'] = phonon_data['freq_unit'].replace('-1', '$^{-1}$')
                        ax_dos.axvline(phonon_data['softest_mode_freq'], ls='--', c='r',
                                       label='$\omega_\mathrm{{min}} = {:5.3f}$ {}'.format(phonon_data['softest_mode_freq'], phonon_data['freq_unit']))
                        ax_dos.legend()
                    else:
                        exit('Failed to read .phonon file')
                else:
                    with open(seed + '.phonon_dos', 'r') as f:
                        flines = f.readlines()
                    for ind, line in enumerate(flines):
                        if 'begin dos' in line.lower():
                            begin = ind+1
                            break
                    data_flines = flines[begin:-1]
                    with open(seed + '.phonon_dos_tmp', 'w') as f:
                        for line in data_flines:
                            f.write(line + '\n')
                    raw_data = np.loadtxt(seed + '.phonon_dos_tmp')
                    energies = raw_data[:, 0]
                    dos = raw_data[:, 1]
                    max_density = np.max(dos)
                    from os import remove
                    remove(seed + '.phonon_dos_tmp')
                    plot_window = [-10, 350]

            if phonons:
                ylabel = 'Phonon DOS'
                xlabel = 'Wavenumber (cm$^{{-1}}$)'
                # xlabel = 'Wavenumber ({})'.format(phonon_data['freq_unit'])
            else:
                ylabel = 'DOS'
                xlabel = 'Energy (eV)'
            if plot_bandstructure:
                ax_dos.set_xticks([0.6*max_density])
                ax_dos.set_xticklabels([ylabel])
                ax_dos.axhline(0, c='grey', ls='--', lw=1)
                ax_dos.set_xlim(0, max_density*1.2)
                ax_dos.set_ylim(plot_window)
                ax_dos.axvline(0, c='k')
                ax_dos.plot(dos, energies, lw=1, c='k', label='matador', ls=ls[seed_ind])
                if seed_ind == 0 and not phonons:
                    ax_dos.fill_betweenx(energies[np.where(energies <= 0)], 0, dos[np.where(energies <= 0)], facecolor=valence, alpha=0.5)
                    ax_dos.fill_betweenx(energies[np.where(energies >= 0)], 0, dos[np.where(energies >= 0)], facecolor=conduction, alpha=0.5)
            else:
                ax_dos.set_xlabel(xlabel)
                ax_dos.set_ylabel(ylabel)
                ax_dos.axvline(0, c='grey', ls='--', lw=1)
                ax_dos.set_ylim(0, max_density*1.2)
                ax_dos.set_xlim(plot_window)
                ax_dos.axhline(0, c='k')
                ax_dos.plot(energies, dos, lw=1, c='k', ls=ls[seed_ind])
                if seed_ind == 0 and not phonons:
                    ax_dos.fill_between(energies[np.where(energies <= 0)], 0, dos[np.where(energies <= 0)], facecolor=valence, alpha=0.5)
                    ax_dos.fill_between(energies[np.where(energies >= 0)], 0, dos[np.where(energies >= 0)], facecolor=conduction, alpha=0.5)
            ax_dos.grid('off')

    if kwargs.get('pdf'):
        plt.savefig(seed.replace('.bands', '').replace('.phonon', '') + '_spectral.pdf', bbox_inches='tight')
    if kwargs.get('png'):
        plt.savefig(seed.replace('.bands', '').replace('.phonon', '') + '_spectral.png', bbox_inches='tight', dpi=300)

    plt.show()


def modes2bands(phonon_dispersion, path, branches):
    """ Reorder phonon eigenvalues such that bands join up correctly,
    based on local gradients.

    Input:

        | phonon_dispersion : np.ndarray, containing eigenvalues as function of q
        | path              : np.ndarray, containing q-point path
        | branches          : list(int), containing branch start/end points of q-point path

    """
    import numpy as np
    return phonon_dispersion
    eigs = phonon_dispersion[0]
    for branch_ind, branch in enumerate(branches[:1]):
        eigs_branch = eigs[:, branch]
        for i in range(1, len(branch)-1):
            guess = 2*eigs_branch[:, i] - eigs_branch[:, i-1]
            if np.any(np.argsort(guess) != np.argsort(eigs_branch[:, i])):
                print(i)
                print('Actuals:')
                print(eigs_branch[:, i+1])
                print('Guesses:')
                print(guess)
                print('!!!!!')
                print(np.argsort(guess))
                for j in range(i+1, len(branch)):
                    eigs_branch[:, j] = eigs_branch[np.argsort(guess), j]
        eigs[:, branch] = eigs_branch

    phonon_dispersion[0] = eigs

    return phonon_dispersion
