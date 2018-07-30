# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot densities of states and
bandstructures for electronic and vibrational calculations.

"""


import os
import numpy as np
from matador.utils.viz_utils import ELEMENT_COLOURS
from matador.plotting.plotting import plotting_function


@plotting_function
def plot_spectral(seeds, **kwargs):
    """ Plot bandstructure and/or optional DOS from `<seed>.bands` and
    `<seed>.adaptive.dat` file.

    Parameters:
        seeds (list): list of filenames of bands/phonon files

    Keyword Arguments:
        plot_bandstructure (bool): whether to plot bandstructure
        plot_dos (bool): whether to plot density of states
        dos (str): separate seed name for pDOS/DOS data
        phonons (bool): whether to plot phonon or electronic data
        labels (list): list of strings for legend labels for multiple bandstructures
        gap (bool): whether to draw on the band gap
        colour_by_seed (bool): plot with a separate colour per bandstructure
        external_efermi (float): replace scraped Fermi energy with this value (eV)
        highlight_bands (list): list of integer indices, colour the bands with these indices in red
        band_colour (str): if passed "occ", bands will be coloured using cmap depending on whether
            they lie above or below the Fermi level. If passed 'random', colour bands randomly from
            the cmap. Otherwise, override all colour options with matplotlib-interpretable colour
            (e.g. hexcode or html colour name) to use for all bands (DEFAULT: 'occ').
        cmap (str): matplotlib colourmap name to use for the bands
        n_colours (int): number of colours to use from cmap (DEFAULT: 4).
        no_stacked_pdos (bool): whether to plot projected DOS as stack or overlapping.
        preserve_kspace_distance (bool): whether to preserve distances in reciprocal space when
            linearising the kpoint path. If False, bandstructures of different lattice parameters
            with the same Bravais lattice can be more easily compared. If True, bandstructures may
            appear rarefied or compressed in particular regions.
        band_reorder (bool): try to reorder bands based on local gradients (DEFAULT: True for phonons, otherwise False).
        title (str): optional plot title
        pdos_hide_tot (bool): whether or not to plot the total DOS on a PDOS plot; this is to hide
            regions where the PDOS is negative (leading to total DOS lower than stacked PDOS) (DEFAULT: False).

    """
    import matplotlib.pyplot as plt
    from cycler import cycler
    from os.path import isfile
    from matador.scrapers.castep_scrapers import bands2dict, cell2dict, phonon2dict, optados2dict
    from matador.utils.cell_utils import doc2spg
    from seekpath import get_path
    # set defaults and update class with desired values
    prop_defaults = {'plot_bandstructure': True, 'plot_dos': True,
                     'phonons': False, 'gap': False,
                     'colour_by_seed': False, 'external_efermi': None,
                     'labels': None, 'cmap': None, 'band_colour': 'occ',
                     'n_colours': 4,
                     'no_stacked_pdos': False, 'preserve_kspace_distance': False,
                     'band_reorder': None, 'title': None,
                     'verbosity': 0, 'highlight_bands': None, 'pdos_hide_tot': False}
    prop_defaults.update(kwargs)
    kwargs = prop_defaults

    if kwargs.get('cmap') is None:
        colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    else:
        if isinstance(kwargs['cmap'], str):
            print('Adjusting colour palette... to {}'.format(kwargs.get('cmap')))
            colours = plt.cm.get_cmap(kwargs.get('cmap')).colors
            plt.rcParams['axes.prop_cycle'] = cycler('color', colours)
        elif isinstance(kwargs['cmap'], list):
            print('Reading list of colours {}...'.format(kwargs.get('cmap')))
            colours = kwargs['cmap']
            plt.rcParams['axes.prop_cycle'] = cycler('color', colours)

    dos_legend = None
    disp_legend = None

    if (kwargs.get('phonons') and kwargs['band_colour'] == 'occ') or kwargs['band_colour'] == 'random':
        kwargs['band_colour'] = None

    if not isinstance(seeds, list):
        seeds = [seeds]

    if kwargs.get('plot_window') is not None:
        if isinstance(kwargs.get('plot_window'), list):
            if len(kwargs.get('plot_window')) != 2:
                exit('plot_window must have length 2 or be a single number')
            plot_window = sorted(kwargs.get('plot_window'))
        else:
            plot_window = (-kwargs.get('plot_window'), kwargs.get('plot_window'))
    else:
        plot_window = None

    if kwargs['plot_bandstructure'] and not kwargs['plot_dos']:
        fig, ax_dispersion = plt.subplots(figsize=(7, 6))
    elif kwargs['plot_bandstructure'] and kwargs['plot_dos']:
        fig, ax_grid = plt.subplots(1, 3, figsize=(8, 6), sharey=True,
                                    gridspec_kw={'width_ratios': [4, 1, 1],
                                                 'wspace': 0.05,
                                                 'left': 0.15})
        ax_dispersion = ax_grid[0]
        ax_dos = ax_grid[1]
        ax_grid[2].axis('off')
    elif not kwargs['plot_bandstructure'] and kwargs['plot_dos']:
        fig, ax_dos = plt.subplots(1, figsize=(8, 4))

    valence = colours[0]
    conduction = colours[-1]
    crossing = colours[int(len(colours) / 2)]

    if len(seeds) > 1 or kwargs.get('colour_by_seed'):
        seed_colours = colours
        ls = ['-'] * len(seeds)
        colour_by_seed = True
        if kwargs.get('labels') is None:
            kwargs['labels'] = [seed.split('/')[-1].split('.')[0] for seed in seeds]
    else:
        ls = []
        colour_by_seed = False
        for i in range(len(seeds)):
            if i % 3 == 0:
                ls.append('-')
            elif i % 3 == 1:
                ls.append('--')
            elif i % 3 == 2:
                ls.append('-.')

    for seed_ind, seed in enumerate(seeds):
        seed = seed.replace('.bands', '').replace('.phonon', '')
        if kwargs['plot_bandstructure']:
            if kwargs['phonons']:
                dispersion, s = phonon2dict(seed + '.phonon', verbosity=kwargs.get('verbosity'))
                branch_key = 'qpoint_branches'
                num_key = 'num_qpoints'
                path_key = 'qpoint_path'
                eig_key = 'eigenvalues_q'
                band_key = 'num_branches'
                dispersion['num_spins'] = 1
                spin_key = 'num_spins'
                if plot_window is None:
                    plot_window = [min(-10, np.min(dispersion[eig_key])-10), np.max(dispersion[eig_key])]
            else:
                dispersion, s = bands2dict(seed + '.bands',
                                           summary=True,
                                           gap=kwargs.get('gap'),
                                           external_efermi=kwargs.get('external_efermi'),
                                           verbosity=kwargs.get('verbosity'))
                if not s:
                    raise RuntimeError(dispersion)
                branch_key = 'kpoint_branches'
                num_key = 'num_kpoints'
                path_key = 'kpoint_path'
                eig_key = 'eigenvalues_k_s'
                band_key = 'num_bands'
                spin_key = 'num_spins'
                if plot_window is None:
                    plot_window = [-10, 10]

            path = [0]
            for branch in dispersion[branch_key]:
                for ind, kpt in enumerate(dispersion[path_key][branch]):
                    if ind != len(branch) - 1:
                        if kwargs['preserve_kspace_distance']:
                            diff = np.sqrt(np.sum((kpt - dispersion[path_key][branch[ind + 1]])**2))
                        else:
                            diff = 1.
                        path.append(path[-1] + diff)
            path = np.asarray(path)
            path /= np.max(path)
            assert len(path) == int(dispersion[num_key]) - len(dispersion[branch_key]) + 1
            if kwargs['band_reorder'] or (kwargs['band_reorder'] is None and kwargs['phonons']):
                print('Reordering bands based on local gradients...')
                dispersion[eig_key] = match_bands(dispersion[eig_key], dispersion[branch_key])

            # seem to have to reset this here for some reason
            for branch_ind, branch in enumerate(dispersion[branch_key]):
                plt.rcParams['axes.prop_cycle'] = cycler('color', colours)
                for ns in range(dispersion[spin_key]):
                    for nb in range(dispersion[band_key]):
                        colour = None
                        alpha = 1
                        label = None
                        if not kwargs['phonons']:
                            if dispersion[spin_key] == 2:
                                if ns == 0:
                                    colour = 'red'
                                    alpha = 0.3
                                else:
                                    colour = 'blue'
                                    alpha = 0.3
                            else:
                                if kwargs.get('band_colour') == 'occ':
                                    band_min = np.min(dispersion[eig_key][ns][nb][branch])
                                    band_max = np.max(dispersion[eig_key][ns][nb][branch])
                                    if band_max < 0:
                                        colour = valence
                                    elif band_min > 0:
                                        colour = conduction
                                    elif band_min < 0 and band_max > 0:
                                        colour = crossing

                        if colour_by_seed:
                            colour = seed_colours[seed_ind]

                        if kwargs.get('band_colour') is not None:
                            if kwargs.get('band_colour') != 'occ':
                                colour = kwargs.get('band_colour')

                        if kwargs.get('highlight_bands') is not None:
                            if nb in kwargs.get('highlight_bands'):
                                colour = 'red'
                            else:
                                alpha = 0.5
                        if branch_ind == 0 and ns == 0 and nb == 0 and kwargs.get('labels') is not None:
                            label = kwargs.get('labels')[seed_ind]

                        ax_dispersion.plot(path[(np.asarray(branch)-branch_ind).tolist()],
                                           dispersion[eig_key][ns][nb][branch], c=colour,
                                           lw=1, ls=ls[seed_ind], alpha=alpha, label=label)

            if len(seeds) > 1:
                disp_legend = ax_dispersion.legend(loc='upper center', facecolor='w',
                                                   frameon=True, fancybox=False, shadow=False, framealpha=1)

            ax_dispersion.axhline(0, ls='--', lw=1, c='grey')
            ax_dispersion.set_ylim(plot_window)
            if kwargs['phonons']:
                ylabel = 'Wavenumber (cm$^{-1}$)'
            else:
                ylabel = r'Energy (eV)'
            ax_dispersion.set_ylabel(ylabel)
            ax_dispersion.set_xlim(0, 1)
            xticks = []
            xticklabels = []
            shear_planes = []
            labelled = []

            # get dispersion path labels
            spg_structure = None
            if kwargs['phonons']:
                spg_structure = doc2spg(dispersion)
            else:
                if not os.path.isfile(seed + '.cell'):
                    print('Failed to find {}.cell, will not be able to generate labels.'.format(seed))

                doc, success = cell2dict(seed + '.cell',
                                         db=False, verbosity=kwargs['verbosity'],
                                         outcell=True, positions=True)
                if success:
                    spg_structure = doc2spg(doc)
                else:
                    print('Failed to scrape {}.cell, will not be able to generate labels.'.format(seed))

            if spg_structure is not False and spg_structure is not None:
                seekpath_results = get_path(spg_structure)
                path_labels = seekpath_results['point_coords']

            for branch_ind, branch in enumerate(dispersion[branch_key]):
                for sub_ind, ind in enumerate(branch):
                    kpt = dispersion[path_key][ind]
                    for label, point in path_labels.items():
                        if np.allclose(point, kpt):
                            if ind - branch_ind not in labelled:
                                label = label.replace('GAMMA', r'\Gamma')
                                label = label.replace('SIGMA', r'\Sigma')
                                label = label.replace('DELTA', r'\Delta')
                                label = label.replace('LAMBDA', r'\Lambda')
                                if sub_ind == len(branch) - 1:
                                    if branch_ind < len(dispersion[branch_key]) - 1:
                                        _tmp = dispersion[path_key]
                                        next_point = _tmp[dispersion[branch_key][branch_ind + 1][0]]
                                        for new_label, new_point in path_labels.items():
                                            new_label = new_label.replace('GAMMA', r'\Gamma')
                                            new_label = new_label.replace('SIGMA', r'\Sigma')
                                            new_label = new_label.replace('DELTA', r'\Delta')
                                            new_label = new_label.replace('LAMBDA', r'\Lambda')
                                            if np.allclose(new_point, next_point):
                                                label = '{}|{}'.format(label, new_label)
                                                ax_dispersion.axvline(path[ind - branch_ind], ls='-', c='grey', zorder=1, lw=0.5)
                                                labelled.append(ind - branch_ind)
                                                shear_planes.append(ind)
                                label = '${}$'.format(label)
                                ax_dispersion.axvline(path[ind - branch_ind], ls='--', c='grey', zorder=0, lw=0.5)
                                xticklabels.append(label)
                                xticks.append(path[ind - branch_ind])
                                break

            # plot band gaps
            if not kwargs['phonons'] and kwargs['gap'] and dispersion['band_gap'] > 0:
                vbm_pos = dispersion['band_gap_path_inds'][1]
                vbm = dispersion['valence_band_min']
                cbm_pos = dispersion['band_gap_path_inds'][0]
                cbm = dispersion['conduction_band_max']
                if vbm_pos != cbm_pos:
                    vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                    cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                    ax_dispersion.plot([path[vbm_pos - vbm_offset], path[cbm_pos - cbm_offset]], [vbm, cbm],
                                       ls=ls[seed_ind],
                                       c='blue',
                                       label='indirect gap {:3.3f} eV'.format(cbm - vbm))

                vbm_pos = dispersion['direct_gap_path_inds'][1]
                vbm = dispersion['direct_valence_band_min']
                cbm_pos = dispersion['direct_gap_path_inds'][0]
                cbm = dispersion['direct_conduction_band_max']
                vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                ax_dispersion.plot([path[vbm_pos - vbm_offset], path[cbm_pos - cbm_offset]], [vbm, cbm],
                                   ls=ls[seed_ind],
                                   c='red',
                                   label='direct gap {:3.3f} eV'.format(cbm - vbm))
                ax_dispersion.legend(loc='upper center',
                                     bbox_to_anchor=(0.5, 1.1),
                                     fancybox=True, shadow=True,
                                     ncol=2, handlelength=1)

            ax_dispersion.set_xticks(xticks)
            ax_dispersion.set_xticklabels(xticklabels)
            ax_dispersion.grid(False)

        if kwargs['plot_dos']:
            if not kwargs['phonons']:
                if kwargs.get('dos') is None:
                    # look for dat files, and just use the first
                    import glob
                    dos_seed = glob.glob('{}*.dat'.format(seed))
                    dos_data, s = optados2dict(dos_seed[0])
                else:
                    dos_data, s = optados2dict(kwargs.get('dos'))
                energies = dos_data['energies']
                dos = dos_data['dos']
                if 'spin_dos' in dos_data:
                    max_density = max(np.max(np.abs(dos_data['spin_dos']['down'][np.where(energies > plot_window[0])])),
                                      np.max(np.abs(dos_data['spin_dos']['up'][np.where(energies > plot_window[0])])))
                else:
                    max_density = np.max(dos[np.where(energies > plot_window[0])])

                if 'pdos' in dos_data:
                    pdos = dos_data['pdos']
            else:
                if not isfile(seed + '.phonon_dos'):
                    phonon_data, s = phonon2dict(seed + '.phonon')
                    if not s:
                        raise RuntimeError(phonon_data)
                    else:
                        if plot_window is None:
                            plot_window = [min(-10, np.min(phonon_data['eigenvalues_q']) - 10), np.max(phonon_data['eigenvalues_q'])]
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
                        new_energies = np.reshape(energies, (1, len(energies)))
                        new_energies -= np.reshape(energies, (1, len(energies))).T
                        dos = np.sum(hist * np.exp(-(new_energies)**2 / gaussian_width), axis=1)
                        dos = np.divide(dos, np.sqrt(2 * np.pi * gaussian_width**2))
                        max_density = np.max(dos)
                        phonon_data['freq_unit'] = phonon_data['freq_unit'].replace('-1', '$^{-1}$')
                        ax_dos.axvline(phonon_data['softest_mode_freq'], ls='--', c='r',
                                       label=(r'$\omega_\mathrm{{min}} = {:5.3f}$ {}'
                                              .format(phonon_data['softest_mode_freq'],
                                                      phonon_data['freq_unit'])))
                else:
                    with open(seed + '.phonon_dos', 'r') as f:
                        flines = f.readlines()
                    for ind, line in enumerate(flines):
                        if 'begin dos' in line.lower():
                            projector_labels = line.split()[5:]
                            projector_labels = [(label, None) for label in projector_labels]
                            begin = ind + 1
                            break
                    data_flines = flines[begin:-1]
                    with open(seed + '.phonon_dos_tmp', 'w') as f:
                        for line in data_flines:
                            f.write(line + '\n')
                    raw_data = np.loadtxt(seed + '.phonon_dos_tmp')
                    energies = raw_data[:, 0]
                    dos = raw_data[:, 1]
                    dos_data = {}
                    dos_data['dos'] = dos
                    max_density = np.max(dos)
                    dos_data['energies'] = energies
                    if plot_window is None:
                        plot_window = [np.min(energies[np.where(dos > 1e-3)]) - 10, np.max(energies[np.where(dos > 1e-3)])]
                    dos_data['pdos'] = dict()
                    for i, label in enumerate(projector_labels):
                        dos_data['pdos'][label] = raw_data[:, i + 2]
                    pdos = dos_data['pdos']
                    from os import remove
                    remove(seed + '.phonon_dos_tmp')

            if kwargs['phonons']:
                ylabel = 'Phonon DOS'
                xlabel = 'Wavenumber (cm$^{{-1}}$)'
            else:
                if kwargs['plot_bandstructure']:
                    ylabel = 'DOS'
                else:
                    ylabel = 'DOS (eV$^{{-1}}$Å$^{{-3}}$)'
                xlabel = 'Energy (eV)'

            if len(seeds) > 1:
                colour = seed_colours[seed_ind]
            else:
                colour = 'k'

            ax_dos.grid(False)

            if kwargs['plot_bandstructure']:
                ax_dos.set_xticks([0.6 * max_density])
                ax_dos.set_xticklabels([ylabel])
                ax_dos.axhline(0, c='grey', ls='-', lw=0.5)
                if 'spin_dos' in dos_data:
                    ax_dos.set_xlim(-max_density*1.2, max_density * 1.2)
                else:
                    ax_dos.set_xlim(0, max_density * 1.2)
                ax_dos.set_ylim(plot_window)
                ax_dos.axvline(0, c='k')

                if 'pdos' not in dos_data and 'spin_dos' not in dos_data:
                    ax_dos.plot(dos, energies, lw=1, ls=ls[seed_ind], color=colour, zorder=1e10, label='Total DOS')
                    if 'pdos' not in dos_data and 'spin_dos' not in dos_data:
                        ax_dos.fill_betweenx(energies, 0, dos, alpha=0.2, color=colour)

            else:
                ax_dos.set_xlabel(xlabel)
                ax_dos.set_ylabel(ylabel)
                ax_dos.axvline(0, c='grey', lw=0.5)
                if 'spin_dos' in dos_data:
                    ax_dos.set_ylim(-max_density*1.2, max_density * 1.2)
                else:
                    ax_dos.set_ylim(0, max_density * 1.2)
                ax_dos.set_xlim(plot_window)
                ax_dos.axhline(0, c='grey', lw=0.5)

                if 'pdos' not in dos_data and 'spin_dos' not in dos_data:
                    ax_dos.plot(energies, dos, lw=1, ls=ls[seed_ind], alpha=1, color=colour, zorder=1e10, label='Total DOS')
                    if 'pdos' not in dos_data and 'spin_dos' not in dos_data:
                        ax_dos.fill_between(energies, 0, dos, alpha=0.2, color=colour)

            if 'pdos' in dos_data and len(seeds) == 1:
                dos_colours = []
                for ind, projector in enumerate(pdos):
                    if ind == 0:
                        stack = np.zeros_like(pdos[projector])

                    if projector[0] is None:
                        projector_label = '${}$-character'.format(projector[1])
                    elif projector[1] is None:
                        projector_label = projector[0]
                    else:
                        projector_label = '{p[0]} (${p[1]}$)'.format(p=projector)

                    # if species-projected only, then use VESTA colours
                    if projector[0] is not None and projector[1] is None:
                        dos_colours.append(ELEMENT_COLOURS.get(projector[0]))
                    # if species_ang-projected, then use VESTA colours but lightened
                    elif projector[0] is not None and projector[1] is not None:
                        from copy import deepcopy
                        dos_colour = deepcopy(ELEMENT_COLOURS.get(projector[0]))
                        multi = ['s', 'p', 'd', 'f'].index(projector[1]) - 1
                        for jind, _ in enumerate(dos_colour):
                            dos_colour[jind] = max(min(dos_colour[jind]+multi*0.2, 1), 0)
                        dos_colours.append(dos_colour)
                    # otherwise if just ang-projected, use colour_cycle
                    else:
                        dos_colours.append(list(plt.rcParams['axes.prop_cycle'].by_key()['color'])[ind])

                    if not kwargs['no_stacked_pdos']:
                        alpha = 0.8
                    else:
                        alpha = 0.7

                    # mask negative contributions with 0
                    pdos[projector] = np.ma.masked_where(pdos[projector] < 0, pdos[projector], copy=True)
                    np.ma.set_fill_value(pdos[projector], 0)
                    pdos[projector] = np.ma.filled(pdos[projector])

                    if not np.max(pdos[projector]) < 1e-8:

                        if kwargs['plot_bandstructure']:
                            ax_dos.plot(stack+pdos[projector], energies, lw=1, zorder=1000, color=dos_colours[-1])
                            ax_dos.fill_betweenx(energies, stack, stack+pdos[projector], alpha=alpha, label=projector_label,
                                                 color=dos_colours[-1])
                        else:
                            ax_dos.plot(energies, stack+pdos[projector], lw=1, zorder=1000, color=dos_colours[-1])
                            ax_dos.fill_between(energies, stack, stack+pdos[projector], alpha=alpha, label=projector_label,
                                                color=dos_colours[-1])

                        if not kwargs['no_stacked_pdos']:
                            stack += pdos[projector]

                if not kwargs['pdos_hide_tot'] and not kwargs['no_stacked_pdos']:
                    if kwargs['plot_bandstructure']:
                        ax_dos.plot(stack, energies, lw=1, ls=ls[seed_ind], alpha=1, color=colour, zorder=1e10, label='Total DOS')
                    else:
                        ax_dos.plot(energies, stack, lw=1, ls=ls[seed_ind], alpha=1, color=colour, zorder=1e10, label='Total DOS')

            elif 'spin_dos' in dos_data:
                print('Plotting spin dos')
                if kwargs['plot_bandstructure']:
                    ax_dos.plot(dos_data['spin_dos']['up'], energies, lw=1, ls=ls[seed_ind], color='r', zorder=1e10, label='spin-up channel', alpha=alpha)
                    ax_dos.plot(dos_data['spin_dos']['down'], energies, lw=1, ls=ls[seed_ind], color='b', zorder=1e10, label='spin-down channel', alpha=alpha)
                    ax_dos.fill_betweenx(energies, 0, dos_data['spin_dos']['up'], alpha=0.2, color='r')
                    ax_dos.fill_betweenx(energies, 0, dos_data['spin_dos']['down'], alpha=0.2, color='b')
                else:
                    ax_dos.plot(energies, dos_data['spin_dos']['up'], lw=1, ls=ls[seed_ind], color='r', zorder=1e10, label='spin-up channel')
                    ax_dos.plot(energies, dos_data['spin_dos']['down'], lw=1, ls=ls[seed_ind], color='b', zorder=1e10, label='spin-down channel')
                    ax_dos.fill_between(energies, 0, dos_data['spin_dos']['up'], alpha=0.2, color='r')
                    ax_dos.fill_between(energies, 0, dos_data['spin_dos']['down'], alpha=0.2, color='b')

            if len(seeds) == 1:
                dos_legend = ax_dos.legend(bbox_to_anchor=(1, 1), facecolor='w', frameon=True, fancybox=False, shadow=False)

    if kwargs.get('title') is not None:
        fig.suptitle(kwargs.get('title'))

    if any([kwargs.get('pdf'), kwargs.get('svg'), kwargs.get('png')]):
        bbox_extra_artists = []
        if dos_legend is not None:
            bbox_extra_artists.append(dos_legend)
        if disp_legend is not None:
            bbox_extra_artists.append(disp_legend)
        if not bbox_extra_artists:
            bbox_extra_artists = None
        if kwargs.get('pdf'):
            plt.savefig(seeds[0].replace('.bands', '').replace('.phonon', '') + '_spectral.pdf',
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
        if kwargs.get('svg'):
            plt.savefig(seeds[0].replace('.bands', '').replace('.phonon', '') + '_spectral.svg',
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
        if kwargs.get('png'):
            plt.savefig(seeds[0].replace('.bands', '').replace('.phonon', '') + '_spectral.png',
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)

    else:
        print('Displaying plot...')
        plt.show()


def match_bands(dispersion, branches):
    """ Recursively reorder eigenvalues such that bands join up correctly,
    based on local gradients.

    Parameters:
        dispersion (numpy.ndarray): array containing eigenvalues as a
            function of q/k
        branches (:obj:`list` of :obj:`int`): list containing branches of
            k/q-point path

    Returns:
        numpy.ndarray: reordered branches.

    """
    from copy import deepcopy

    for channel_ind, channel in enumerate(dispersion):
        eigs = channel
        for branch_ind, branch in enumerate(branches):
            eigs_branch = eigs[:, branch]
            converged = False
            counter = 0
            i_cached = 0
            while not converged and counter < len(branch):
                counter += 1
                for i in range(i_cached+1, len(branch) - 1):
                    guess = (2 * eigs_branch[:, i] - eigs_branch[:, i-1])
                    argsort_guess = np.argsort(guess)
                    if np.any(np.argsort(guess) != np.argsort(eigs_branch[:, i+1])):
                        tmp_copy = deepcopy(eigs)
                        for ind, mode in enumerate(np.argsort(eigs_branch[:, i]).tolist()):
                            eigs_branch[mode, i+1:] = tmp_copy[:, branch][argsort_guess[ind], i+1:]
                        for other_branch in branches[branch_ind:]:
                            eigs_other_branch = eigs[:, other_branch]
                            for ind, mode in enumerate(np.argsort(eigs_branch[:, i]).tolist()):
                                eigs_other_branch[mode] = tmp_copy[:, other_branch][argsort_guess[ind]]
                            eigs[:, other_branch] = eigs_other_branch
                        eigs[:, branch] = eigs_branch
                        i_cached = i
                        break

                    if i == len(branch) - 2:
                        converged = True

        dispersion[channel_ind] = eigs.reshape(1, len(eigs), len(eigs[0]))

    return dispersion
