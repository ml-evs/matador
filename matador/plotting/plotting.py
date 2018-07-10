# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements several routines for dispersion plots, phase
diagrams, as well as voltage and volume expansion plots.

"""


import os
import numpy as np
from matador.utils.chem_utils import get_formula_from_stoich
from matador.utils.viz_utils import ELEMENT_COLOURS


def plotting_function(function):
    """ Wrapper for plotting functions to safely fail on X-forwarding
    errors.
    """

    from functools import wraps
    from matador.utils.print_utils import print_warning, print_failure
    from matador.config import load_custom_settings

    @wraps(function)
    def wrapped_plot_function(*args, **kwargs):
        """ Wrap and return the plotting function. """
        from tkinter import TclError
        result = None
        # if we're going to be saving a figure, switch to Agg to avoid X-forwarding
        saving = False

        try:
            for arg in args:
                if arg.savefig:
                    import matplotlib
                    matplotlib.use('Agg')
                    saving = True
                    break
        except AttributeError:
            pass
        if not saving:
            if any([kwargs.get('pdf'), kwargs.get('svg'), kwargs.get('png')]):
                import matplotlib
                matplotlib.use('Agg')
                saving = True

        settings = load_custom_settings(kwargs.get('config_fname'), quiet=False, override=kwargs.get('override'))
        try:
            import matplotlib.pyplot as plt
            style = settings.get('plotting', {}).get('default_style')
            if style is None or style == 'matador':
                style = '/'.join(__file__.split('/')[:-1]) + '/../config/matador.mplstyle'
            if kwargs.get('debug'):
                print('Using style {}'.format(style))
            plt.style.use(style)
            result = function(*args, **kwargs)
        except TclError as exc:
            print_failure('Caught exception: {}'.format(type(exc).__name__))
            print_warning('Error message was: {}'.format(exc))
            print_warning('This is probably an X-forwarding error')
            print_failure('Skipping plot...')

        return result

    return wrapped_plot_function


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
            print('Reading list of colours...'.format(kwargs.get('cmap')))
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

                if not kwargs['pdos_hide_tot'] and 'spin_dos' not in dos_data:
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

                if not kwargs['pdos_hide_tot'] and 'spin_dos' not in dos_data:
                    ax_dos.plot(energies, dos, lw=1, ls=ls[seed_ind], alpha=1, color=colour, zorder=1e10, label='Total DOS')
                    if 'pdos' not in dos_data and 'spin_dos' not in dos_data:
                        ax_dos.fill_between(energies, 0, dos, alpha=0.2, color=colour)

            if 'pdos' in dos_data and len(seeds) == 1:
                dos_colours = []
                for ind, projector in enumerate(pdos):
                    if ind == 0:
                        stack = np.zeros_like(pdos[projector])

                    if projector[0] is None:
                        projector_label = projector[1]
                    elif projector[1] is None:
                        projector_label = projector[0]
                    else:
                        projector_label = '{p[0]} ({p[1]})'.format(p=projector)

                    # if species-projected only, then use VESTA colours
                    if projector[0] is not None and projector[1] is None:
                        dos_colours.append(ELEMENT_COLOURS.get(projector[0]))
                    # if species_ang-projected only, then use VESTA colours but lightened
                    elif projector[0] is not None and projector[1] is not None:
                        from copy import deepcopy
                        dos_colour = deepcopy(ELEMENT_COLOURS.get(projector[0]))
                        multi = ['s', 'p', 'd', 'f'].index(projector[1]) - 1
                        for jind, _ in enumerate(dos_colour):
                            dos_colour[jind] = max(min(dos_colour[jind]+multi*0.2, 1), 0)
                        dos_colours.append(dos_colour)
                    else:
                        dos_colours.append(None)

                    if not kwargs['no_stacked_pdos']:
                        alpha = 0.8
                    else:
                        alpha = 0.7

                    # mask negative contributions with 0
                    pdos[projector] = np.ma.masked_where(pdos[projector] < 0, pdos[projector], copy=True)
                    np.ma.set_fill_value(pdos[projector], 0)
                    pdos[projector] = np.ma.filled(pdos[projector])

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


def get_hull_labels(hull, label_cutoff=0.0, num_species=2):
    """ Return list of structures to labels on phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): phase diagram to plot.

    Keyword arguments:
        num_species (int): structures containing this number of species
            will be labelled.

    Returns:
        label_cursor (list(dict)): list of matador documents to label.

    """
    eps = 1e-9
    if isinstance(label_cutoff, list) and len(label_cutoff) == 2:
        label_cutoff = sorted(label_cutoff)
        # first, only apply upper limit as we need to filter by stoich aftewards
        label_cursor = [doc for doc in hull.hull_cursor if doc['hull_distance'] <= label_cutoff[1]]
    else:
        if isinstance(label_cutoff, list):
            assert len(label_cutoff) == 1, 'Incorrect number of label_cutoff values passed, should be 1 or 2.'
            label_cutoff = label_cutoff[0]
        label_cursor = [doc for doc in hull.hull_cursor if doc['hull_distance'] <= label_cutoff + eps]

    num_labels = len(set([get_formula_from_stoich(doc['stoichiometry']) for doc in label_cursor]))
    if num_labels < len(label_cursor):
        tmp_cursor = []
        for doc in label_cursor:
            if doc['stoichiometry'] not in [_doc['stoichiometry'] for _doc in tmp_cursor]:
                tmp_cursor.append(doc)
            else:
                label_cursor = tmp_cursor
    if isinstance(label_cutoff, list) and len(label_cutoff) == 2:
        # now apply lower limit
        label_cursor = [doc for doc in label_cursor if label_cutoff[0] <= doc['hull_distance'] <= label_cutoff[1]]
    # remove chemical potentials and unwanted e.g. binaries
    label_cursor = [doc for doc in label_cursor if len(doc['stoichiometry']) == num_species]

    return label_cursor


@plotting_function
def plot_voltage_curve(hull, ax=None, show=False):
    """ Plot voltage curve calculated for phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        ax (matplotlib.axes.Axes): an existing axis on which to plot.
        show (bool): whether to show plot in an X window.

    """
    import matplotlib.pyplot as plt
    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    if ax is None:
        if hull.savefig:
            if len(hull.voltage_data['voltages']) != 1:
                fig = plt.figure(facecolor=None, figsize=(4, 3.5))
            else:
                fig = plt.figure(facecolor=None, figsize=(4, 3.5))
        else:
            fig = plt.figure(facecolor=None)
        axQ = fig.add_subplot(111)
    else:
        axQ = ax
    if hull.args.get('expt') is not None:
        expt_data = np.loadtxt(hull.args.get('expt'), delimiter=',')
        if hull.args.get('expt_label'):
            axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label=hull.args.get('expt_label'))
        else:
            axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label='Experiment')
    for ind, voltage in enumerate(hull.voltage_data['voltages']):
        for i in range(1, len(voltage) - 1):
            if i == 1 and hull.args.get('expt'):
                axQ.plot([hull.voltage_data['Q'][ind][i - 1], hull.voltage_data['Q'][ind][i]], [voltage[i], voltage[i]],
                         marker='*',
                         lw=2,
                         c=hull.colours[ind],
                         label='DFT (this work)')
            elif i == 1 and len(hull.voltage_data['voltages']) != 1:
                axQ.plot([hull.voltage_data['Q'][ind][i - 1], hull.voltage_data['Q'][ind][i]], [voltage[i], voltage[i]],
                         # marker='o',
                         markersize=5,
                         lw=2,
                         c=hull.colours[ind],
                         label=get_formula_from_stoich(hull.endstoichs[ind], tex=True))
            else:
                axQ.plot([hull.voltage_data['Q'][ind][i - 1], hull.voltage_data['Q'][ind][i]], [voltage[i], voltage[i]],
                         # marker='o',
                         markersize=5,
                         lw=2,
                         c=hull.colours[ind])
            if i != len(voltage) - 2:
                axQ.plot([hull.voltage_data['Q'][ind][i], hull.voltage_data['Q'][ind][i]], [voltage[i], voltage[i + 1]],
                         # marker='o',
                         markersize=5,
                         lw=2,
                         c=hull.colours[ind])
    if hull.args.get('labels') or hull.args.get('label_cutoff') is not None:
        label_cursor = get_hull_labels(hull, num_species=2)
        for i in range(len(label_cursor)):
            axQ.annotate(get_formula_from_stoich(label_cursor[i]['stoichiometry'],
                                                 elements=hull.elements, tex=True),
                         xy=(hull.voltage_data['Q'][0][i+1]+0.02*max(hull.voltage_data['Q'][0]),
                             hull.voltage_data['voltages'][0][i+1]+0.02*max(hull.voltage_data['voltages'][0])),
                         textcoords='data',
                         ha='center',
                         zorder=9999)
    if hull.args.get('expt') or len(hull.voltage_data['voltages']) != 1:
        axQ.legend(loc=1)
    axQ.set_ylabel('Voltage (V) vs {}$^+/${}'.format(hull.elements[0], hull.elements[0]))
    axQ.set_xlabel('Gravimetric cap. (mAh/g)')
    _, end = axQ.get_ylim()
    from matplotlib.ticker import MultipleLocator
    axQ.yaxis.set_major_locator(MultipleLocator(0.2))
    axQ.set_ylim(0, 1.1 * end)
    _, end = axQ.get_xlim()
    axQ.set_xlim(0, 1.1 * end)
    axQ.grid(False)
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)

    if hull.savefig:
        if hull.args.get('pdf'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_voltage.pdf', dpi=500, transparent=True)
        if hull.args.get('svg'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_voltage.svg', dpi=500, transparent=True)
        if hull.args.get('png'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_voltage.png', dpi=500, transparent=True)
    elif show:
        plt.show()

    return axQ


@plotting_function
def plot_thermo_curves(seed, show=False, **kwargs):
    """ Plot free energies and heat capacity as a function of temperature, as
    calculated from a CASTEP thermodynamics run in <seed>.castep.

    Parameters:
        seed (str): filename of thermodynamics <seed>.castep file

    Keyword Arguments:
        plot_energy (bool): whether to plot T vs Energy
        plot_heat_cap (bool): whether to plot T vs Cv

    """
    import matplotlib.pyplot as plt
    from matador.scrapers.castep_scrapers import castep2dict
    from matador.utils.chem_utils import AVOGADROS_NUMBER, ELECTRON_CHARGE

    # set defaults
    prop_defaults = {'plot_energy': True, 'plot_heat_cap': True}
    prop_defaults.update(kwargs)
    kwargs = prop_defaults

    if kwargs['plot_energy'] and not kwargs['plot_heat_cap']:
        _, ax_energy = plt.subplots(figsize=(5, 5))
    elif kwargs['plot_energy'] and kwargs['plot_heat_cap']:
        _, ax_grid = plt.subplots(1, 2, figsize=(10, 5), sharey=False,
                                  gridspec_kw={'width_ratios': [1, 1],
                                               'wspace': 0.15,
                                               'left': 0.15})
        ax_energy = ax_grid[0]
        ax_Cv = ax_grid[1]
    elif not kwargs['plot_energy'] and kwargs['plot_heat_cap']:
        _, ax_Cv = plt.subplots(figsize=(5, 5))

    data, success = castep2dict(seed, db=False)
    if not success:
        raise RuntimeError('Failed to read {}: {}'.format(seed, data))

    temps = data['thermo_temps']
    Cv = data['thermo_heat_cap_Cv']
    free_energy = data['thermo_free_energy_F']
    entropy = data['thermo_entropy_S']
    enthalpy = data['thermo_enthalpy_E']
    # E0 = data['zero_point_E']

    entropyT = []
    for i, val in enumerate(entropy):
        entropyT.append(temps[i] * val * (1. / ELECTRON_CHARGE) * (1. / AVOGADROS_NUMBER))  # J/mol/K --> eV

    if kwargs['plot_energy']:
        ax_energy.plot(temps, free_energy, marker='o', lw=1, ls='-', alpha=1, label='Free Energy')
        ax_energy.plot(temps, entropyT, marker='o', lw=1, ls='-', alpha=1, label='Entropy*T')
        ax_energy.plot(temps, enthalpy, marker='o', lw=1, ls='-', alpha=1, label='Enthalpy')
        ax_energy.set_xlabel('Temperature (K)')
        ax_energy.set_title('Energy (eV)')
        ax_energy.legend()

    if kwargs['plot_heat_cap']:
        ax_Cv.plot(temps, Cv, marker='o', lw=1, ls='-', alpha=1, label='Heat Capacity (J/mol/K)')
        ax_Cv.set_title('Heat Capacity (J/mol/K)')
        ax_Cv.set_xlabel('Temperature (K)')

    if any([kwargs.get('pdf'), kwargs.get('svg'), kwargs.get('png')]):
        if kwargs.get('pdf'):
            plt.savefig(seed.replace('.castep', '') + '_thermoplots.pdf', dpi=300, transparent=True)
        if kwargs.get('svg'):
            plt.savefig(seed.replace('.castep', '') + '_thermoplots.svg', dpi=300, transparent=True)
        if kwargs.get('png'):
            plt.savefig(seed.replace('.castep', '') + '_thermoplots.png', dpi=300, transparent=True)

    elif show:
        plt.show()


def get_linear_cmap(colours, num_colours=100, list_only=False):
    """ Create a linear colormap from a list of colours.

    Parameters:
        colours (:obj:`list` of :obj:`str`): list of fractional RGB/hex
            values of colours

    Keyword arguments:
        num_colours (int): number of colours in resulting cmap
        list_only (bool): return only a list of colours

    Returns:
        :obj:`matplotlib.colors.LinearSegmentedColormap` or :obj:`list`:
            returns list of colours if `list_only` is True, otherwise
            :obj:`matplotlib.colors.LinearSegmentedColormap`.

    """
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    colours = [to_rgb(colour) for colour in colours]

    uniq_colours = []
    _colours = [tuple(colour) for colour in colours]
    for colour in _colours:
        if colour not in uniq_colours:
            uniq_colours.append(colour)
    _colours = uniq_colours
    linear_cmap = []
    repeat = int(num_colours / len(_colours))
    for ind, colour in enumerate(_colours):
        if ind == len(_colours) - 1:
            break
        diff = np.asarray(_colours[ind + 1]) - np.asarray(_colours[ind])
        diff_norm = diff / repeat
        for i in range(repeat):
            linear_cmap.append(np.asarray(colour) + i * diff_norm)
    if list_only:
        return linear_cmap

    return LinearSegmentedColormap.from_list('linear_cmap', linear_cmap, N=num_colours)


@plotting_function
def plot_volume_curve(hull, show=False):
    """ Plot volume curve calculated for phase diagram.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        show (bool): whether or not to display plot in X-window.

    """
    import matplotlib.pyplot as plt
    from matador.utils.cursor_utils import get_array_from_cursor
    from matador.utils.chem_utils import get_generic_grav_capacity
    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    if hull.savefig:
        fig = plt.figure(facecolor=None, figsize=(4, 3.5))
    else:
        fig = plt.figure(facecolor=None)
    ax = fig.add_subplot(111)
    stable_hull_dist = get_array_from_cursor(hull.hull_cursor, 'hull_distance')

    hull_vols = []
    hull_comps = []
    for i in range(len(hull.volume_data['vol_per_y'])):
        if stable_hull_dist[i] <= 0 + 1e-16:
            hull_vols.append(hull.volume_data['volume_ratio_with_bulk'][i])
            hull_comps.append(hull.volume_data['x'][i])
            s = 40
            zorder = 1000
            markeredgewidth = 1.5
            c = hull.colours[1]
            alpha = 1
        else:
            s = 30
            zorder = 900
            alpha = 0.3
            markeredgewidth = 0
            c = 'grey'

        ax.scatter(hull.volume_data['x'][i] / (1 + hull.volume_data['x'][i]),
                   hull.volume_data['volume_ratio_with_bulk'][i],
                   marker='o', s=s, edgecolor='k',
                   lw=markeredgewidth, c=c, zorder=zorder, alpha=alpha)

    hull_comps, hull_vols = np.asarray(hull_comps), np.asarray(hull_vols)
    ax.plot(hull_comps / (1 + hull_comps), hull_vols, marker='o', lw=4, c=hull.colours[0], zorder=100)

    ax.set_xlabel(r'$x$ in ' + hull.elements[0] + '$_x$' + hull.elements[1] + '$_{1-x}$')
    ax.set_ylabel('Volume ratio with bulk {}'.format(hull.volume_data['bulk_species']))
    ax.set_ylim(0, np.max(hull.volume_data['volume_ratio_with_bulk']))
    ax.set_xlim(-0.05, 1.05)
    ax.yaxis.set_label_position('left')
    ax.grid(False)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    tick_locs = [0, 0.2, 0.4, 0.6, 0.8, 1]
    ax2.set_xticks(tick_locs)
    new_tick_labels = [int(get_generic_grav_capacity([loc, 1-loc], [hull.elements[0], hull.elements[1]]))
                       for loc in tick_locs[:-1]]
    new_tick_labels[0] = 0
    new_tick_labels.append(r'$\infty$')
    ax2.set_xlabel('Gravimetric capacity (mAh/g)')
    ax2.set_xticklabels(new_tick_labels)
    ax2.grid(False)
    dark_grey = '#262626'
    for spine in ['left', 'top', 'right', 'bottom']:
        ax.spines[spine].set_color(dark_grey)
        ax2.spines[spine].set_color(dark_grey)
        ax.spines[spine].set_linewidth(0.5)
        ax2.spines[spine].set_linewidth(0.5)
    # ax.yaxis.set_ticks(range(0, int(end)+1, 5))
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
    if hull.savefig:
        if hull.args.get('pdf'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_volume.pdf', dpi=300)
        if hull.args.get('svg'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_volume.svg', dpi=300)
        if hull.args.get('png'):
            plt.savefig(hull.elements[0] + hull.elements[1] + '_volume.png', dpi=300, bbox_inches='tight')
    elif show:
        plt.show()


@plotting_function
def plot_2d_hull(hull, ax=None, show=False, plot_points=True,
                 plot_hull_points=True, labels=None, label_cutoff=None, colour_by_source=False,
                 sources=None, source_labels=None, title=True,
                 **kwargs):
    """ Plot calculated hull, returning ax and fig objects for further editing.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        ax (matplotlib.axes.Axes): an existing axis on which to plot,
        show (bool): whether or not to display the plot in an X window,
        plot_points (bool): whether or not to display off-hull structures,
        plot_hull_points (bool): whether or not to display on-hull structures,
        labels (bool): whether to label formulae of hull structures, also read from
            hull.args.
        label_cutoff (float/:obj:`tuple` of :obj:`float`): draw labels less than or
            between these distances form the hull, also read from hull.args.
        colour_by_source (bool): plot and label points by their sources
        alpha (float): alpha value of points when colour_by_source is True
        sources (list): list of possible provenances to colour when colour_by_source
            is True (others will be grey)
        title (str/bool): whether to include a plot title.

    Returns:
        matplotlib.axes.Axes: matplotlib axis with plot.

    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as colours

    if ax is None:
        fig = plt.figure(facecolor=None, figsize=(8, 6))
        ax = fig.add_subplot(111)

    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    hull.default_cmap_list = get_linear_cmap(hull.colours[1:4], list_only=True)
    hull.default_cmap = get_linear_cmap(hull.colours[1:4], list_only=False)

    if labels is None:
        labels = hull.args.get('labels', False)
    if label_cutoff is None:
        label_cutoff = hull.args.get('label_cutoff')

    scale = 1
    scatter = []
    x_elem = [hull.elements[0]]
    one_minus_x_elem = list(hull.elements[1:])
    tie_line = hull.structure_slice[hull.hull.vertices]

    # plot hull structures
    if plot_hull_points:
        ax.scatter(tie_line[:, 0], tie_line[:, 1],
                   c=hull.colours[1],
                   marker='o', zorder=99999, edgecolor='k',
                   s=scale*40, lw=1.5)
    ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1],
            c=hull.colours[0], zorder=1)

    if hull.hull_cutoff > 0:
        ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1] + hull.hull_cutoff,
                '--', c=hull.colours[1], alpha=0.5, zorder=1, label='')

    # annotate hull structures
    if labels or label_cutoff is not None:
        label_cursor = get_hull_labels(hull, num_species=2)
        for ind, doc in enumerate(label_cursor):
            arrowprops = dict(arrowstyle="-|>", lw=2, alpha=1, zorder=1, shrinkA=2, shrinkB=4)
            if (ind + 2) < np.argmin(tie_line[:, 1]):
                position = (0.8 * tie_line[ind + 2, 0], 1.15 * (tie_line[ind + 2, 1]) - 0.05)
            elif (ind + 2) == np.argmin(tie_line[:, 1]):
                position = (tie_line[ind + 2, 0], 1.15 * (tie_line[ind + 2, 1]) - 0.05)
            else:
                position = (min(1.1 * tie_line[ind + 2, 0] + 0.15, 0.95), 1.15 * (tie_line[ind + 2, 1]) - 0.05)
            ax.annotate(get_formula_from_stoich(doc['stoichiometry'],
                                                elements=hull.elements,
                                                latex_sub_style=r'\mathregular',
                                                tex=True),
                        xy=(tie_line[ind+2, 0], tie_line[ind+2, 1]),
                        xytext=position,
                        textcoords='data',
                        ha='right',
                        va='bottom',
                        arrowprops=arrowprops,
                        zorder=1)

    # points for off hull structures; we either colour by source or by energy
    if plot_points and not colour_by_source:

        if hull.hull_cutoff == 0:
            # if no specified hull cutoff, ignore labels and colour by hull distance
            cmap = hull.default_cmap
            if plot_points:
                scatter = ax.scatter(hull.structures[np.argsort(hull.hull_dist), 0][::-1],
                                     hull.structures[np.argsort(hull.hull_dist), -1][::-1],
                                     s=scale*40,
                                     c=np.sort(hull.hull_dist)[::-1],
                                     zorder=10000,
                                     cmap=cmap, norm=colours.LogNorm(0.02, 2))

                cbar = plt.colorbar(scatter, aspect=30, pad=0.02,
                                    ticks=[0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                cbar.ax.tick_params(length=0)
                cbar.ax.set_yticklabels([0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                cbar.ax.yaxis.set_ticks_position('right')
                cbar.ax.set_frame_on(False)
                cbar.outline.set_visible(False)
                cbar.set_label('Distance from hull (eV)')

        elif hull.hull_cutoff != 0:
            # if specified hull cutoff colour those below
            c = hull.colours[1]
            for ind in range(len(hull.structures)):
                if hull.hull_dist[ind] <= hull.hull_cutoff or hull.hull_cutoff == 0:
                    if plot_points:
                        scatter.append(ax.scatter(hull.structures[ind, 0], hull.structures[ind, 1],
                                                  s=scale*40,
                                                  alpha=0.9, c=c,
                                                  zorder=300))
            if plot_points:
                ax.scatter(hull.structures[1:-1, 0], hull.structures[1:-1, 1],
                           s=scale*30, lw=0,
                           alpha=0.3, c=hull.colours[-2],
                           edgecolor='k', zorder=10)

    elif colour_by_source:
        from matador.utils.cursor_utils import get_guess_doc_provenance
        if sources is None:
            sources = ['AIRSS', 'GA', 'OQMD', 'SWAPS', 'ICSD']
        if source_labels is None:
            source_labels = sources
        else:
            assert len(source_labels) == len(sources)

        colour_choices = {source: hull.colours[ind + 1] for ind, source in enumerate(sources)}
        colours = []
        concs = []
        energies = []
        zorders = []
        for doc in hull.cursor:
            source = get_guess_doc_provenance(doc['source'])
            if source not in sources:
                # use grey for undesired sources
                colours.append(hull.colours[-2])
                if 'Other' not in sources:
                    sources.append('Other')
                    labels.append('Other')
                    colour_choices['Other'] = hull.colours[-2]
            else:
                colours.append(colour_choices[source])
            zorders.append(sources.index(source))
            concs.append(doc['concentration'])
            energies.append(doc['formation_enthalpy_per_atom'])

        alpha = kwargs.get('alpha')
        if alpha is None:
            alpha = 0.2

        for ind, conc in enumerate(concs):
            if hull.cursor[ind]['hull_distance'] <= 0 + 1e-9 and not plot_hull_points:
                ax.scatter(conc, energies[ind],
                           c=colours[ind], alpha=alpha, s=scale*40,
                           zorder=zorders[ind]+1e5, lw=1.5)
            else:
                ax.scatter(conc, energies[ind],
                           c=colours[ind], alpha=alpha, s=scale*20,
                           zorder=zorders[ind]+100)

        for ind, source in enumerate(sources):
            ax.scatter(1e10, 1e10, c=colour_choices[source], label=source_labels[ind], alpha=alpha, lw=1)

        legend = ax.legend(loc=9, facecolor='w', frameon=True, fancybox=False, shadow=False)
        legend.set_zorder(1e20)

    eform_limits = (np.min(hull.structures[:, 1]), np.max(hull.structures[:, 1]))
    lims = (-0.1 if eform_limits[0] >= 0 else 1.4*eform_limits[0],
            eform_limits[1] if eform_limits[0] >= 0 else 0.1)
    ax.set_ylim(lims)

    if isinstance(title, bool) and title:
        if len(one_minus_x_elem) == 1:
            ax.set_title(x_elem[0] + r'$_\mathrm{x}$' + one_minus_x_elem[0] + r'$_\mathrm{1-x}$')
        if hull._non_binary:
            ax.set_title(r'{d[0]}$_\mathrm{{x}}$({d[1]})$_\mathrm{{1-x}}$'.format(d=hull.chempot_search))
    elif title:
        ax.set_title(title)

    plt.locator_params(nbins=3)
    ax.set_xlabel(r'x in {}$_\mathrm{{x}}${}$_\mathrm{{1-x}}$'.format(x_elem[0], one_minus_x_elem[0]))
    ax.grid(False)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.25, 0.33, 0.5, 0.66, 0.75, 1])
    ax.set_xticklabels(ax.get_xticks())
    # ax.set_yticks(np.arange(, np.min(hull.structure_slice[hull.hull.vertices, 1]) - 0.15, -0.1))
    # ax.set_yticklabels(['{:.1f}'.format(val) for val in ax.get_yticks()])
    ax.set_ylabel('Formation energy (eV/atom)')

    if hull.savefig:
        fname = ''.join(hull.elements) + '_hull'
        exts = ['pdf', 'svg', 'png']
        for ext in exts:
            if hull.args.get(ext):
                plt.savefig('{}.{}'.format(fname, ext),
                            dpi=500, bbox_inches='tight', transparent=True)
                print('Wrote {}.{}'.format(fname, ext))
    elif show:
        plt.show()

    return ax


@plotting_function
def plot_ternary_hull(hull, axis=None, show=False, plot_points=True, hull_cutoff=None, label_cutoff=None, expecting_cbar=True, labels=None):
    """ Plot calculated ternary hull as a 2D projection.

    Parameters:
        hull (matador.hull.QueryConvexHull): matador hull object.

    Keyword arguments:
        axis (matplotlib.axes.Axes): matplotlib axis object on which to plot.
        show (bool): whether or not to show plot in X window.
        plot_points (bool): whether or not to plot each structure as a point.
        label_cutoff (float/:obj:`tuple` of :obj:`float`): draw labels less than or
            between these distances form the hull, also read from hull.args.
        expecting_cbar (bool): whether or not to space out the plot to preserve
            aspect ratio if a colourbar is present.
        labels (bool): whether or not to label on-hull structures

    Returns:
        matplotlib.axes.Axes: matplotlib axis with plot.

    """
    import ternary
    import matplotlib.pyplot as plt
    import matplotlib.colors as colours
    from matador.utils.chem_utils import get_generic_grav_capacity

    plt.rcParams['axes.linewidth'] = 0
    plt.rcParams['xtick.major.size'] = 0
    plt.rcParams['ytick.major.size'] = 0
    plt.rcParams['xtick.minor.size'] = 0
    plt.rcParams['ytick.minor.size'] = 0

    if labels is None:
        labels = hull.args.get('labels')
    if label_cutoff is None:
        label_cutoff = hull.args.get('label_cutoff')
        if label_cutoff is None:
            label_cutoff = 0
    else:
        labels = True

    if hull_cutoff is None and hull.hull_cutoff is None:
        hull_cutoff = 0
    else:
        hull_cutoff = hull.hull_cutoff

    print('Plotting ternary hull...')
    if hull.args.get('capmap') or hull.args.get('efmap'):
        scale = 100
    elif hull.args.get('sampmap'):
        scale = 20
    else:
        scale = 1
    fontsize = plt.rcParams['font.size']

    if axis is not None:
        fig, ax = ternary.figure(scale=scale, ax=axis)
    else:
        fig, ax = ternary.figure(scale=scale)
    if hull.args.get('capmap') or hull.args.get('efmap') or hull.args.get('sampmap'):
        fig.set_size_inches(8, 5)
    elif not expecting_cbar:
        fig.set_size_inches(5, 5)
    else:
        fig.set_size_inches(6.67, 5)
    ax.boundary(linewidth=2.0, zorder=99)
    ax.gridlines(color='black', multiple=scale * 0.1, linewidth=0.5)

    ax.clear_matplotlib_ticks()
    ticks = [float(val) for val in np.linspace(0.0, 1.0, 6)]
    ax.ticks(axis='lbr', linewidth=1, multiple=scale*0.2, offset=0.02, fontsize=fontsize-2,
             ticks=ticks, tick_formats='%.1f')

    ax.set_title('-'.join(hull.elements), fontsize=fontsize + 2, y=1.02)
    ax.left_axis_label(hull.elements[2], fontsize=fontsize + 2)
    ax.right_axis_label(hull.elements[1], fontsize=fontsize + 2)
    ax.bottom_axis_label(hull.elements[0], fontsize=fontsize + 2)

    concs = np.zeros((len(hull.structures), 3))

    concs[:, :-1] = hull.structures[:, :-1]
    for i in range(len(concs)):
        # set third triangular coordinate
        concs[i, -1] = 1 - concs[i, 0] - concs[i, 1]

    stable = np.asarray([concs[ind] for ind in hull.hull.vertices])

    # sort by hull distances so things are plotting the right order
    concs = concs[np.argsort(hull.hull_dist)].tolist()
    hull_dist = np.sort(hull.hull_dist)

    filtered_concs = []
    filtered_hull_dists = []
    for ind, conc in enumerate(concs):
        if conc not in filtered_concs:
            if hull_dist[ind] <= hull.hull_cutoff or (hull.hull_cutoff == 0 and hull_dist[ind] < 0.1):
                filtered_concs.append(conc)
                filtered_hull_dists.append(hull_dist[ind])
    if hull.args.get('debug'):
        print('Trying to plot {} points...'.format(len(filtered_concs)))

    concs = np.asarray(filtered_concs)
    hull_dist = np.asarray(filtered_hull_dists)

    min_cut = 0.0
    max_cut = 0.2

    hull.colours = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    hull.default_cmap_list = get_linear_cmap(hull.colours[1:4], list_only=True)
    hull.default_cmap = get_linear_cmap(hull.colours[1:4], list_only=False)
    n_colours = len(hull.default_cmap_list)
    colours_hull = hull.default_cmap_list

    cmap = hull.default_cmap
    cmap_full = plt.cm.get_cmap('Pastel2')
    pastel_cmap = colours.LinearSegmentedColormap.from_list('Pastel2', cmap_full.colors)

    for plane in hull.hull.planes:
        plane.append(plane[0])
        plane = np.asarray(plane)
        ax.plot(scale * plane, c=hull.colours[0], lw=1.5, alpha=1, zorder=98)

    if hull.args.get('pathways'):
        for phase in stable:
            if phase[0] == 0 and phase[1] != 0 and phase[2] != 0:
                ax.plot([scale * phase, [scale, 0, 0]], c='r', alpha=0.2, lw=6, zorder=99)

    # add points
    if plot_points:
        colours_list = []
        colour_metric = hull_dist
        for i, _ in enumerate(colour_metric):
            if colour_metric[i] >= max_cut:
                colours_list.append(n_colours - 1)
            elif colour_metric[i] <= min_cut:
                colours_list.append(0)
            else:
                colours_list.append(int((n_colours - 1) * (colour_metric[i] / max_cut)))
        colours_list = np.asarray(colours_list)
        ax.scatter(scale*stable, marker='o', color=hull.colours[1], edgecolors='black', zorder=9999999,
                   s=150, lw=1.5)
        ax.scatter(scale*concs, colormap=cmap, colorbar=True, cbarlabel='Distance from hull (eV/atom)',
                   c=hull_dist, vmax=max_cut, vmin=min_cut, zorder=1000, s=40, alpha=0)
        for i, _ in enumerate(concs):
            ax.scatter(
                scale * concs[i].reshape(1, 3),
                color=colours_hull[colours_list[i]],
                marker='o',
                zorder=10000 - colours_list[i],
                s=70 * (1 - float(colours_list[i]) / n_colours) + 15,
                lw=1,
                edgecolors='black'
            )

    # add colourmaps
    if hull.args.get('capmap'):
        capacities = dict()
        from ternary.helpers import simplex_iterator
        for (i, j, k) in simplex_iterator(scale):
            capacities[(i, j, k)] = get_generic_grav_capacity([
                float(i) / scale, float(j) / scale, float(scale - i - j) / scale
            ], hull.elements)
        ax.heatmap(capacities, style="hexagonal", cbarlabel='Gravimetric capacity (maH/g)',
                   vmin=0, vmax=3000, cmap=pastel_cmap)
    elif hull.args.get('efmap'):
        energies = dict()
        fake_structures = []
        from ternary.helpers import simplex_iterator
        for (i, j, k) in simplex_iterator(scale):
            fake_structures.append([float(i) / scale, float(j) / scale, 0.0])
        fake_structures = np.asarray(fake_structures)
        plane_energies, _, _ = hull.get_hull_distances(fake_structures)
        ind = 0
        for (i, j, k) in simplex_iterator(scale):
            energies[(i, j, k)] = -1 * plane_energies[ind]
            ind += 1
        ax.heatmap(energies, style="hexagonal", cbarlabel='Formation energy (eV/atom)', vmax=0, cmap='bone')
    elif hull.args.get('sampmap'):
        sampling = dict()
        from ternary.helpers import simplex_iterator
        eps = 1.0 / float(scale)
        for (i, j, k) in simplex_iterator(scale):
            sampling[(i, j, k)] = np.size(np.where((concs[:, 0] <= float(i)/scale + eps) *
                                                   (concs[:, 0] >= float(i)/scale - eps) *
                                                   (concs[:, 1] <= float(j)/scale + eps) *
                                                   (concs[:, 1] >= float(j)/scale - eps) *
                                                   (concs[:, 2] <= float(k)/scale + eps) *
                                                   (concs[:, 2] >= float(k)/scale - eps)))
        ax.heatmap(sampling, style="hexagonal", cbarlabel='Number of structures', cmap='afmhot')

    # add labels
    if labels:
        label_cursor = get_hull_labels(hull, label_cutoff=label_cutoff, num_species=3)
        if len(label_cursor) == 1:
            label_coords = [[0.25, 0.5]]
        else:
            label_coords = [[0.1+(val-0.5)*0.3, val] for val in np.linspace(0.5, 0.8, int(round(len(label_cursor)/2.)+1))]
            label_coords += [[0.9-(val-0.5)*0.3, val+0.2] for val in np.linspace(0.5, 0.8, int(round(len(label_cursor)/2.)))]
        from matador.utils.hull_utils import barycentric2cart
        for ind, doc in enumerate(label_cursor):
            conc = np.asarray(doc['concentration'] + [1 - sum(doc['concentration'])])
            formula = get_formula_from_stoich(doc['stoichiometry'], tex=True, elements=hull.elements, latex_sub_style=r'\mathregular')
            arrowprops = dict(arrowstyle="-|>", color='k', lw=2, alpha=0.5, zorder=1, shrinkA=2, shrinkB=4)
            cart = barycentric2cart([doc['concentration'] + [0]])[0][:2]
            min_dist = 1e20
            closest_label = 0
            for coord_ind, coord in enumerate(label_coords):
                dist = np.sqrt((cart[0] - coord[0])**2 + (cart[1] - coord[1])**2)
                if dist < min_dist:
                    min_dist = dist
                    closest_label = coord_ind
            ax.annotate(formula, scale*conc,
                        textcoords='data',
                        xytext=[scale*val for val in label_coords[closest_label]],
                        ha='right',
                        va='bottom',
                        arrowprops=arrowprops)
            del label_coords[closest_label]

    plt.tight_layout(w_pad=0.2)

    if hull.savefig:
        if hull.args.get('png'):
            plt.savefig(''.join(hull.elements) + '_hull.png', dpi=400, transparent=True, bbox_inches='tight')
        if hull.args.get('svg'):
            plt.savefig(''.join(hull.elements) + '_hull.svg', dpi=400, transparent=True, bbox_inches='tight')
        if hull.args.get('pdf'):
            plt.savefig(''.join(hull.elements) + '_hull.pdf', dpi=400, transparent=True, bbox_inches='tight')
    elif show:
        ax.show()

    return ax
