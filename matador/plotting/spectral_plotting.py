# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot densities of states and
bandstructures for electronic and vibrational calculations.

"""


import os
import numpy as np
from matador.utils.viz_utils import get_element_colours
from matador.plotting.plotting import plotting_function
from matador.scrapers import optados2dict, phonon2dict, bands2dict
from matador.scrapers import cell2dict, res2dict
from matador.orm.spectral import ElectronicDispersion, VibrationalDispersion


@plotting_function
def plot_spectral(seeds, **kwargs):
    """ This function wraps all of the spectral plotting capability of matador.
    When provided with a seed, or seeds, several files will be checked:

        - <seed>.bands: CASTEP bandstructure expected, not DOS,
        - <seed>.adaptive.dat: OptaDOS total DOS,
        - <seed>.pdos.dat: OptaDOS pDOS,
        - <seed>.pdis.dat: OptaDOS projected dispersion,

    or, if the "phonons" flag is passed, or if a <seed>.phonon file is detected,

        - <seed>.phonon: CASTEP phonon dispersion curves expected,
        - <seed>.phonon_dos: CASTEP phonon DOS.

    This function will then automatically choose what to plot, favouring a bandstructure
    with "the-most-projected" DOS it can find.

    Parameters:
        seeds (list): list of filenames of bands/phonon files

    Keyword Arguments:
        plot_bandstructure (bool): whether to plot bandstructure, if available
        plot_dos (bool): whether to plot density of states, if available
        plot_pdos (bool): whether or not to plot projected DOS, if available
        plot_pdis (bool): whether to plot projected dispersion, if available
        dos (str): separate seed name for pDOS/DOS data
        phonons (bool): whether to plot phonon or electronic data
        labels (list): list of strings for legend labels for multiple bandstructures
        gap (bool): whether to draw on the band gap
        colour_by_seed (bool): plot with a separate colour per bandstructure
        external_efermi (float or list): replace scraped Fermi energy with this value (eV) (can be
            specified per spin channel).
        highlight_bands (list): list of integer indices, colour the bands with these indices in red
        band_colour (str): if passed "occ", bands will be coloured using cmap depending on whether
            they lie above or below the Fermi level. If passed 'random', colour bands randomly from
            the cmap. Otherwise, override all colour options with matplotlib-interpretable colour
            (e.g. hexcode or html colour name) to use for all bands (DEFAULT: 'occ').
        cmap (str): matplotlib colourmap name to use for the bands
        n_colours (int): number of colours to use from cmap (DEFAULT: 4).
        no_stacked_pdos (bool): whether to plot projected DOS as stack or overlapping.
        spin_only (str): either 'up' or 'down' to only plot one spin channel.
        preserve_kspace_distance (bool): whether to preserve distances in reciprocal space when
            linearising the kpoint path. If False, bandstructures of different lattice parameters
            with the same Bravais lattice can be more easily compared. If True, bandstructures may
            appear rarefied or compressed in particular regions.
        pdis_interpolation_factor (float): multiple by which to interpolate pDIS bands
        pdis_point_scale (float): size of points in pDIS (DEFAULT: 50).
        band_reorder (bool): try to reorder bands based on local gradients (DEFAULT: True for phonons, otherwise False).
        title (str): optional plot title
        pdos_hide_tot (bool): whether or not to plot the total DOS on a PDOS plot; this is to hide
            regions where the PDOS is negative (leading to total DOS lower than stacked PDOS) (DEFAULT: False).

    """
    import matplotlib.pyplot as plt
    from cycler import cycler
    # set defaults and update class with desired values
    prop_defaults = {'plot_bandstructure': True, 'plot_dos': True, 'plot_pdos': True, 'plot_pdis': True,
                     'phonons': False, 'gap': False,
                     'colour_by_seed': False, 'external_efermi': None,
                     'labels': None, 'cmap': None, 'band_colour': 'occ',
                     'n_colours': 4, 'spin_only': None, 'figsize': None,
                     'pdis_interpolation_factor': 2, 'pdis_point_scale': 25,
                     'no_stacked_pdos': False, 'preserve_kspace_distance': False,
                     'band_reorder': None, 'title': None,
                     'verbosity': 0, 'highlight_bands': None, 'pdos_hide_tot': True}
    for key in kwargs:
        if kwargs[key] is not None:
            prop_defaults[key] = kwargs[key]
    kwargs = prop_defaults

    if kwargs.get('cmap') is None:
        kwargs['colours'] = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        plt.rcParams['axes.prop_cycle'] = cycler('color', kwargs['colours'])
    else:
        if isinstance(kwargs['cmap'], str):
            print('Adjusting colour palette... to {}'.format(kwargs.get('cmap')))
            kwargs['colours'] = plt.cm.get_cmap(kwargs.get('cmap')).colors
            plt.rcParams['axes.prop_cycle'] = cycler('color', kwargs['colours'])
        elif isinstance(kwargs['cmap'], list):
            print('Reading list of colours {}...'.format(kwargs.get('cmap')))
            kwargs['colours'] = kwargs['cmap']
            plt.rcParams['axes.prop_cycle'] = cycler('color', kwargs['colours'])

    if (kwargs.get('phonons') and kwargs['band_colour'] == 'occ') or kwargs['band_colour'] == 'random':
        kwargs['band_colour'] = None

    if not isinstance(seeds, list):
        seeds = [seeds]

    if len(seeds) > 1:
        if kwargs['plot_pdis'] or kwargs['plot_dos']:
            kwargs['plot_pdos'] = False
            kwargs['plot_pdis'] = False
            print('Disabling projections as mutiple seeds requested.')

    if kwargs.get('plot_window') is not None:
        if isinstance(kwargs.get('plot_window'), list):
            if len(kwargs.get('plot_window')) != 2:
                exit('plot_window must have length 2 or be a single number')
            kwargs['plot_window'] = sorted(kwargs.get('plot_window'))
        else:
            kwargs['plot_window'] = (-kwargs.get('plot_window'), kwargs.get('plot_window'))
    else:
        kwargs['plot_window'] = None

    if kwargs['plot_dos']:
        # check an optados file exists
        exts = ['pdos.dat', 'adaptive.dat', 'fixed.dat', 'linear.dat', 'jdos.dat', 'phonon_dos', 'bands_dos']
        kwargs['plot_dos'] = any([any([os.path.isfile('{}.{}'.format(seed, ext)) for ext in exts]) for seed in seeds])

    if kwargs['plot_pdos']:
        exts = ['pdos.dat', 'phonon_dos']
        kwargs['plot_pdos'] = any([any([os.path.isfile('{}.{}'.format(seed, ext)) for ext in exts]) for seed in seeds])

    figsize = kwargs['figsize']
    if kwargs['plot_bandstructure'] and not kwargs['plot_dos']:
        if figsize is None:
            figsize = (7, 6)
        fig, ax_dispersion = plt.subplots(figsize=figsize)
    elif kwargs['plot_bandstructure'] and kwargs['plot_dos']:
        if figsize is None:
            figsize = (10, 6)
        fig, ax_grid = plt.subplots(1, 3, figsize=figsize, sharey=True,
                                    gridspec_kw={'width_ratios': [4, 2, 1],
                                                 'wspace': 0.1,
                                                 'left': 0.15})
        ax_dispersion = ax_grid[0]
        ax_dos = ax_grid[1]
        ax_grid[2].axis('off')
    elif not kwargs['plot_bandstructure'] and kwargs['plot_dos']:
        if figsize is None:
            figsize = (6, 3)
        fig, ax_dos = plt.subplots(1, figsize=figsize)

    kwargs['valence'] = kwargs['colours'][0]
    kwargs['conduction'] = kwargs['colours'][-1]
    kwargs['crossing'] = kwargs['colours'][int(len(kwargs['colours']) / 2)]

    if len(seeds) > 1 or kwargs.get('colour_by_seed'):
        kwargs['seed_colours'] = kwargs['colours']
        kwargs['ls'] = ['-'] * len(seeds)
        kwargs['colour_by_seed'] = True
        if kwargs.get('labels') is None:
            kwargs['labels'] = [seed.split('/')[-1].split('.')[0] for seed in seeds]

    kwargs['ls'] = []
    for i in range(len(seeds)):
        if i % 3 == 0:
            kwargs['ls'].append('-')
        elif i % 3 == 1:
            kwargs['ls'].append('--')
        elif i % 3 == 2:
            kwargs['ls'].append('-.')

    bbox_extra_artists = []
    if kwargs['plot_bandstructure']:
        ax_dispersion = dispersion_plot(seeds, ax_dispersion, kwargs, bbox_extra_artists)

    if kwargs['plot_dos']:
        ax_dos = dos_plot(seeds, ax_dos, kwargs, bbox_extra_artists)

    if kwargs.get('title') is not None:
        fig.suptitle(kwargs.get('title'))

    if any([kwargs.get('pdf'), kwargs.get('svg'), kwargs.get('png')]):
        if not bbox_extra_artists:
            bbox_extra_artists = None
        filename = seeds[0].split('/')[-1].replace('.bands', '').replace('.phonon', '') + '_spectral'
        if kwargs.get('pdf'):
            plt.savefig('{}.pdf'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
        if kwargs.get('svg'):
            plt.savefig('{}.svg'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)
        if kwargs.get('png'):
            plt.savefig('{}.png'.format(filename),
                        bbox_inches='tight', transparent=True, bbox_extra_artists=bbox_extra_artists)

    else:
        plt.tight_layout()
        print('Displaying plot...')
        plt.show()


def dispersion_plot(seeds, ax_dispersion, kwargs, bbox_extra_artists):
    """ Plot a dispersion/bandstructure on the given axis. Will detect
    and projected dispersion data automatically.

    Parameters:
        seeds (list): the seednames of the data to plot.
        ax_dispersion (matplotlib.Axes): the axis to plot on.
        kwargs (dict): any plotting keywords (from e.g. dispersion script).
        bbox_extra_artists (list): a list to which to append legends.

    Returns:
        matplotlib.Axes: the axis that was plotted on.

    """
    from cycler import cycler
    import matplotlib.pyplot as plt
    plotted_pdis = False
    for seed_ind, seed in enumerate(seeds):
        seed = seed.replace('.bands', '').replace('.phonon', '')
        if os.path.isfile('{}.phonon'.format(seed)):
            dispersion, s = phonon2dict(seed + '.phonon', verbosity=kwargs.get('verbosity'))
            if not s:
                raise RuntimeError(dispersion)

            dispersion = VibrationalDispersion(dispersion)

            if kwargs['plot_window'] is None:
                kwargs['plot_window'] = [min(-10, np.min(dispersion.eigs) - 10), np.max(dispersion.eigs)]

        elif os.path.isfile('{}.bands'.format(seed)):
            dispersion, s = bands2dict(seed + '.bands',
                                       verbosity=kwargs.get('verbosity'))
            if not s:
                raise RuntimeError(dispersion)

            if os.path.isfile('{}.pdis.dat'.format(seed)) and kwargs['plot_pdis']:
                pdis_data, s = optados2dict('{}.pdis.dat'.format(seed))
                if not s:
                    raise RuntimeError(pdis_data)
            else:
                pdis_data = None

            dispersion = ElectronicDispersion(dispersion, pdis_data)

            if kwargs['plot_window'] is None:
                kwargs['plot_window'] = [-10, 10]

        else:
            raise RuntimeError('{}.bands/.phonon not found.'.format(seed))

        path = dispersion.linearise_path(preserve_kspace_distance=kwargs['preserve_kspace_distance'])

        # try to match bands if requested
        if kwargs['band_reorder'] or (kwargs['band_reorder'] is None and kwargs['phonons']):
            print('Reordering bands based on local gradients...')
            dispersion.reorder_bands()

        if dispersion.projectors and len(seeds) == 1 and kwargs['plot_pdis']:
            ax_dispersion = projected_bandstructure_plot(dispersion, ax_dispersion, path,
                                                         bbox_extra_artists,
                                                         **kwargs)
            kwargs['band_colour'] = 'grey'
            plotted_pdis = True

        # loop over branches and plot
        if not plotted_pdis:
            if kwargs.get('external_efermi') is None:
                spin_fermi_energy = dispersion.spin_fermi_energy
            else:
                spin_fermi_energy = kwargs.get('external_efermi')
            if len(spin_fermi_energy) == 1 and dispersion.num_spins != 1:
                spin_fermi_energy = [spin_fermi_energy] * dispersion.num_spins

            for branch_ind, branch in enumerate(dispersion.kpoint_branches):

                # seem to have to reset colours here for some reason
                plt.rcParams['axes.prop_cycle'] = cycler('color', kwargs['colours'])

                for ns in range(dispersion.num_spins):

                    if ns == 1 and kwargs.get('spin_only') == 'up':
                        continue
                    elif ns == 0 and kwargs.get('spin_only') == 'down':
                        continue

                    for nb in range(dispersion.num_bands):
                        colour, alpha, label = _get_lineprops(dispersion, spin_fermi_energy, nb, ns, branch, branch_ind, seed_ind, kwargs)

                        ax_dispersion.plot(path[(np.asarray(branch)-branch_ind).tolist()],
                                           dispersion.eigs[ns][nb][branch] - spin_fermi_energy[ns],
                                           c=colour, ls=kwargs['ls'][seed_ind], alpha=alpha, label=label)

    if len(seeds) > 1:
        disp_legend = ax_dispersion.legend(loc='upper center',
                                           frameon=True, fancybox=False, shadow=False, framealpha=1)
        bbox_extra_artists.append(disp_legend)

    ax_dispersion.axhline(0, ls='--', lw=1, c='grey')
    ax_dispersion.set_ylim(kwargs['plot_window'])
    if kwargs['phonons']:
        ylabel = 'Wavenumber (cm$^{-1}$)'
    else:
        ylabel = r'Energy (eV)'
    ax_dispersion.set_ylabel(ylabel)
    ax_dispersion.set_xlim(0, 1)
    _add_path_labels(seed, dispersion, ax_dispersion, path, seed_ind, kwargs)

    return ax_dispersion


def dos_plot(seeds, ax_dos, kwargs, bbox_extra_artists):
    """ Plot a density of states on the given axis. Will detect
    pDOS and spin-dependent DOS data automatically.

    Parameters:
        seeds (list): the seednames of the data to plot.
        ax_dos (matplotlib.Axes): the axis to plot on.
        kwargs (dict): any plotting keywords (from e.g. dispersion script).
        bbox_extra_artists (list): a list to which to append legends.

    Returns:
        matplotlib.Axes: the axis that was plotted on.

    """
    for seed_ind, seed in enumerate(seeds):
        seed = seed.replace('.bands', '').replace('.phonon', '')
        if kwargs['plot_dos']:
            if not kwargs['phonons']:
                if kwargs.get('dos') is None:
                    # look for dat files, and just use the first
                    exts = ['adaptive.dat', 'fixed.dat', 'linear.dat', 'bands_dos']
                    for ext in exts:
                        if os.path.isfile('{}.{}'.format(seed, ext)):
                            dos_seed = '{}.{}'.format(seed, ext)
                            break
                    else:
                        raise SystemExit('No total DOS files found.')
                else:
                    dos_seed = kwargs.get('dos')

                if dos_seed.endswith('.bands_dos'):
                    # if bands_dos exists, do some manual broadening
                    dos_data, s = bands2dict(dos_seed)
                    raw_weights = []
                    space_size = 1000
                    gaussian_width = kwargs.get('gaussian_width')
                    if gaussian_width is None:
                        gaussian_width = 0.1

                    raw_eigenvalues = []
                    for kind, qpt in enumerate(dos_data['eigenvalues_k_s']):
                        weight = dos_data['kpoint_weights'][kind]
                        for eig in qpt:
                            raw_weights.append(weight)
                            raw_eigenvalues.append(eig)
                    raw_eigenvalues = np.asarray(raw_eigenvalues)
                    hist, energies = np.histogram(raw_eigenvalues, bins=space_size)
                    # shift bin edges to bin centres
                    energies -= energies[1] - energies[0]
                    energies = energies[:-1]
                    new_energies = np.reshape(energies, (1, len(energies)))
                    new_energies = new_energies - np.reshape(energies, (1, len(energies))).T
                    dos = np.sum(hist * np.exp(-(new_energies)**2 / gaussian_width), axis=1)
                    dos = np.divide(dos, np.sqrt(2 * np.pi * gaussian_width**2))
                    dos_data['dos'] = dos
                    dos_data['energies'] = energies
                else:
                    dos_data, s = optados2dict(dos_seed, verbosity=0)

                if not s:
                    raise RuntimeError(dos_data)

                energies = dos_data['energies']
                dos = dos_data['dos']

                if kwargs['plot_window'] is None:
                    kwargs['plot_window'] = [-10, 10]

                if 'spin_dos' in dos_data:
                    max_density = max(np.max(np.abs(dos_data['spin_dos']['down'][np.where(energies > kwargs['plot_window'][0])])),
                                      np.max(np.abs(dos_data['spin_dos']['up'][np.where(energies > kwargs['plot_window'][0])])))
                else:
                    max_density = np.max(dos[np.where(np.logical_and(energies < kwargs['plot_window'][1], energies > kwargs['plot_window'][0]))])

                if kwargs['plot_pdos']:
                    pdos_seed = '{}.pdos.dat'.format(seed)
                    pdos_data = {}
                    if os.path.isfile(pdos_seed):
                        pdos_data, s = optados2dict(pdos_seed, verbosity=0)
                        if not s:
                            raise RuntimeError(pdos_data)
                        dos_data['pdos'] = pdos_data
            else:
                if not os.path.isfile(seed + '.phonon_dos'):
                    phonon_data, s = phonon2dict(seed + '.phonon')
                    if not s:
                        raise RuntimeError(phonon_data)
                    else:
                        if kwargs['plot_window'] is None:
                            kwargs['plot_window'] = [min(-10, np.min(phonon_data['eigenvalues_q']) - 10), np.max(phonon_data['eigenvalues_q'])]
                        space_size = 1000
                        gaussian_width = kwargs.get('gaussian_width', 100)
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
                    if kwargs['plot_window'] is None:
                        kwargs['plot_window'] = [np.min(energies[np.where(dos > 1e-3)]) - 10, np.max(energies[np.where(dos > 1e-3)])]

                    if kwargs['plot_pdos']:
                        dos_data['pdos'] = dict()
                        for i, label in enumerate(projector_labels):
                            dos_data['pdos'][label] = raw_data[:, i + 2]
                        pdos = dos_data['pdos']
                        pdos_data = dos_data

                    from os import remove
                    remove(seed + '.phonon_dos_tmp')

            if kwargs['phonons']:
                ylabel = 'Phonon DOS'
                xlabel = 'Wavenumber (cm$^{{-1}}$)'
            else:
                if 'dos_unit_label' in dos_data:
                    ylabel = dos_data['dos_unit_label'].replace('A^3', 'Å$^{3}$')
                else:
                    if kwargs['plot_bandstructure']:
                        ylabel = 'DOS'
                    else:
                        ylabel = 'DOS (eV$^{{-1}}$Å$^{{-3}}$)'
                xlabel = 'Energy (eV)'

            if kwargs['plot_bandstructure']:
                ax_dos.set_xlabel(ylabel)
                ax_dos.axhline(0, c='grey', ls='--', lw=1)
                if 'spin_dos' in dos_data:
                    ax_dos.set_xlim(-max_density*1.2, max_density * 1.2)
                else:
                    ax_dos.set_xlim(0, max_density * 1.2)
                ax_dos.set_ylim(kwargs['plot_window'])
                ax_dos.axvline(0, c='grey', lw=1)
                ax_dos.xaxis.set_ticks_position('none')

                if 'spin_dos' not in dos_data:
                    ax_dos.plot(dos, energies, ls=kwargs['ls'][seed_ind],
                                color='grey', zorder=1e10, label='Total DOS')
                    if not kwargs['plot_dos']:
                        ax_dos.fill_betweenx(energies[np.where(energies > 0)], 0, dos[np.where(energies > 0)], alpha=0.2, color=kwargs['conduction'])
                        ax_dos.fill_betweenx(energies[np.where(energies <= 0)], 0, dos[np.where(energies <= 0)], alpha=0.2, color=kwargs['valence'])

            else:
                ax_dos.set_xlabel(xlabel)
                ax_dos.set_ylabel(ylabel)
                ax_dos.axvline(0, c='grey', lw=1, ls='--')
                if 'spin_dos' in dos_data:
                    ax_dos.set_ylim(-max_density*1.2, max_density * 1.2)
                else:
                    ax_dos.set_ylim(0, max_density * 1.2)
                ax_dos.set_xlim(kwargs['plot_window'])
                ax_dos.axhline(0, c='grey', lw=1)

                if 'spin_dos' not in dos_data:
                    ax_dos.plot(energies, dos, ls=kwargs['ls'][seed_ind], alpha=1,
                                c='grey', zorder=1e10, label='Total DOS')
                    if not kwargs['plot_dos']:
                        ax_dos.fill_between(energies[np.where(energies > 0)], 0, dos[np.where(energies > 0)], alpha=0.2, color=kwargs['conduction'])
                        ax_dos.fill_between(energies[np.where(energies <= 0)], 0, dos[np.where(energies <= 0)], alpha=0.2, color=kwargs['valence'])

            if kwargs['plot_pdos'] and len(seeds) == 1:
                pdos = pdos_data['pdos']
                energies = pdos_data['energies']
                stack = np.zeros_like(pdos[list(pdos.keys())[0]])
                projector_labels, dos_colours = _get_projector_info([projector for projector in pdos])
                for ind, projector in enumerate(pdos):

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
                            label = None
                            if not kwargs['no_stacked_pdos']:
                                ax_dos.fill_betweenx(energies, stack, stack+pdos[projector],
                                                     alpha=alpha, label=projector_labels[ind],
                                                     color=dos_colours[ind])
                            else:
                                label = projector_labels[ind]
                            ax_dos.plot(stack + pdos[projector], energies,
                                        alpha=1, color=dos_colours[ind], label=label)
                        else:
                            label = None
                            if not kwargs['no_stacked_pdos']:
                                ax_dos.fill_between(energies, stack, stack+pdos[projector],
                                                    alpha=alpha, label=projector_labels[ind],
                                                    color=dos_colours[ind])
                            else:
                                label = projector_labels[ind]
                            ax_dos.plot(energies, stack + pdos[projector],
                                        alpha=1, color=dos_colours[ind], label=label)

                        if not kwargs['no_stacked_pdos']:
                            stack += pdos[projector]

            elif 'spin_dos' in dos_data:
                if kwargs['plot_bandstructure']:
                    if kwargs.get('spin_only') in [None, 'down']:
                        print('Plotting only spin down channel...')
                        ax_dos.fill_betweenx(energies, 0, dos_data['spin_dos']['down'], alpha=0.2, color='b')
                        ax_dos.plot(dos_data['spin_dos']['down'], energies, ls=kwargs['ls'][seed_ind], color='b', zorder=1e10, label='spin-down channel')
                    if kwargs.get('spin_only') in [None, 'up']:
                        print('Plotting only spin up channel...')
                        ax_dos.fill_betweenx(energies, 0, dos_data['spin_dos']['up'], alpha=0.2, color='r')
                        ax_dos.plot(dos_data['spin_dos']['up'], energies, ls=kwargs['ls'][seed_ind], color='r', zorder=1e10, label='spin-up channel')
                else:
                    if kwargs.get('spin_only') in [None, 'down']:
                        print('Plotting only spin down channel...')
                        ax_dos.plot(energies, dos_data['spin_dos']['down'], ls=kwargs['ls'][seed_ind], color='b', zorder=1e10, label='spin-down channel')
                        ax_dos.fill_between(energies, 0, dos_data['spin_dos']['down'], alpha=0.2, color='b')
                    if kwargs.get('spin_only') in [None, 'up']:
                        print('Plotting only spin up channel...')
                        ax_dos.plot(energies, dos_data['spin_dos']['up'], ls=kwargs['ls'][seed_ind], color='r', zorder=1e10, label='spin-up channel')
                        ax_dos.fill_between(energies, 0, dos_data['spin_dos']['up'], alpha=0.2, color='r')

            if len(seeds) == 1:
                if kwargs['plot_bandstructure']:
                    dos_legend = ax_dos.legend(bbox_to_anchor=(1, 1),
                                               frameon=True, fancybox=False, shadow=False)
                else:
                    dos_legend = ax_dos.legend(frameon=True, fancybox=False, shadow=False)
                bbox_extra_artists.append(dos_legend)

    return ax_dos


def projected_bandstructure_plot(dispersion, ax, path, bbox_extra_artists, interpolation_factor=2, point_scale=25, **kwargs):
    """ Plot projected bandstructure with weightings from OptaDOS pdis.dat file.

    Parameters:
        dispersion (matador.orm.spectral.ElectronicDispersion): scraped
            data for bandstructure and pdis.
        seed (str): seed name for files to scrape.
        ax (matplotlib.pyplot.Axes): axis to plot on.
        bbox_extra_artists (list): list to append any legends too.

    Keyword arguments:
        interpolation_factor (float): amount by which to interpolate bands.
        point_scale (float): rescale points by this amount

    Returns:
        matplotlib.pyplot.Axes: the axis that was plotted on.

    """

    pdis = dispersion.projector_weights
    projectors = dispersion.projectors

    pdis[pdis < 0] = 0
    pdis[pdis > 1] = 1

    keep_inds = []
    for ind, _ in enumerate(projectors):
        if np.max(pdis[:, :, ind]) > 1e-8:
            keep_inds.append(ind)

    projector_labels, dos_colours = _get_projector_info(projectors)

    fermi_energy = kwargs.get('external_efermi') or dispersion.fermi_energy

    _ordered_scatter(path, dispersion.eigs_s_k[0].T - fermi_energy, pdis, dispersion.kpoint_branches,
                     interpolation_factor=interpolation_factor, point_scale=point_scale,
                     ax=ax, colours=dos_colours)

    if not kwargs['plot_pdos']:
        for ind, _ in enumerate(projectors):
            if ind in keep_inds:
                ax.scatter(1e20, 0, facecolor=dos_colours[ind],
                           label=projector_labels[ind], lw=0)

        legend = ax.legend(loc=1, frameon=True, fancybox=False, shadow=False)
        legend.set_zorder(1e20)
        bbox_extra_artists.append(legend)

    return ax


def _ordered_scatter(path, eigs, pdis, branches, ax=None, colours=None, interpolation_factor=2, point_scale=25):
    """ Plots an ordered scatter plot of a projected bandstructure.

    Parameters:
        path (np.ndarray): linearised [0, 1] kpoint path array.
        eigs (np.ndarray): (num_kpoints x num_bands) array containing eigenvalues
        pdis (np.ndarray): (num_kpoints x num_bands x num_projectors) array containing
            projector weights.
        branches (list): list of branch indices, e.g. for two branches [[0,1,2], [3, 4]].

    Keyword arguments:
        ax (matplotlib.Axes): axis to plot on
        colours (list): colours assigned for each projector.
        interpolation_factor (float): multiplier for fineness of band interpolation.

    """
    from scipy.interpolate import interp1d
    flat_pts_k = []
    flat_pts_e = []
    flat_sizes = []
    flat_colours = []
    flat_zorders = []

    for nb in range(len(eigs[0])):
        for branch_ind, branch in enumerate(branches):
            k = path[(np.asarray(branch) - branch_ind).tolist()]
            e = eigs[branch, nb]
            projections = pdis[branch, nb]

            ek_fn = interp1d(k, e)
            k_interp = np.linspace(np.min(k), np.max(k), num=int(interpolation_factor*len(k)))
            ek_interp = ek_fn(k_interp)
            projections = projections.T
            interp_projections = []
            for i, _ in enumerate(projections):
                interp_projections.append(interp1d(k, projections[i])(k_interp))
            projections = np.asarray(interp_projections).T
            pts = np.array([k_interp, ek_interp]).T.reshape(-1, 1, 2)

            if colours is not None:
                plot_colours = [colours[i] for i in range(len(projections[0]))]
            else:
                plot_colours = [None for i in range(len(projections[0]))]
            for i, _ in enumerate(projections):
                # zeros mess up zorder, so add small shift then subtract
                # before calculating real sizes
                projections[projections <= 1e-9] = 1e-9
                sizes = np.cumsum(projections[i])
                zorders = 1000*(-nb + 1-sizes)
                projections[projections <= 1e-9] = 0
                sizes = np.cumsum(projections[i])
                for j in range(len(projections[i])):
                    flat_pts_k.append(pts[i, 0, 0])
                    flat_pts_e.append(pts[i, 0, 1])
                    size = sizes[j]
                    flat_sizes.append(point_scale*(size)**2)
                    flat_colours.append(plot_colours[j])
                    flat_zorders.append(zorders[j])
            ax.plot(pts[:, 0, 0], pts[:, 0, 1], lw=0.5, alpha=0.5, c='grey', zorder=0)

    flat_zorders = np.asarray(flat_zorders)
    flat_pts_k = np.asarray(flat_pts_k)[np.argsort(flat_zorders)]
    flat_pts_e = np.asarray(flat_pts_e)[np.argsort(flat_zorders)]
    flat_sizes = np.asarray(flat_sizes)[np.argsort(flat_zorders)]
    flat_colours = np.asarray(flat_colours)[np.argsort(flat_zorders)]

    ax.scatter(flat_pts_k, flat_pts_e, edgecolor=flat_colours, s=flat_sizes, lw=0, marker='o', facecolor=flat_colours)


def _get_lineprops(dispersion, spin_fermi_energy, nb, ns, branch, branch_ind, seed_ind, kwargs):
    """ Get the properties of the line to plot. """
    colour = None
    alpha = 1
    label = None
    if isinstance(dispersion, ElectronicDispersion):
        if dispersion.num_spins == 2:
            if ns == 0:
                colour = 'red'
                alpha = 0.3
            else:
                colour = 'blue'
                alpha = 0.3
        else:
            if kwargs.get('band_colour') == 'occ':
                band_min = np.min(dispersion.eigs[ns][nb][branch]) - spin_fermi_energy[ns]
                band_max = np.max(dispersion.eigs[ns][nb][branch]) - spin_fermi_energy[ns]

                if band_max < 0:
                    colour = kwargs.get('valence')
                elif band_min > 0:
                    colour = kwargs.get('conduction')
                elif band_min < 0 < band_max:
                    colour = kwargs.get('crossing')

    if kwargs['colour_by_seed']:
        colour = kwargs.get('seed_colours')[seed_ind]

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

    return colour, alpha, label


def _add_path_labels(seed, dispersion, ax_dispersion, path, seed_ind, kwargs):
    """ Scrape k-point path labels from cell file and seekpath, then add them to the plot. """
    from matador.utils.cell_utils import doc2spg
    from seekpath import get_path
    xticks = []
    xticklabels = []
    shear_planes = []
    labelled = []
    path_labels = dict()

    # first, try to grab them from the cell file
    if os.path.isfile(seed + '.cell'):
        doc, success = cell2dict(seed + '.cell',
                                 db=False, verbosity=kwargs.get('verbosity', 0),
                                 lattice=True, positions=True)
        if kwargs['phonons']:
            key = 'phonon_fine_kpoint_path'
        else:
            key = 'spectral_kpoints_path'
        if key in doc and key + '_labels' in doc:
            for label, point in zip(doc[key + '_labels'], doc[key]):
                path_labels[label] = point
            print('Detected path labels from cell file')

    if not path_labels:
        # try to get dispersion path labels from spglib/seekpath
        spg_structure = None
        if isinstance(dispersion, VibrationalDispersion):
            spg_structure = doc2spg(dispersion._data)
        else:
            res = False
            cell = False
            if os.path.isfile(seed + '.res'):
                res = True
            elif os.path.isfile(seed + '.cell'):
                cell = True
            else:
                print('Failed to find {}.cell or {}.res, will not be able to generate labels.'.format(seed, seed))

            success = False
            if cell:
                doc, success = cell2dict(seed + '.cell',
                                         db=False, verbosity=kwargs.get('verbosity', 0),
                                         lattice=True, positions=True)
            if res and not success:
                doc, success = res2dict(seed + '.res',
                                        db=False, verbosity=kwargs.get('verbosity', 0))

            if cell or res:
                if success:
                    spg_structure = doc2spg(doc)
                else:
                    print('Failed to scrape {}.cell/.res, will not be able to generate labels.'.format(seed))

        if spg_structure:
            seekpath_results = get_path(spg_structure)
            path_labels = seekpath_results['point_coords']

    for branch_ind, branch in enumerate(dispersion.kpoint_branches):
        for sub_ind, ind in enumerate(branch):
            kpt = dispersion.kpoint_path[ind]
            for label, point in path_labels.items():
                if np.allclose(point, kpt):
                    if ind - branch_ind not in labelled:
                        label = label.replace('GAMMA', r'\Gamma')
                        label = label.replace('SIGMA', r'\Sigma')
                        label = label.replace('DELTA', r'\Delta')
                        label = label.replace('LAMBDA', r'\Lambda')
                        if sub_ind == len(branch) - 1:
                            if branch_ind < len(dispersion.kpoint_branches) - 1:
                                _tmp = dispersion.kpoint_path
                                next_point = _tmp[dispersion.kpoint_branches[branch_ind + 1][0]]
                                for new_label, new_point in path_labels.items():
                                    new_label = new_label.replace('GAMMA', r'\Gamma')
                                    new_label = new_label.replace('SIGMA', r'\Sigma')
                                    new_label = new_label.replace('DELTA', r'\Delta')
                                    new_label = new_label.replace('LAMBDA', r'\Lambda')
                                    # import matplotlib
                                    if np.allclose(new_point, next_point):
                                        label = '\\dfrac{{{}}}{{{}}}'.format(label, new_label)
                                        ax_dispersion.axvline(path[ind - branch_ind], ls='-', c='grey', zorder=1, lw=0.5)
                                        labelled.append(ind - branch_ind)
                                        shear_planes.append(ind)
                        label = '${}$'.format(label.replace('$', ''))
                        ax_dispersion.axvline(path[ind - branch_ind], ls='--', c='grey', zorder=0, lw=0.5)
                        xticklabels.append(label)
                        xticks.append(path[ind - branch_ind])
                        break

    if isinstance(dispersion, ElectronicDispersion) and kwargs['gap']:
        if dispersion.num_spins != 1:
            raise NotImplementedError('Band gap summary not implemented for multiple spin channels.')
        if dispersion.band_gap > 0:
            vbm_pos = dispersion['band_gap_path_inds'][1]
            vbm = dispersion['valence_band_min'] - dispersion.fermi_energy
            cbm_pos = dispersion['band_gap_path_inds'][0]
            cbm = dispersion['conduction_band_max'] - dispersion.fermi_energy
            if vbm_pos != cbm_pos:
                vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                ax_dispersion.plot([path[vbm_pos - vbm_offset], path[cbm_pos - cbm_offset]], [vbm, cbm],
                                   ls=kwargs['ls'][seed_ind],
                                   c='blue',
                                   label='indirect gap {:3.3f} eV'.format(cbm - vbm))

            vbm_pos = dispersion['direct_gap_path_inds'][1]
            vbm = dispersion['direct_valence_band_min'] - dispersion.fermi_energy
            cbm_pos = dispersion['direct_gap_path_inds'][0]
            cbm = dispersion['direct_conduction_band_max'] - dispersion.fermi_energy
            vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
            cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
            ax_dispersion.plot([path[vbm_pos - vbm_offset], path[cbm_pos - cbm_offset]], [vbm, cbm],
                               ls=kwargs['ls'][seed_ind],
                               c='red',
                               label='direct gap {:3.3f} eV'.format(cbm - vbm))
            ax_dispersion.legend(loc='upper center',
                                 bbox_to_anchor=(0.5, 1.1),
                                 fancybox=True, shadow=True,
                                 ncol=2, handlelength=1)

    if seed_ind == 0:
        ax_dispersion.set_xticks(xticks)
        ax_dispersion.set_xticklabels(xticklabels)
        ax_dispersion.grid(False)


def _get_projector_info(projectors):
    """ Grab appropriate colours and labels from a list of projectors.

    Parameters:
        projectors (list): list containing (element_str, l_channel) tuples.

    Returns:
        list: list of projector labels, e.g. {element_str}-${l_channel}$.
        list: list of colours for density of states, derived from vesta colours.

    """

    import matplotlib.pyplot as plt
    element_colours = get_element_colours()
    projector_labels = []
    dos_colours = []
    for ind, projector in enumerate(projectors):
        if projector[0] is None:
            projector_label = '${}$-character'.format(projector[1])
        elif projector[1] is None:
            projector_label = projector[0]
        else:
            projector_label = '{p[0]} (${p[1]}$)'.format(p=projector)
        projector_labels.append(projector_label)

        # if species-projected only, then use VESTA colours
        if projector[0] is not None and projector[1] is None:
            dos_colours.append(element_colours.get(projector[0]))
        # if species_ang-projected, then use VESTA colours but lightened
        elif projector[0] is not None and projector[1] is not None:
            from copy import deepcopy
            dos_colour = deepcopy(element_colours.get(projector[0]))
            multi = ['s', 'p', 'd', 'f'].index(projector[1]) - 1
            for jind, _ in enumerate(dos_colour):
                dos_colour[jind] = max(min(dos_colour[jind]+multi*0.2, 1), 0)
            dos_colours.append(dos_colour)
        # otherwise if just ang-projected, use colour_cycle
        else:
            dos_colours.append(list(plt.rcParams['axes.prop_cycle'].by_key()['color'])[ind])

    return projector_labels, dos_colours
