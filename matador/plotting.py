# coding: utf-8
""" This file implements several routines for
dispersion plots, phase diagrams, as well as
voltage and volume expansion plots.
"""

from traceback import print_exc
import numpy as np
import matplotlib.pyplot as plt
from matador.utils.chem_utils import get_formula_from_stoich
from matador.utils.print_utils import print_warning


def plot_spectral(seeds, **kwargs):
    """ Plot bandstructure and optional DOS from <seed>.bands and
    <seed>.adaptive.dat file.

    Input:

        | seeds: list, list of filenames of bands files

    Args:

        | plot_bandstructure : bool, whether to plot bandstructure
        | plot_dos           : bool, whether to plot density of states
        | dos                : str, separate seed name for pDOS/DOS data
        | phonons            : bool, whether to plot phonon or electronic data
        | labels             : list(str), legend labels for multiple bandstructures
        | cell               : bool, whether to work out correct labels from structure in cell file
        | gap                : bool, draw on the band gap
        | colour_by_seed     : bool, plot with a separate colour per bandstructure
        | external_efermi    : float, replace scraped Fermi energy with this value (eV)
        | highlight_bands    : list(int), colour the bands with these indices in red

    """
    from matador.scrapers.castep_scrapers import bands2dict, cell2dict, phonon2dict, optados2dict
    from matador.utils.cell_utils import doc2spg
    from seekpath import get_path
    import seaborn as sns
    from os.path import isfile
    # set defaults and update class with desired values
    prop_defaults = {'plot_bandstructure': True, 'plot_dos': False,
                     'phonons': False, 'cell': False, 'gap': False,
                     'colour_by_seed': False, 'external_efermi': None,
                     'labels': None,
                     'verbosity': 0, 'highlight_bands': None}
    prop_defaults.update(kwargs)
    kwargs = prop_defaults

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

    if len(seeds) > 1 or kwargs.get('colour_by_seed'):
        seed_colours = colours
        ls = ['-']*len(seeds)
        colour_by_seed = True
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

    if kwargs.get('plot_window') is not None:
        if isinstance(kwargs.get('plot_window'), list):
            if len(kwargs.get('plot_window')) != 2:
                exit('plot_window must have length 2 or be a single number')
            plot_window = sorted(kwargs.get('plot_window'))
        else:
            plot_window = (-kwargs.get('plot_window'), kwargs.get('plot_window'))
    else:
        plot_window = (-5, 5)

    if kwargs['plot_bandstructure'] and not kwargs['plot_dos']:
        fig, ax_dispersion = plt.subplots(figsize=(5, 5))
    elif kwargs['plot_bandstructure'] and kwargs['plot_dos']:
        fig, ax_grid = plt.subplots(1, 2, figsize=(8.5, 5), sharey=True,
                                    gridspec_kw={'width_ratios': [4, 1],
                                                 'wspace': 0.05,
                                                 'left': 0.15})
        ax_dispersion = ax_grid[0]
        ax_dos = ax_grid[1]
    elif not kwargs['plot_bandstructure'] and kwargs['plot_dos']:
        fig, ax_dos = plt.subplots(figsize=(7, 5))

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
                plot_window = [np.min(dispersion[eig_key]), np.max(dispersion[eig_key])]
            else:
                dispersion, s = bands2dict(seed + '.bands',
                                           summary=True,
                                           gap=kwargs.get('gap'),
                                           external_efermi=kwargs.get('external_efermi'),
                                           verbosity=kwargs.get('verbosity'))
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
            if kwargs['phonons']:
                dispersion[eig_key] = modes2bands(dispersion[eig_key], path, dispersion[branch_key])
            for branch_ind, branch in enumerate(dispersion[branch_key]):
                for ns in range(dispersion[spin_key]):
                    for nb in range(dispersion[band_key]):
                        if not kwargs['phonons']:
                            if dispersion[spin_key] == 2:
                                if ns == 0:
                                    colour = 'red'
                                else:
                                    colour = 'blue'
                            else:
                                # if np.max(dispersion[eig_key][ns][nb][branch]) < 0:
                                    # if colour_by_seed:
                                        # colour = seed_colours[seed_ind]
                                    # else:
                                        # colour = valence
                                # elif np.min(dispersion[eig_key][ns][nb][branch]) > 0:
                                    # if colour_by_seed:
                                        # colour = seed_colours[seed_ind]
                                    # else:
                                        # colour = conduction
                                # elif np.min(dispersion[eig_key][ns][nb][branch]) < 0 and np.max(dispersion[eig_key][ns][nb][branch]) > 0:
                                    # if colour_by_seed:
                                        # colour = seed_colours[seed_ind]
                                    # else:
                                        # colour = crossing
                                # else:
                                if colour_by_seed:
                                    colour = seed_colours[seed_ind]
                                else:
                                    colour = 'black'
                        else:
                            colour = 'black'
                        if kwargs.get('highlight_bands') is not None and nb in kwargs.get('highlight_bands'):
                            colour = 'red'
                            alpha = 0.5
                        else:
                            alpha = 1
                        if branch_ind == 0 and ns == 0 and nb == 0 and kwargs.get('labels') is not None:
                            label = kwargs.get('labels')[seed_ind]
                        else:
                            label = None
                        ax_dispersion.plot(path[(np.asarray(branch)-branch_ind).tolist()], dispersion[eig_key][ns][nb][branch],
                                           c=colour, lw=1, marker=None, ls=ls[seed_ind], alpha=alpha, label=label)
                if branch_ind != len(dispersion[branch_key])-1:
                    ax_dispersion.axvline(path[branch[-1]-branch_ind], ls='-.', lw=1, c='grey')
            if kwargs.get('labels') is not None:
                ax_dispersion.legend()
            ax_dispersion.axhline(0, ls='--', lw=1, c='grey')
            ax_dispersion.set_ylim(plot_window)
            if kwargs['phonons']:
                ylabel = 'Wavenumber (cm$^{-1}$)'
            else:
                ylabel = '$\epsilon_k$ (eV)'
            ax_dispersion.set_ylabel(ylabel)
            ax_dispersion.set_xlim(-0.05, 1.05)
            if kwargs['phonons']:
                dispersion['freq_unit'] = dispersion['freq_unit'].replace('-1', '$^{-1}$')
                ax_dispersion.axhline(np.min(dispersion['softest_mode_freq']), ls=ls[seed_ind], c='r',
                                      label='$\omega_\mathrm{{min}} = {:5.3f}$ {}'.format(np.min(dispersion['softest_mode_freq']), dispersion['freq_unit']))
                ax_dispersion.legend()
            xticks = []
            xticklabels = []
            shear_planes = []
            labelled = []

            if kwargs['cell']:
                doc, success = cell2dict(seed + '.cell', db=False, verbosity=kwargs['verbosity'], outcell=True, positions=True)
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
            if not kwargs['phonons'] and kwargs['gap'] and dispersion['band_gap'] > 0:
                vbm_pos = dispersion['band_gap_path_inds'][1]
                vbm = dispersion['valence_band_min']
                cbm_pos = dispersion['band_gap_path_inds'][0]
                cbm = dispersion['conduction_band_max']
                if vbm_pos != cbm_pos:
                    vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                    cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                    ax_dispersion.plot([path[vbm_pos-vbm_offset], path[cbm_pos-cbm_offset]], [vbm, cbm], ls=ls[seed_ind], c='blue', label='indirect gap {:3.3f} eV'.format(cbm-vbm))
                vbm_pos = dispersion['direct_gap_path_inds'][1]
                vbm = dispersion['direct_valence_band_min']
                cbm_pos = dispersion['direct_gap_path_inds'][0]
                cbm = dispersion['direct_conduction_band_max']
                vbm_offset = sum([vbm_pos > ind for ind in shear_planes])
                cbm_offset = sum([cbm_pos > ind for ind in shear_planes])
                ax_dispersion.plot([path[vbm_pos-vbm_offset], path[cbm_pos-cbm_offset]], [vbm, cbm], ls=ls[seed_ind], c='red', label='direct gap {:3.3f} eV'.format(cbm-vbm))
                ax_dispersion.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True, ncol=2, handlelength=1)
                # else:
                    # ax_dispersion.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True, ncol=1, handlelength=1)
            ax_dispersion.set_xticks(xticks)
            ax_dispersion.set_xticklabels(xticklabels)
        if kwargs['plot_dos']:
            if not kwargs['phonons']:
                if kwargs.get('dos') is None:
                    dos_data, s = optados2dict(seed + '.adaptive.dat')
                else:
                    dos_data, s = optados2dict(kwargs.get('dos'))
                energies = dos_data['energies']
                dos = dos_data['dos']
                max_density = np.max(dos[np.where(energies > plot_window[0])])
                if 'pdos' in dos_data:
                    pdos = dos_data['pdos']
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
                            projector_labels = line.split()[5:]
                            num_projectors = len(projector_labels)
                            begin = ind+1
                            break
                    print('Found {} projectors'.format(num_projectors))
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
                    dos_data['pdos'] = dict()
                    for i, label in enumerate(projector_labels):
                        dos_data['pdos'][label] = raw_data[:, i+2]
                    pdos = dos_data['pdos']
                    from os import remove
                    remove(seed + '.phonon_dos_tmp')
                    plot_window = [-10, 350]

            if kwargs['phonons']:
                ylabel = 'Phonon DOS'
                xlabel = 'Wavenumber (cm$^{{-1}}$)'
                # xlabel = 'Wavenumber ({})'.format(phonon_data['freq_unit'])
            else:
                ylabel = 'DOS'
                xlabel = 'Energy (eV)'
            if kwargs['plot_bandstructure']:
                ax_dos.set_xticks([0.6*max_density])
                ax_dos.set_xticklabels([ylabel])
                ax_dos.axhline(0, c='grey', ls='--', lw=1)
                ax_dos.set_xlim(0, max_density*1.2)
                ax_dos.set_ylim(plot_window)
                ax_dos.axvline(0, c='k')
                ax_dos.plot(dos, energies, lw=1, c='k', ls=ls[seed_ind])
                if 'pdos' in dos_data:
                    for ind, projector in enumerate(pdos):
                        ax_dos.plot(pdos[projector], energies, lw=1, zorder=1000)
                        ax_dos.fill_betweenx(energies, 0, pdos[projector], alpha=0.3, label=projector)
                ax_dos.legend()

                if seed_ind == 0 and not kwargs['phonons']:
                    if 'pdos' not in dos_data:
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
                if 'pdos' in dos_data:
                    for ind, projector in enumerate(pdos):
                        ax_dos.plot(energies, pdos[projector], lw=1, zorder=1000)
                        ax_dos.fill_between(energies, 0, pdos[projector], alpha=0.3, label=projector)
                ax_dos.legend()
                if seed_ind == 0 and not kwargs['phonons']:
                    if 'pdos' not in dos_data:
                        ax_dos.fill_between(energies[np.where(energies <= 0)], 0, dos[np.where(energies <= 0)], facecolor=valence, alpha=0.5)
                        ax_dos.fill_between(energies[np.where(energies >= 0)], 0, dos[np.where(energies >= 0)], facecolor=conduction, alpha=0.5)
            ax_dos.grid('off')

    if kwargs.get('pdf'):
        plt.savefig(seed.replace('.bands', '').replace('.phonon', '') + '_spectral.pdf', bbox_inches='tight')
    if kwargs.get('svg'):
        plt.savefig(seed.replace('.bands', '').replace('.phonon', '') + '_spectral.svg', bbox_inches='tight', transparent=True)
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


def plot_voltage_curve(hull, show=False):
    """ Plot calculated voltage curve. """
    import matplotlib.pyplot as plt
    import numpy as np
    from traceback import print_exc
    if hull.savefig:
        if len(hull.voltages) != 1:
            fig = plt.figure(facecolor=None, figsize=(4, 3.5))
        else:
            fig = plt.figure(facecolor=None, figsize=(4, 3.5))
    else:
        fig = plt.figure(facecolor=None)
    axQ = fig.add_subplot(111)
    if hull.args.get('expt') is not None:
        try:
            expt_data = np.loadtxt(hull.args.get('expt'), delimiter=',')
        except:
            print_exc()
            pass
        if hull.args.get('expt_label'):
            axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label=hull.args.get('expt_label'))
        else:
            axQ.plot(expt_data[:, 0], expt_data[:, 1], c='k', lw=2, ls='-', label='Experiment')
    for ind, voltage in enumerate(hull.voltages):
        for i in range(len(voltage)-1):
            if i == 0 and hull.args.get('expt'):
                axQ.plot([hull.Q[ind][i-1], hull.Q[ind][i]], [voltage[i], voltage[i]], marker='*',
                         lw=2, c=hull.colours[ind], label='DFT (this work)')
            elif i == 0 and len(hull.voltages) != 1:
                axQ.plot([hull.Q[ind][i-1], hull.Q[ind][i]], [voltage[i], voltage[i]], marker='o', markersize=5,
                         lw=2, c=hull.colours[ind], label=get_formula_from_stoich(hull.endstoichs[ind], tex=True))
            else:
                axQ.plot([hull.Q[ind][i-1], hull.Q[ind][i]], [voltage[i], voltage[i]], marker='o', markersize=5,
                         lw=2, c=hull.colours[ind])
            if i != len(voltage)-2:
                axQ.plot([hull.Q[ind][i], hull.Q[ind][i]], [voltage[i], voltage[i+1]], marker='o', markersize=5,
                         lw=2, c=hull.colours[ind])
    if hull.args.get('labels'):
        eps = 1e-9
        hull.label_cursor = [doc for doc in hull.hull_cursor if doc['hull_distance'] <= 0+eps]
        hull.label_cursor = hull.label_cursor[1:-1]
        for i in range(len(hull.label_cursor)):
            axQ.annotate(get_formula_from_stoich(hull.label_cursor[i]['stoichiometry'],
                                                 elements=hull.elements, tex=True),
                         xy=(hull.Q[0][i+1]+0.02*max(hull.Q[0]), hull.voltages[0][i+1]+0.02*max(hull.voltages[0])),
                         textcoords='data',
                         ha='center',
                         zorder=9999)
    if hull.args.get('expt') or len(hull.voltages) != 1:
        axQ.legend(loc=1)
    axQ.set_ylabel('Voltage (V) vs {}$^+$/{}'.format(hull.elements[0], hull.elements[0]))
    axQ.set_xlabel('Gravimetric cap. (mAh/g)')
    start, end = axQ.get_ylim()
    axQ.set_ylim(0, 1.1*end)
    start, end = axQ.get_xlim()
    axQ.set_xlim(0, 1.1*end)
    axQ.grid('off')
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
    try:
        import seaborn as sns
        sns.despine()
        dark_grey = '#262626'
        for spine in ['left', 'bottom']:
            axQ.spines[spine].set_linewidth(0.5)
            axQ.spines[spine].set_color(dark_grey)
    except:
        pass

    if hull.savefig:
        if hull.args.get('pdf'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_voltage.pdf',
                        dpi=500, transparent=True)
        if hull.args.get('svg'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_voltage.svg',
                        dpi=500, transparent=True)
        if hull.args.get('png'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_voltage.png',
                        dpi=500, transparent=True)
    elif show:
        plt.show()


def plot_2d_hull_bokeh(hull):
    """ Plot interactive hull with Bokeh. """
    from bokeh.plotting import figure, save, output_file
    from bokeh.models import ColumnDataSource, HoverTool, Range1d
    import matplotlib.pyplot as plt
    import numpy as np
    from .utils.glmol_wrapper import get_glmol_placeholder_string
    from .export import generate_hash, generate_relevant_path
    # grab tie-line structures
    tie_line_data = dict()
    tie_line_data['composition'] = list()
    tie_line_data['energy'] = list()
    for ind in range(len(hull.hull.vertices)):
        if hull.structure_slice[hull.hull.vertices[ind], 1] <= 0:
            tie_line_data['composition'].append(hull.structure_slice[hull.hull.vertices[ind], 0])
            tie_line_data['energy'].append(hull.structure_slice[hull.hull.vertices[ind], 1])
    tie_line_data['energy'] = np.asarray(tie_line_data['energy'])
    tie_line_data['composition'] = np.asarray(tie_line_data['composition'])
    tie_line_data['energy'] = tie_line_data['energy'][np.argsort(tie_line_data['composition'])]
    tie_line_data['composition'] = np.sort(tie_line_data['composition'])

    # points for off hull structures
    hull_data = dict()
    hull_data['composition'] = hull.structures[:, 0]
    hull_data['energy'] = hull.structures[:, 1]
    hull_data['hull_distance'] = hull.hull_dist
    hull_data['formula'], hull_data['text_id'] = [], []
    hull_data['space_group'], hull_data['hull_dist_string'] = [], []
    for structure in hull.info:
        hull_data['formula'].append(structure[0])
        hull_data['text_id'].append(structure[1])
        hull_data['space_group'].append(structure[2])
        hull_data['hull_dist_string'].append(structure[3])
    cmap_limits = [0, 0.5]
    colormap = plt.cm.get_cmap('Dark2')
    cmap_input = np.interp(hull_data['hull_distance'], cmap_limits, [0.15, 0.4], left=0.15, right=0.4)
    colours = colormap(cmap_input, 1, True)
    bokeh_colours = ["#%02x%02x%02x" % (r, g, b) for r, g, b in colours[:, 0:3]]
    fixed_colours = colormap([0.0, 0.15], 1, True)
    tie_line_colour, on_hull_colour = ["#%02x%02x%02x" % (r, g, b) for r, g, b in fixed_colours[:, 0:3]]

    tie_line_source = ColumnDataSource(data=tie_line_data)
    hull_source = ColumnDataSource(data=hull_data)

    hover = HoverTool(tooltips="""
                      <div>
                          <div>
                              <span style="font-size: 16px; font-family: "Fira Sans", sans-serif">
                                  Formula: @formula <br>
                                  ID: @text_id <br>
                                  Space group: @space_group <br>
                                  Distance from hull: @hull_dist_string
                              </span>
                          </div>
                      </div>
                      """)

    tools = ['pan', 'wheel_zoom', 'reset', 'save']
    tools.append(hover)

    fig = figure(tools=tools)

    fig.xaxis.axis_label = 'x'
    fig.yaxis.axis_label = 'Formation energy (eV/atom)'
    fig.xaxis.axis_label_text_font_size = '20pt'
    fig.xaxis.axis_label_text_font = "Fira Sans, sans-serif"
    fig.yaxis.axis_label_text_font_size = '20pt'
    fig.yaxis.axis_label_text_font = "Fira Sans, sans-serif"
    fig.yaxis.axis_label_text_font_style = 'normal'
    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.background_fill_alpha = 0
    fig.border_fill_alpha = 0
    fig.title.text_font_size = '20pt'
    fig.title.align = 'center'

    ylim = [-0.1 if np.min(hull.structure_slice[hull.hull.vertices, 1]) > 0
            else np.min(hull.structure_slice[hull.hull.vertices, 1])-0.15,
            0.1 if np.max(hull.structure_slice[hull.hull.vertices, 1]) > 1
            else np.max(hull.structure_slice[hull.hull.vertices, 1])+0.1]
    fig.x_range = Range1d(-0.1, 1.1)
    fig.y_range = Range1d(ylim[0], ylim[1])

    fig.line('composition', 'energy',
             source=tie_line_source,
             line_width=4,
             line_color=tie_line_colour)
    hull_scatter = fig.scatter('composition', 'energy',
                               source=hull_source,
                               alpha=1,
                               size=10,
                               fill_color=bokeh_colours,
                               line_color=None)
    fig.tools[0].renderers.append(hull_scatter)
    fig.square('composition', 'energy',
               source=tie_line_source,
               line_color='black',
               color=on_hull_colour,
               line_width=2,
               alpha=1,
               size=10)
    fig.plot_width = 800
    fig.plot_height = 600
    path = '/u/fs1/me388/data/hulls/'
    fname = generate_relevant_path(hull.args) + '_' + generate_hash() + '.html'
    output_file(path+fname, title='Convex hull')
    print('Hull will be available shortly at http://www.tcm.phy.cam.ac.uk/~me388/hulls/' + fname)
    save(fig)
    glmol = False
    if glmol:
        html_string, js_string = get_glmol_placeholder_string()
        with open(path+fname) as f:
            flines = f.readlines()
            for ind, line in enumerate(flines):
                if "<div class=\"bk-root\">" in line:
                    flines.insert(ind - 1, html_string)
                    break
            flines.append(js_string)
        with open(path+fname, 'w') as f:
            f.write('\n'.join(map(str, flines)))


def plot_3d_ternary_hull(hull):
    """ Plot calculated ternary hull in 3D. """
    from mpl_toolkits.mplot3d import axes3d
    import numpy as np
    from matador.utils.hull_utils import barycentric2cart
    # avoids annoying flake8 warning
    del axes3d
    print('WARNING: deprecated')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    coords = barycentric2cart(hull.structures)
    stable = coords[np.where(hull.hull_dist < 0 + 1e-9)]
    stable = np.asarray(stable)
    ax.plot_trisurf(stable[:, 0], stable[:, 1], stable[:, 2], cmap=plt.cm.gnuplot, linewidth=1, color='grey', alpha=0.2)
    ax.scatter(stable[:, 0], stable[:, 1], stable[:, 2], s=100, c='k', marker='o')
    if len(hull.failed_structures) > 0:
        ax.scatter(coords[hull.failed_structures, 0], coords[hull.failed_structures, 1], coords[hull.failed_structures, 2], c='r')
    ax.set_zlim(-1, 1)
    ax.view_init(-90, 90)


def get_text_info(cursor, hull=False, html=False):
    """ Grab textual info for Bokeh plot labels. """
    info = []
    if hull:
        stoich_strings = []
    for ind, doc in enumerate(cursor):
        stoich_string = ''
        for elem in doc['stoichiometry']:
            stoich_string += elem[0]
            stoich_string += '$_{' + str(elem[1]) + '}$' if elem[1] != 1 else ''
        if hull:
            if stoich_string not in stoich_strings:
                stoich_strings.append(stoich_string)
        info_string = "{0:^10}\n{1:^24}\n{2:^5s}\n{3:.3f} eV".format(stoich_string,
                                                                     doc['text_id'][0] + ' ' + doc['text_id'][1],
                                                                     doc['space_group'],
                                                                     doc['hull_distance'])
        if html:
            for char in ['$', '_', '{', '}']:
                info_string = info_string.replace(char, '')
            info_string = info_string.split('\n')
        info.append(info_string)
    if hull:
        info = stoich_strings
        return info


def get_linear_cmap(colours, N=100, list_only=False):
    """ Create a linear colormap from a list of colours.

    Input:

        | colours: list(str), list of fractional RGB values of colours

    Args:

        | N        : int, number of colours in resulting cmap
        | list_only: bool, return only a list of colours

    Returns:

        | cmap: LinearSegmentedColormap or list if list_only.

    """
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
    uniq_colours = []
    _colours = [tuple(colour) for colour in colours]
    for colour in _colours:
        if colour not in uniq_colours:
            uniq_colours.append(colour)
    _colours = uniq_colours
    linear_cmap = []
    repeat = int(N/len(_colours))
    for ind, colour in enumerate(_colours):
        if ind == len(_colours) - 1:
            break
        diff = np.asarray(_colours[ind+1]) - np.asarray(_colours[ind])
        diff_norm = diff / repeat
        for i in range(repeat):
            linear_cmap.append(np.asarray(colour) + i*diff_norm)
    if list_only:
        return linear_cmap
    else:
        return LinearSegmentedColormap.from_list('linear_cmap', linear_cmap, N=N)


def plot_volume_curve(hull, show=False):
    """ Plot calculated volume curve. """
    from matador.utils.cursor_utils import get_array_from_cursor
    from matador.utils.chem_utils import get_generic_grav_capacity
    if hull.savefig:
        fig = plt.figure(facecolor=None, figsize=(4, 3.5))
    else:
        fig = plt.figure(facecolor=None)
    ax = fig.add_subplot(111)
    stable_hull_dist = get_array_from_cursor(hull.hull_cursor, 'hull_distance')
    hull_vols = []
    hull_comps = []
    bulk_vol = hull.vol_per_y[-1]
    for i in range(len(hull.vol_per_y)):
        if stable_hull_dist[i] <= 0 + 1e-16:
            hull_vols.append(hull.vol_per_y[i])
            hull_comps.append(hull.x[i])
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
        ax.scatter(hull.x[i]/(1+hull.x[i]), hull.vol_per_y[i]/bulk_vol, marker='o', s=s, edgecolor='k', lw=markeredgewidth,
                   c=c, zorder=zorder, alpha=alpha)
    hull_comps, hull_vols = np.asarray(hull_comps), np.asarray(hull_vols)
    ax.plot(hull_comps/(1+hull_comps), hull_vols/bulk_vol, marker='o', lw=4,
            c=hull.colours[0], zorder=100)
    ax.set_xlabel('$\mathrm{x}$ in $\mathrm{'+hull.elements[0]+'_x'+hull.elements[1]+'}_{1-x}$')
    ax.set_ylabel('Volume ratio with bulk')
    ax.set_ylim(0, 5*np.sort(hull_vols)[-2]/bulk_vol)
    ax.set_xlim(-0.05, 1.05)
    ax.yaxis.set_label_position('left')
    ax.set_xticklabels(ax.get_xticks())
    ax.grid('off')
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    tick_locs = np.linspace(0, 1, 6, endpoint=True).tolist()
    ax2.set_xticks(tick_locs)
    new_tick_labels = ['{}'.format(int(get_generic_grav_capacity([loc, 1-loc], [hull.elements[0], hull.elements[1]]))) for loc in tick_locs[:-1]]
    new_tick_labels.append('$\infty$')
    ax2.set_xlabel('Gravimetric capacity (mAh/g)')
    ax2.set_xticklabels(new_tick_labels)
    ax2.grid('off')
    try:
        import seaborn as sns
        sns.despine(top=False, right=False)
        dark_grey = '#262626'
        for spine in ['left', 'top', 'right', 'bottom']:
            ax.spines[spine].set_color(dark_grey)
            ax2.spines[spine].set_color(dark_grey)
            ax.spines[spine].set_linewidth(0.5)
            ax2.spines[spine].set_linewidth(0.5)
    except:
        pass
    # ax.yaxis.set_ticks(range(0, int(end)+1, 5))
    plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
    if hull.savefig:
        if hull.args.get('pdf'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_volume.pdf',
                        dpi=300)
        if hull.args.get('svg'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_volume.svg',
                        dpi=300)
        if hull.args.get('png'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_volume.png',
                        dpi=300, bbox_inches='tight')
    elif show:
        plt.show()


def plot_2d_hull(hull, ax=None, dis=False, show=False, plot_points=True,
                 plot_hull_points=True, labels=False, colour_by_source=False,
                 **kwargs):
    """ Plot calculated hull, returning ax and fig objects for further editing.

    Args:

        | ax               : matplotlib axis object, an existing axis on which to plot,
        | show             : bool, whether or not to display the plot in an X window,
        | plot_points      : bool, whether or not to display off-hull structures,
        | plot_hull_points : bool, whether or not to display on-hull structures,
        | labels           : bool, whether to label formulae of hull structures.
        | colour_by_source : bool, plot and label points by their sources

    Other kwargs:

        | alpha            : float, alpha value of points when colour_by_source is True
        | sources          : list, list of possible provenances to colour when colour_by_source
                             is True (others will be grey)

    """

    import matplotlib.pyplot as plt
    import matplotlib.colors as colours
    if ax is None:
        if hull.savefig:
            fig = plt.figure(facecolor=None, figsize=(8, 6))
        else:
            fig = plt.figure(facecolor=None)
        ax = fig.add_subplot(111)
    if not hull.plot_params:
        hull.set_plot_param()
    scatter = []
    x_elem = [hull.elements[0]]
    one_minus_x_elem = list(hull.elements[1:])
    tie_line = hull.structure_slice[hull.hull.vertices]
    plt.draw()
    # star structures on hull
    if len(hull.structure_slice) != 2:
        if plot_hull_points:
            ax.scatter(tie_line[:, 0], tie_line[:, 1],
                       c=hull.colours[1],
                       marker='o', zorder=99999, edgecolor='k',
                       s=hull.scale*40, lw=1.5, alpha=1)
        ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1],
                c=hull.colours[0], lw=2, alpha=1, zorder=1000)
        if hull.hull_cutoff > 0:
            ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1] + hull.hull_cutoff,
                    '--', c=hull.colours[1], lw=1, alpha=0.5, zorder=1000, label='')
        # annotate hull structures
        if hull.args.get('labels') or labels:
            eps = 1e-9
            hull.label_cursor = [doc for doc in hull.hull_cursor if doc['hull_distance'] <= 0+eps]
            if len(set([get_formula_from_stoich(doc['stoichiometry']) for doc in hull.label_cursor])) > len(hull.label_cursor):
                tmp_cursor = []
                for doc in hull.label_cursor:
                    if doc['stoichiometry'] not in [_doc['stoichiometry'] for _doc in tmp_cursor]:
                        tmp_cursor.append(doc)
                if len(tmp_cursor) != len(set([doc['stoichiometry'] for doc in hull.label_cursor])):
                    print_warning('Something has gone wrong with labels...')
                else:
                    hull.label_cursor = tmp_cursor
            # remove chemical potentials
            hull.label_cursor = hull.label_cursor[1:-1]
            for ind, doc in enumerate(hull.label_cursor):
                arrowprops = dict(arrowstyle="-|>", color='k')
                if (ind+2) < np.argmin(tie_line[:, 1]):
                    position = (0.8*tie_line[ind+2, 0], 1.15*(tie_line[ind+2, 1])-0.05)
                elif (ind+2) == np.argmin(tie_line[:, 1]):
                    position = (tie_line[ind+2, 0], 1.15*(tie_line[ind+2, 1])-0.05)
                else:
                    position = (min(1.1*tie_line[ind+2, 0]+0.15, 0.95), 1.15*(tie_line[ind+2, 1])-0.05)
                ax.annotate(get_formula_from_stoich(doc['stoichiometry'], elements=hull.elements, tex=True),
                            xy=(tie_line[ind+2, 0], tie_line[ind+2, 1]),
                            xytext=position,
                            textcoords='data',
                            ha='right',
                            va='bottom',
                            arrowprops=arrowprops,
                            zorder=1)
        lw = hull.scale * 0 if hull.mpl_new_ver else 1
        if plot_points and not colour_by_source:
            # points for off hull structures
            if hull.hull_cutoff == 0:
                # if no specified hull cutoff, ignore labels and colour
                # by distance from hull
                cmap_full = plt.cm.get_cmap('Dark2')
                cmin = 0.15
                cmax = 0.5
                cindmin, cindmax = int(cmin*len(cmap_full.colors)), int(cmax*len(cmap_full.colors))
                cmap = colours.LinearSegmentedColormap.from_list('Dark2', cmap_full.colors[cindmin:cindmax])
                if plot_points:
                    scatter = ax.scatter(hull.structures[np.argsort(hull.hull_dist), 0][::-1],
                                         hull.structures[np.argsort(hull.hull_dist), -1][::-1],
                                         s=hull.scale*40, lw=lw, alpha=1, c=np.sort(hull.hull_dist)[::-1],
                                         edgecolor='k', zorder=10000, cmap=cmap, norm=colours.LogNorm(0.02, 2))
                    cbar = plt.colorbar(scatter, aspect=30, pad=0.02, ticks=[0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                    cbar.ax.tick_params(length=0)
                    cbar.ax.set_yticklabels([0, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28])
                    cbar.ax.yaxis.set_ticks_position('right')
                    cbar.ax.set_frame_on(False)
                    cbar.outline.set_visible(False)
                    cbar.set_label('Distance from hull (eV)')
            elif hull.hull_cutoff != 0:
                # if specified hull cutoff, label and colour those below
                c = hull.colours[1]
                for ind in range(len(hull.structures)):
                    if hull.hull_dist[ind] <= hull.hull_cutoff or hull.hull_cutoff == 0:
                        if plot_points:
                            scatter.append(ax.scatter(hull.structures[ind, 0], hull.structures[ind, 1],
                                           s=hull.scale*40, lw=lw, alpha=0.9, c=c, edgecolor='k',
                                           zorder=300))
                if plot_points:
                    ax.scatter(hull.structures[1:-1, 0], hull.structures[1:-1, 1], s=hull.scale*30, lw=lw,
                               alpha=0.3, c=hull.colours[-2],
                               edgecolor='k', zorder=10)

        elif colour_by_source:
            from matador.utils.cursor_utils import get_guess_doc_provenance
            if kwargs.get('sources') is None:
                sources = ['AIRSS', 'GA', 'ICSD', 'OQMD', 'SWAPS']
            else:
                sources = kwargs.get('sources')
            colour_choices = {source: hull.colours[ind+1] for ind, source in enumerate(sources)}
            colours = []
            concs = []
            energies = []
            for doc in hull.cursor:
                source = get_guess_doc_provenance(doc['source'])
                if source not in sources:
                    # use grey for undesired sources
                    colours.append(hull.colours[-2])
                else:
                    colours.append(colour_choices[source])
                concs.append(doc['concentration'])
                energies.append(doc['formation_enthalpy_per_atom'])

            alpha = kwargs.get('alpha')
            if alpha is None:
                alpha = 0.2

            ax.scatter(concs, energies, c=colours, alpha=alpha, s=hull.scale*20)
            for source in sources:
                ax.scatter(1e10, 1e10, c=colour_choices[source], label=source, alpha=alpha)
            ax.scatter(1e10, 1e10, c=hull.colours[-2], label='Other', alpha=alpha)
            ax.legend()

        # tie lines
        ax.set_ylim(-0.1 if np.min(hull.structure_slice[hull.hull.vertices, 1]) > 0
                    else np.min(hull.structure_slice[hull.hull.vertices, 1])-0.15,
                    0.1 if np.max(hull.structure_slice[hull.hull.vertices, 1]) > 1
                    else np.max(hull.structure_slice[hull.hull.vertices, 1])+0.1)
    else:
        scatter = []
        print_exc()
        c = hull.colours[1]
        lw = hull.scale * 0 if hull.mpl_new_ver else 1
        for ind in range(len(hull.hull_cursor)):
            if type(hull.hull_cursor[ind]['concentration']) is list:
                if plot_points:
                    scatter.append(ax.scatter(hull.hull_cursor[ind]['concentration'][0], hull.hull_cursor[ind]['formation_enthalpy_per_atom'],
                                   s=hull.scale*40, lw=1.5, alpha=1, c=c, edgecolor='k',
                                   zorder=1000))
            else:
                if plot_points:
                    scatter.append(ax.scatter(hull.hull_cursor[ind]['concentration'], hull.hull_cursor[ind]['formation_enthalpy_per_atom'],
                                   s=hull.scale*40, lw=1.5, alpha=1, c=c, edgecolor='k',
                                   zorder=1000))
            ax.plot([0, 1], [0, 0], lw=2, c=hull.colours[0], zorder=900)
        for ind in range(len(hull.structures)):
            if plot_points:
                scatter.append(ax.scatter(hull.structures[ind, 0], hull.structures[ind, 1],
                               s=hull.scale*40, lw=lw, alpha=0.9, c=c, edgecolor='k',
                               zorder=300))

    if len(one_minus_x_elem) == 1:
        ax.set_title(x_elem[0] + '$_\mathrm{x}$' + one_minus_x_elem[0] + '$_\mathrm{1-x}$')
    if hull.non_binary:
        ax.set_title(hull.chempot_search[0] + '$_\mathrm{x}$(' + hull.chempot_search[1] + ')$_\mathrm{1-x}$')
    plt.locator_params(nbins=3)
    ax.set_xlabel('x in {}$_\mathrm{{x}}${}$_\mathrm{{1-x}}$'.format(x_elem[0], one_minus_x_elem[0]))
    ax.grid(False)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.33, 0.5, 0.66, 1])
    ax.set_xticklabels(ax.get_xticks())
    ax.set_yticks(np.arange(0, np.min(hull.structure_slice[hull.hull.vertices, 1])-0.15, -0.2))
    ax.set_yticklabels(ax.get_yticks())
    ax.set_ylabel('Formation energy (eV/atom)')
    try:
        import seaborn as sns
        sns.despine(ax=ax, left=False, bottom=False)
    except:
        pass

    if hull.savefig:
        if hull.args.get('pdf'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_hull.pdf',
                        dpi=500, bbox_inches='tight', transparent=True)
        if hull.args.get('svg'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_hull.svg',
                        dpi=500, bbox_inches='tight', transparent=True)
        if hull.args.get('png'):
            plt.savefig(hull.elements[0]+hull.elements[1]+'_hull.png',
                        dpi=500, bbox_inches='tight', transparent=True)
    elif show:
        plt.show()

    return ax


def plot_ternary_hull(hull, axis=None, show=False, plot_points=True, expecting_cbar=True, **kwargs):
    """ Plot calculated ternary hull as a 2D projection.

    Takes optional matplotlib subplot axis as a parameter, and returns
    python-ternary subplot axis object.

    """
    import ternary
    import matplotlib.pyplot as plt
    import matplotlib
    import matplotlib.colors as colours
    from matador.utils.chem_utils import get_generic_grav_capacity
    try:
        import seaborn as sns
        sns.set_style({
            'axes.facecolor': 'white', 'figure.facecolor': 'white',
            'xtick.major.size': 0, 'xtick.minor.size': 0,
            'ytick.major.size': 0, 'ytick.minor.size': 0,
            'axes.linewidth': 0.0})
    except:
        print_exc()
        pass

    print('Plotting ternary hull...')
    if hull.args.get('capmap') or hull.args.get('efmap'):
        scale = 100
    elif hull.args.get('sampmap'):
        scale = 20
    else:
        scale = 1
    fontsize = matplotlib.rcParams['font.size']

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
    ax.gridlines(color='black', multiple=scale*0.1, linewidth=0.5)

    ax.clear_matplotlib_ticks()
    ticks = [float(val) for val in np.linspace(0.0, 1.0, 6)]
    ax.ticks(axis='lbr', linewidth=1, multiple=scale*0.2, offset=0.02, fontsize=fontsize-2,
             ticks=ticks, tick_formats='%.1f')

    ax.set_title(''.join(hull.elements), fontsize=fontsize+2, y=1.02)
    ax.left_axis_label(hull.elements[2], fontsize=fontsize+2)
    ax.right_axis_label(hull.elements[1], fontsize=fontsize+2)
    ax.bottom_axis_label(hull.elements[0], fontsize=fontsize+2)

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

    Ncolours = len(hull.default_cmap_list)
    min_cut = 0.01
    max_cut = 0.2
    colours_hull = hull.default_cmap_list

    cmap = hull.default_cmap
    cmap_full = plt.cm.get_cmap('Pastel2')
    pastel_cmap = colours.LinearSegmentedColormap.from_list('Pastel2', cmap_full.colors)

    for plane in hull.hull.planes:
        plane.append(plane[0])
        plane = np.asarray(plane)
        ax.plot(scale*plane, c=hull.colours[0], lw=1.5, alpha=1, zorder=98)

    if hull.args.get('pathways'):
        for phase in stable:
            if phase[0] == 0 and phase[1] != 0 and phase[2] != 0:
                ax.plot([scale*phase, [scale, 0, 0]], c='r', alpha=0.2, lw=6, zorder=99)

    # add points
    if plot_points:
        colours_list = []
        colour_metric = hull_dist
        for i in range(len(colour_metric)):
            if colour_metric[i] >= max_cut:
                colours_list.append(Ncolours-1)
            elif colour_metric[i] <= min_cut:
                colours_list.append(0)
            else:
                colours_list.append(int((Ncolours-1)*(colour_metric[i] / max_cut)))
        colours_list = np.asarray(colours_list)
        ax.scatter(scale*stable, marker='o', color=colours_hull[0], edgecolors='black', zorder=9999999,
                   s=150, lw=1.5)
        ax.scatter(scale*concs, colormap=cmap, colorbar=True, cbarlabel='Distance from hull (eV/atom)',
                   c=hull_dist, vmax=max_cut, vmin=min_cut, zorder=1000, s=40, alpha=0)
        for i in range(len(concs)):
            ax.scatter(scale*concs[i].reshape(1, 3),
                       color=colours_hull[colours_list[i]],
                       marker='o',
                       zorder=10000-colours_list[i],
                       # alpha=max(0.1, 1-2*hull_dist[i]),
                       s=70*(1-float(colours_list[i])/Ncolours)+15,
                       lw=1, edgecolors='black')

    # add colourmaps
    if hull.args.get('capmap'):
        capacities = dict()
        from ternary.helpers import simplex_iterator
        for (i, j, k) in simplex_iterator(scale):
            capacities[(i, j, k)] = get_generic_grav_capacity([float(i)/scale, float(j)/scale, float(scale-i-j)/scale], hull.elements)
        ax.heatmap(capacities, style="hexagonal", cbarlabel='Gravimetric capacity (maH/g)',
                   vmin=0, vmax=3000, cmap=pastel_cmap)
    elif hull.args.get('efmap'):
        energies = dict()
        fake_structures = []
        from ternary.helpers import simplex_iterator
        for (i, j, k) in simplex_iterator(scale):
            fake_structures.append([float(i)/scale, float(j)/scale, 0.0])
        fake_structures = np.asarray(fake_structures)
        plane_energies, _, _ = hull.get_hull_distances(fake_structures)
        ind = 0
        for (i, j, k) in simplex_iterator(scale):
            energies[(i, j, k)] = -1*plane_energies[ind]
            ind += 1
        ax.heatmap(energies, style="hexagonal", cbarlabel='Formation energy (eV/atom)',
                   vmax=0, cmap='bone')
    elif hull.args.get('sampmap'):
        sampling = dict()
        from ternary.helpers import simplex_iterator
        eps = 1.0/float(scale)
        for (i, j, k) in simplex_iterator(scale):
            sampling[(i, j, k)] = np.size(np.where((concs[:, 0] <= float(i)/scale + eps) *
                                                   (concs[:, 0] >= float(i)/scale - eps) *
                                                   (concs[:, 1] <= float(j)/scale + eps) *
                                                   (concs[:, 1] >= float(j)/scale - eps) *
                                                   (concs[:, 2] <= float(k)/scale + eps) *
                                                   (concs[:, 2] >= float(k)/scale - eps)))
        ax.heatmap(sampling, style="hexagonal", cbarlabel='Number of structures',
                   cmap='afmhot')
    plt.tight_layout()

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
