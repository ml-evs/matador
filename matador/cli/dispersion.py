# coding: utf-8
# Distributed under the terms of the MIT License.

""" This script mimics the dispersion.pl script bundled
with CASTEP. For the given bands file, a bandstructure
is created. If a <seed>.adaptive.dat file exists
then a combined BS/DOS plot will be created.
"""

from os.path import isfile
import argparse
import glob


def main():
    """ Parse args and run the script. """
    parser = argparse.ArgumentParser(
        prog='dispersion',
        description='simple plotting script for bandstructures/DOS based on matador')
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--pdf', action='store_true',
                        help='save pdf rather than showing plot in X')
    parser.add_argument('--png', action='store_true',
                        help='save png rather than showing plot in X')
    parser.add_argument('--svg', action='store_true',
                        help='save svg rather than showing plot in X')
    parser.add_argument('--labels', type=str, nargs='*',
                        help='list of legend labels')
    parser.add_argument('--dos_only', action='store_true', help='only plot DOS')
    parser.add_argument('--bs_only', action='store_true', help='only plot dispersion')
    parser.add_argument('--preserve_kspace_distance', action='store_true',
                        help='when linearising kpoint path, ensure distance in reciprocal space is conserved')
    parser.add_argument('--cmap', type=str,
                        help='matplotlib colourmap name to use')
    parser.add_argument('--n_colours', type=int,
                        help='number of colours to use from colourmap (DEFAULT: 6)')
    parser.add_argument('--no_stacked_pdos', action='store_true',
                        help='plot PDOS as overlap rather than stack')
    parser.add_argument('--no_band_reorder', action='store_true',
                        help='don\'t reorder bands based on local gradients')
    parser.add_argument('--band_reorder', action='store_true',
                        help='try to reorder bands based on local gradients')
    parser.add_argument('--pdos_hide_tot', action='store_true',
                        help='plot PDOS without total DOS, i.e. if PDOS is negative in parts')
    parser.add_argument('--band_colour', type=str,
                        help='override all other colour options with a matploblib-interpretable colour string, '
                             'or choose one of "occ" or "random".')
    parser.add_argument('-g', '--gap', action='store_true',
                        help='plot position and size of band gap')
    parser.add_argument('-ph', '--phonons', action='store_true', default=False,
                        help='plot phonon calculation, rather than electronic')
    parser.add_argument('--highlight_bands', nargs='+', type=int,
                        help='specify band numbres to highlight in plot')
    parser.add_argument('-v', '--verbosity', type=int, default=0,
                        help='control verbosity of output')
    parser.add_argument('-pw', '--plot_window', type=float,
                        help='energy window to plot either side of E_F (eV)\
                             (DEFAULT: 5 eV)')
    parser.add_argument('seed', type=str, nargs='*',
                        help='seedname or bands file path')
    kwargs = vars(parser.parse_args())
    if 'gap' in kwargs:
        gap = kwargs['gap']
        del kwargs['gap']
    assert len(kwargs.get('seed')) >= 1, 'Seed name needed!'

    seed = kwargs.get('seed')
    test_seed = kwargs.get('seed')[0]
    verbosity = kwargs.get('verbosity')
    phonons = kwargs.get('phonons')
    labels = kwargs.get('labels')
    cmap = kwargs.get('cmap')
    band_colour = kwargs.get('band_colour')
    if band_colour is None:
        band_colour = 'occ'
    if kwargs['no_band_reorder']:
        band_reorder = False
    elif kwargs['band_reorder']:
        band_reorder = True
    else:
        band_reorder = None

    from matador.plotting import plot_spectral

    del kwargs['seed']
    del kwargs['verbosity']
    del kwargs['labels']
    del kwargs['cmap']
    del kwargs['band_colour']
    del kwargs['band_reorder']
    del kwargs['no_band_reorder']

    from matador.utils.print_utils import print_failure

    if not phonons:
        bs_seed = test_seed.replace('.bands', '') + '.bands'
        bandstructure = isfile(bs_seed)
        dos_seeds = glob.glob(test_seed.replace('.bands', '') + '*.dat')
        dos = any([isfile(dos_seed) for dos_seed in dos_seeds])
        cell_seed = test_seed.replace('.bands', '') + '.cell'
        cell = isfile(cell_seed)

    elif phonons:
        phonon_seed = test_seed.replace('.phonon', '') + '.phonon'
        bandstructure = isfile(phonon_seed)
        dos_seed = test_seed.replace('.phonon', '') + '.phonon_dos'
        dos = isfile(dos_seed)
        cell_seed = test_seed.replace('.phonon', '') + '.cell'
        cell = isfile(cell_seed)

    if not dos and not bandstructure:
        print_failure('Could not find files for specified seed {}.'.format(test_seed))
        exit()

    if kwargs.get('dos_only') and dos:
        bandstructure = False

    if kwargs.get('bs_only') and bandstructure:
        dos = False

    plot_spectral(seed,
                  plot_bandstructure=bandstructure,
                  plot_dos=dos,
                  cell=cell,
                  gap=gap,
                  verbosity=verbosity,
                  labels=labels,
                  cmap=cmap,
                  band_reorder=band_reorder,
                  band_colour=band_colour,
                  **kwargs)


if __name__ == '__main__':
    main()
