# coding: utf-8
# Distributed under the terms of the MIT License.

""" Run calculations in a folder with BatchRun
such that there are no clashes.
"""
import os
import argparse
from matador import __version__, script_epilog
from matador.utils.print_utils import print_notify
from matador.utils.errors import InputError
from matador.compute import BatchRun


def main():
    """ Parse args and run any remaining jobs. """
    parser = argparse.ArgumentParser(
        prog='run3',
        description='Run multiple calculations from a series of .res \
                     files and single cell and param files, typically CASTEP geometry \
                     optimisations. The geometry optimization will \
                     be split into "--rough" (default: 4) chunks of "--rough_iter" (default: 2) \
                     iterations, followed by chunks of "--fine_iter" (default: 20) \
                     iterations, until geom_max_iter is reached in the param file. \
                     Successful runs will be moved to "completed", crashes/failures will go to \
                     "bad_castep" and initial res files will go into "input". Running jobs will \
                     be listed in jobs.txt and those that completed cleanly will be listed \
                     in finished_cleanly.txt.',
        epilog=script_epilog
    )
    parser.add_argument('--version', action='version', version='matador version ' + __version__ + '.')
    parser.add_argument('seed', type=str, nargs='+',
                        help='cell and param seed to use as template for CASTEP calculations OR list of files\
                              to apply run executable on')
    parser.add_argument('-nc', '--ncores', type=int,
                        help='number of cores per node per job [DEFAULT=cpu_count/nprocesses]')
    parser.add_argument('-np', '--nprocesses', type=int, default=1,
                        help='number of concurrent calculations, i.e. number \
                              of concurrent mpiruns [DEFAULT=1]')
    parser.add_argument('-nn', '--nnodes', type=int, default=1,
                        help='number of nodes per job, i.e. number of nodes \
                              using -nc cores [DEFAULT=1].')

    parser.add_argument('-t', '--max_walltime', type=int,
                        help='maximum walltime in seconds (job will quit early to clean up if specified)')
    parser.add_argument('-exec', '--executable', type=str,
                        help='specify path to or name of executable (DEFAULT: castep)')
    parser.add_argument('--no_reopt', action='store_true', default=False,
                        help='do not run geometry optimisation again after first success')
    parser.add_argument('--redirect', type=str,
                        help='filename to redirect output to, can use $seed macro')
    parser.add_argument('--mode', type=str, default='castep',
                        help='either castep or generic')
    parser.add_argument('--noise', action='store_true',
                        help=('add 0.1 A of random noise to positions on every cell, '
                              'useful for converging forces (DEFAULT: off)'))
    parser.add_argument('--squeeze', action='store_true',
                        help='add external pressure to the rough steps of geom opts')
    parser.add_argument('--ignore_jobs_file', action='store_true',
                        help='whether to use the jobs.txt file to avoid clashes')
    parser.add_argument('-d', '--debug', action='store_true', default=False,
                        help='debug output')
    parser.add_argument('-cust', '--custom_params', action='store_true', default=False,
                        help='use custom param file per structure')
    parser.add_argument('-v', '--verbosity', type=int, default=2,
                        help='integer to set level of verbosity')
    parser.add_argument('--archer', action='store_true', default=False,
                        help='use aprun over mpirun')
    parser.add_argument('--slurm', action='store_true', default=False,
                        help='use srun over mpirun')
    parser.add_argument('--intel', action='store_true', default=False,
                        help='use Intel\'s mpirun')
    parser.add_argument('--conv_cutoff', action='store_true', default=False,
                        help='run all res files at cutoff defined in cutoff.conv file')
    parser.add_argument('--conv_kpt', action='store_true', default=False,
                        help='run all res files at kpoint spacings defined in kpt.conv file')
    parser.add_argument('--memcheck', action='store_true', default=False,
                        help='enable memcheck via castep dryrun')
    parser.add_argument('--scratch_prefix', type=str,
                        help='specify absolute path prefix for compute dir e.g. '
                             '--scratch_prefix /scratch/user/ will set the compute directory to /scratch/user/$hostname. '
                             'default value is taken from .matadorrc.')
    parser.add_argument('--maxmem', type=int,
                        help='override max memory for memcheck')
    parser.add_argument('--killcheck', action='store_true', default=True,
                        help='check for $seed.kill file and quit job if present')
    parser.add_argument('--kpts_1D', action='store_true', default=False,
                        help='recalculate a 1D kpoint mesh of spacing specified in template cell')
    parser.add_argument('--spin', type=int, nargs='?', const=5, default=None,
                        help=('if not specified in .cell file, break spin symmetry on first atom using the spin specified by '
                              'the user [DEFAULT: 5]'))
    parser.add_argument('--rough', type=int, default=4,
                        help='choose how many <rough_iter> geometry optimizations \
                              to perform, decrease if lattice is nearly correct. [DEFAULT: 4].')
    parser.add_argument('--rough_iter', type=int, default=2,
                        help='choose how many relaxation steps per rough geometry optimization\
                              to perform, [DEFAULT: 2].')
    parser.add_argument('--fine_iter', type=int, default=20,
                        help='choose how many relaxation steps per fine geometry optimization\
                              to perform, [DEFAULT: 20].')
    parser.add_argument('-l', '--limit', type=int, default=None,
                        help='limit to n structures per run')
    parser.add_argument('--profile', action='store_true',
                        help='profile code with cProfile')
    args = parser.parse_args()

    seed = vars(args)['seed']

    kwargs = vars(args)
    del kwargs['seed']

    from matador.config import load_custom_settings
    settings = load_custom_settings(debug=kwargs.get('debug')).get('run3')
    if settings is not None:
        if kwargs['scratch_prefix'] is None:
            kwargs['scratch_prefix'] = settings.get('scratch_prefix')
            if kwargs['scratch_prefix'] == '.':
                kwargs['scratch_prefix'] = None

        if kwargs['executable'] is None:
            kwargs['executable'] = settings.get('castep_executable', 'castep')

        kwargs['run3_settings'] = settings

    if sum([vars(args)['slurm'], vars(args)['archer'], vars(args)['intel']]) > 1:
        exit('Incompatible MPI arguments specified, please use at most one of --archer/--intel/--slurm.')

    if vars(args).get('profile'):
        import cProfile
        import pstats
        from sys import version_info
        hostname = os.uname()[1]
        pr = cProfile.Profile()
        pr.enable()

    try:
        runner = BatchRun(seed, **kwargs)
        runner.spawn()
    except InputError as exc:
        print_notify(exc)
    except RuntimeError:
        print_notify('Some jobs failed, exiting...')
    except Exception as exc:
        raise exc

    if vars(args).get('profile'):
        pr.disable()
        fname = 'run3-{}-{}-{}.{}.{}'.format(__version__, hostname, version_info.major,
                                             version_info.minor, version_info.micro)
        pr.dump_stats(fname + '.prof')
        with open(fname + '.pstats', 'w') as s:
            sortby = 'cumulative'
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            ps.print_stats()


if __name__ == '__main__':
    main()
