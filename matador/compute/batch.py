# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the BatchRun class for chaining
ComputeTask instances across several structures with
high-throughput.

"""


from collections import defaultdict
import multiprocessing as mp
import os
import glob
import time
import random
import psutil
from matador.utils.print_utils import print_failure, print_warning
from matador.compute.queueing import get_queue_manager
from matador.scrapers.castep_scrapers import cell2dict, param2dict
from matador.compute.compute import ComputeTask
from matador.utils.errors import (
    InputError, CalculationError,
    MaxMemoryEstimateExceeded, NodeCollisionError
)


class BatchRun:
    """ A class that implements the running of multiple generic jobs on
    a series of files without collisions with other nodes using the
    ComputeTask class. Jobs that have been started are listed in
    `jobs.txt`, failed jobs are moved to `bad_castep/`, completed jobs
    are moved to `completed/`.

    Interface initially inspired by on run.pl, run2.pl and PyAIRSS class
    CastepRunner.

    """

    def __init__(self, seed, **kwargs):
        """ Check directory has valid contents and prepare log files
        and directories if not already prepared, then begin running
        calculations.

        Note:
            This class is usually initialised by the run3 script, which
            has a full description of possible arguments.

        Parameters:
            seed (:obj:`list` of :obj:`str`): single entry of param/cell
                file seed for CASTEP geometry optimisations of res
                files, or a list of filenames of $seed to run arbitrary
                executables on. e.g. ['LiAs'] if LiAs.cell and LiAs.param
                exist in cwd full of res files, e.g.2. ['LiAs_1', 'LiAs_2']
                if LiAs_1.in/LiAs_2.in exist, and executable = 'pw6.x < $seed.in'.

        Keyword arguments:
            Exhaustive list found in argparse parser inside `matador/cli/run3.py`.

        """
        # parse args, then co-opt them for passing directly into ComputeTask
        prop_defaults = {'ncores': None, 'nprocesses': 1, 'nnodes': 1,
                         'executable': 'castep', 'no_reopt': False, 'mode': None,
                         'redirect': None, 'debug': False, 'custom_params': False,
                         'verbosity': 0, 'archer': False, 'slurm': False,
                         'intel': False, 'conv_cutoff': False, 'conv_kpt': False,
                         'memcheck': False, 'maxmem': None, 'killcheck': True, 'scratch_prefix': None,
                         'kpts_1D': False, 'spin': None, 'ignore_jobs_file': False,
                         'rough': 4, 'rough_iter': 2, 'fine_iter': 20, 'max_walltime': None,
                         'limit': None, 'profile': False, 'polltime': 30}
        self.args = {}
        self.args.update(prop_defaults)
        self.args.update(kwargs)
        self.debug = self.args.get('debug')

        self.seed = seed
        # if only one seed, check if it is a file, and if so treat
        # this run as a generic run, not a CASTEP cell/param run
        if len(self.seed) == 1 and isinstance(self.seed, list):
            if '*' in self.seed[0]:
                self.seed = glob.glob(self.seed[0])
            elif not os.path.isfile(self.seed[0]):
                self.seed = self.seed[0]

        self.compute_dir = os.uname()[1]
        if self.args.get('scratch_prefix') not in [None, '.']:
            self.compute_dir = '{}/{}'.format(self.args['scratch_prefix'], self.compute_dir).replace('//', '/')
        elif self.args.get('scratch_prefix') == '.':
            self.compute_dir = None

        if self.args.get('mode') is not None:
            self.mode = self.args.get('mode')
        else:
            if isinstance(self.seed, str):
                self.mode = 'castep'
            else:
                self.mode = 'generic'
        del self.args['mode']

        if self.args.get('no_reopt'):
            self.args['reopt'] = False
        else:
            self.args['reopt'] = True
        if 'no_reopt' in self.args:
            del self.args['no_reopt']
        self.nprocesses = int(self.args['nprocesses'])
        del self.args['nprocesses']
        self.limit = self.args.get('limit')
        del self.args['limit']
        self.maxmem = self.args.get('maxmem')
        del self.args['maxmem']
        self.max_walltime = self.args.get('max_walltime')

        # detect and scrape queue settings
        self.queue_mgr = get_queue_manager()

        if self.queue_mgr is not None:
            if self.maxmem is None:
                self.maxmem = self.queue_mgr.max_memory
            if self.max_walltime is None:
                self.max_walltime = self.queue_mgr.walltime

        self.start_time = None
        if self.max_walltime is not None:
            self.start_time = time.time()

        # assign number of cores
        self.all_cores = psutil.cpu_count(logical=False)
        if self.args.get('ncores') is None:
            if self.queue_mgr is None:
                self.args['ncores'] = int(self.all_cores / self.nprocesses)
            else:
                self.args['ncores'] = int(self.queue_mgr.ntasks / self.nprocesses)

        if self.args['nnodes'] < 1 or self.args['ncores'] < 1 or self.nprocesses < 1:
            raise InputError('Invalid number of cores, nodes or processes.')
        if self.all_cores < self.nprocesses:
            raise InputError('Requesting more processes than available cores: {} vs {}'
                             .format(self.all_cores, self.nprocesses))

        # scrape input cell/param/other files
        if self.mode == 'castep':
            self.castep_setup()
        else:
            self.generic_setup()

        # prepare folders and text files
        self.paths = dict()
        if self.args.get('conv_cutoff'):
            self.paths['completed_dir'] = 'completed_cutoff'
        elif self.args.get('conv_kpt'):
            self.paths['completed_dir'] = 'completed_kpts'
        else:
            self.paths['completed_dir'] = 'completed'
        self.paths['failed_dir'] = 'bad_castep'
        self.paths['jobs_fname'] = 'jobs.txt'
        self.paths['completed_fname'] = 'finished_cleanly.txt'
        self.paths['failures_fname'] = 'failures.txt'
        self.paths['memory_fname'] = 'memory_exceeded.txt'
        if not os.path.isfile(self.paths['jobs_fname']):
            with open(self.paths['jobs_fname'], 'a'):
                pass
        if not os.path.isfile(self.paths['completed_fname']):
            with open(self.paths['completed_fname'], 'a'):
                pass
        if not os.path.isfile(self.paths['failures_fname']):
            with open(self.paths['failures_fname'], 'a'):
                pass
        if self.args.get('memcheck'):
            if not os.path.isfile(self.paths['memory_fname']):
                with open(self.paths['memory_fname'], 'a'):
                    pass

    def spawn(self, join=False):
        """ Spawn processes to perform calculations.

        Keyword arguments:
            join (bool): whether or not to attach to ComputeTask
                process. Useful for testing.

        """
        procs = []
        error_queue = mp.Queue()
        for proc_id in range(self.nprocesses):
            procs.append(
                mp.Process(
                    target=self.perform_new_calculations,
                    args=(
                        random.sample(self.file_lists['res'],
                                      len(self.file_lists['res'])),
                        error_queue,
                        proc_id
                    )
                )
            )
        for proc in procs:
            proc.start()
            if join:
                proc.join()

        errors = []
        failed_seeds = []

        # wait for each proc to write to error queue
        try:
            for _, proc in enumerate(procs):
                result = error_queue.get()
                if isinstance(result[1], Exception):
                    errors.append(result)
                    failed_seeds.append(result[2])

            if errors:
                error_message = ''
                for error in errors:
                    error_message += 'Process {} raised error(s): {}. '.format(error[0], error[1])
                    if len({type(error[1]) for error in errors}) == 1:
                        raise errors[0][1]
                    raise type(errors[0][1])(error_message)
                raise BundledErrors(error_message)
        # the only errors that reach here are fatal, e.g. WalltimeError, CriticalError, InputError, KeyboardInterrupt
        except RuntimeError as err:
            result = [proc.join(timeout=2) for proc in procs]
            result = [proc.terminate() for proc in procs if proc.is_alive()]
            print_failure('Fatal error(s) reported:')
            print_warning(err)
            raise err

        print('Nothing left to do.')

    def perform_new_calculations(self, res_list, error_queue, proc_id):
        """ Perform all calculations that have not already
        failed or finished to completion.

        Parameters:
            res_list (:obj:`list` of :obj:`str`): list of structure filenames.
            error_queue (multiprocessing.Queue): queue to push exceptions to
            proc_id (int): process id for logging

        """
        job_count = 0
        if isinstance(res_list, str):
            res_list = [res_list]
        for res in res_list:
            try:
                if not os.path.isfile(res):
                    continue
                # probe once then sleep for a random amount up to 5 seconds
                # before checking again for a lock file, just to protect
                # against collisions in large array jobs on slower parallel file systems
                _ = os.path.isfile('{}.lock'.format(res))
                time.sleep(2 * random.random())
                # wait some additional time if this is a slurm array job
                if self.queue_mgr is not None:
                    extra_wait = self.queue_mgr.array_id % 10 if self.queue_mgr.array_id else 0
                    time.sleep(extra_wait)
                locked = os.path.isfile('{}.lock'.format(res))
                if not self.args.get('ignore_jobs_file'):
                    listed = self._check_jobs_file(res)
                else:
                    listed = []
                running = any([listed, locked])
                if not running:
                    # check we haven't reached job limit
                    if self.limit is not None and job_count >= self.limit:
                        error_queue.put((proc_id, job_count, res))
                        return

                    # check 3 more times if a lock exists with random up to 1 second
                    # waits each time
                    for _ in range(3):
                        time.sleep(random.random())
                        if os.path.isfile('{}.lock'.format(res)):
                            raise NodeCollisionError('Another node wrote this file when I wanted to, skipping...')
                    with open(res + '.lock', 'a') as job_file:
                        pass

                    # write to jobs file
                    with open(self.paths['jobs_fname'], 'a') as job_file:
                        job_file.write(res + '\n')

                    # create full relaxer object for creation and running of job
                    job_count += 1
                    relaxer = ComputeTask(node=None, res=res,
                                          param_dict=self.param_dict,
                                          cell_dict=self.cell_dict,
                                          mode=self.mode, paths=self.paths, compute_dir=self.compute_dir,
                                          timings=(self.max_walltime, self.start_time), maxmem=self.maxmem,
                                          **self.args)
                    # if memory check failed, let other nodes have a go
                    if not relaxer.enough_memory:
                        with open(self.paths['memory_fname'], 'a') as job_file:
                            job_file.write(res + '\n')
                        if os.path.isfile('{}.lock'.format(res)):
                            os.remove('{}.lock'.format(res))
                        with open(self.paths['jobs_fname'], 'r+') as job_file:
                            flines = job_file.readlines()
                            job_file.seek(0)
                        for line in flines:
                            if res not in line:
                                job_file.write(line)
                            job_file.truncate()

                    elif relaxer.success:
                        with open(self.paths['completed_fname'], 'a') as job_file:
                            job_file.write(res + '\n')
                    else:
                        with open(self.paths['failures_fname'], 'a') as job_file:
                            job_file.write(res + '\n')

            # catch memory errors and reset so another node can try
            except MaxMemoryEstimateExceeded:
                reset_single_seed(res)
                continue

            # ignore any other individual calculation errors or node collisions that were caught here
            except CalculationError:
                continue

            # reset txt/lock for an input error, but throw it to prevent other calcs
            except InputError as err:
                reset_single_seed(res)
                error_queue.put((proc_id, err, res))
                return
            # push globally-fatal errors to queue, and return to prevent further calcs
            except RuntimeError as err:
                error_queue.put((proc_id, err, res))
                return
            # finally catch any other generic error in e.g. the above code, normally caused by me
            except Exception as err:
                error_queue.put((proc_id, err, res))
                return

        error_queue.put((proc_id, job_count, ''))

    def generic_setup(self):
        """ Undo things that are set ready for CASTEP jobs... """
        self.cell_dict = None
        self.param_dict = None

        # scan directory for files to run
        self.file_lists = defaultdict(list)
        self.file_lists['res'] = glob.glob('*.res')

    def castep_setup(self):
        """ Set up CASTEP jobs from res files, and $seed.cell/param. """
        # read cell/param files
        exts = ['cell', 'param']
        for ext in exts:
            if not os.path.isfile('{}.{}'.format(self.seed, ext)):
                raise InputError('Failed to find {ext} file, {seed}.{ext}'.format(ext=ext, seed=self.seed))
        self.cell_dict, cell_success = cell2dict(self.seed + '.cell',
                                                 db=False, lattice=False, positions=True)
        if not cell_success:
            print(self.cell_dict)
            raise InputError('Failed to parse cell file')
        self.param_dict, param_success = param2dict(self.seed + '.param', db=False)
        if not param_success:
            print(self.param_dict)
            raise InputError('Failed to parse param file')

        # scan directory for files to run
        self.file_lists = defaultdict(list)
        self.file_lists['res'] = glob.glob('*.res')
        if any(self.seed == file.replace('.res', '') for file in self.file_lists['res']):
            error = ("Found .res file with same name as seed: {}.res. This will wreak havoc on your calculations!\n"
                     .format(self.seed)
                     + "Please rename either your seed.cell/seed.param files, or rename the offending {}.res"
                     .format(self.seed))
            raise InputError(error)

        if not self.file_lists['res']:
            error = (
                'run3 in CASTEP mode requires at least 1 res file in folder, found {}'
                .format(len(self.file_lists['res']))
            )
            raise InputError(error)

        if (len(self.file_lists['res']) < self.nprocesses
                and not any([self.args.get('conv_cutoff'), self.args.get('conv_kpt')])):
            raise InputError('Requested more processes than there are jobs to run!')

        # do some prelim checks of parameters
        if self.param_dict['task'].upper() in ['GEOMETRYOPTIMISATION', 'GEOMETRYOPTIMIZATION']:
            if 'geom_max_iter' not in self.param_dict:
                raise InputError('geom_max_iter is unset, please fix this.')
            if int(self.param_dict['geom_max_iter']) <= 0:
                raise InputError('geom_max_iter is only {}!'.format(self.param_dict['geom_max_iter']))

        # parse convergence args and set them up
        self.convergence_run_setup()

        # delete source from cell and param
        del self.cell_dict['source']
        del self.param_dict['source']

    def convergence_run_setup(self):
        """ Set the correct args for a convergence run. """
        # check if we're doing a conv run
        if self.args.get('conv_cutoff'):
            if os.path.isfile('cutoff.conv'):
                with open('cutoff.conv', 'r') as f:
                    flines = f.readlines()
                    self.args['conv_cutoff'] = []
                    for line in flines:
                        if not line.startswith('#'):
                            self.args['conv_cutoff'].append(int(line))
            else:
                raise InputError('Missing cutoff.conv file')
        else:
            self.args['conv_cutoff'] = None

        if self.args.get('conv_kpt'):
            if os.path.isfile('kpt.conv'):
                with open('kpt.conv', 'r') as f:
                    flines = f.readlines()
                    self.args['conv_kpt'] = []
                    for line in flines:
                        if not line.startswith('#'):
                            self.args['conv_kpt'].append(float(line))
            else:
                raise InputError('Missing with conv.kpt file')
        else:
            self.args['conv_kpt'] = None

    def _check_jobs_file(self, res):
        """ Check if structure is listed in jobs.txt file.

        Parameters:
            res (str): structure name.

        Returns:
            bool: True if already listed in jobs file.

        """
        with open(self.paths['jobs_fname'], 'r') as job_file:
            flines = job_file.readlines()
            for line in flines:
                if res in line:
                    return True
        return False


class BundledErrors(Exception):
    """ Raise this after collecting all exceptions from
    processes.
    """


def reset_job_folder(debug=False):
    """ Remove all lock files and clean up jobs.txt
    ready for job restart.

    Note:
        This should be not called by a ComputeTask instance, in case
        other instances are running.

    Returns:
        num_remaining (int): number of structures left to relax

    """
    res_list = glob.glob('*.res')
    if debug:
        print(res_list)
    for f in res_list:
        root = f.replace('.res', '')
        exts_to_rm = ['res.lock', 'kill']
        for ext in exts_to_rm:
            if os.path.isfile('{}.{}'.format(root, ext)):
                if debug:
                    print('Deleting {}.{}'.format(root, ext))
                os.remove('{}.{}'.format(root, ext))

    # also remove from jobs file
    if os.path.isfile('jobs.txt'):
        with open('jobs.txt', 'r+') as f:
            flines = f.readlines()
            if debug:
                print('Initially {} jobs in jobs.txt'.format(len(flines)))
            f.seek(0)
            for line in flines:
                line = line.strip()
                if line in res_list:
                    print('Excluding {}'.format(line))
                    continue
                f.write(line)
            f.truncate()
            flines = f.readlines()
            if debug:
                print('{} jobs remain in jobs.txt'.format(len(flines)))

    return len(res_list)


def reset_single_seed(seed):
    """ Remove the file lock and jobs.txt entry
    for a single seed.

    Parameters:
        seed (str): the seedname to remove.

    """
    if os.path.isfile('jobs.txt'):
        with open('jobs.txt', 'r+') as f:
            flines = f.readlines()
            f.seek(0)
            for line in flines:
                line = line.strip()
                if seed in line:
                    continue
                f.write(line)
            f.truncate()
            flines = f.readlines()
    if os.path.isfile(seed + '.lock'):
        os.remove(seed + '.lock')
