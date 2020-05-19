# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the :class:`ComputeTask` class for handling
calculations on a single structure.

"""

import os
import sys
import shutil
import subprocess as sp
import glob
import time
import logging
import traceback as tb
from copy import deepcopy
from math import ceil
from psutil import virtual_memory

from matador import __version__
from matador.scrapers.castep_scrapers import cell2dict
from matador.scrapers.castep_scrapers import res2dict, castep2dict
from matador.calculators import CastepCalculator
from matador.crystal import Crystal
from matador.export import doc2cell, doc2param, doc2res
from matador.utils.errors import CriticalError, WalltimeError, InputError, CalculationError, MaxMemoryEstimateExceeded
from matador.utils.chem_utils import get_root_source
import matador.utils.print_utils

MATADOR_CUSTOM_TASKS = ['bulk_modulus', 'projected_bandstructure', 'pdispersion', 'all']

LOG = logging.getLogger('run3')
LOG.setLevel(logging.DEBUG)

__all__ = ["ComputeTask"]


class ComputeTask:
    """ The main use of this class is to call an executable on a given
    structure. The various parameters are passed to this class by the
    common entrypoints, run3 and ilustrado. It is unlikely that you
    will want to use this class directly. Each keyword is saved as
    an attribute for later use.

    Attributes:
        self.final_result (dict): stores the final result of the
            calculation, if it was successful.

    Note:
        By default, calculations are run inside a folder with the same
        name as the host (e.g. node12, or whatever). This decreases the
        load on parallel file systems such as Lustre.

    """

    def __init__(self, res, ncores, nnodes, node, **kwargs):
        """ Make the files to run the calculation and call the desired program.

        Parameters:
            res (str/dict): filename or input structure dict
            ncores (int): number of cores *per node* for mpirun call
            nnodes (int): number of nodes for mpirun call (if None, use 1)
            node (str): node name to run on (if None, run on localhost)

        Keyword arguments:
            param_dict (dict): dictionary of CASTEP parameters
            cell_dict (dict): dictionary of CASTEP cell options
            executable (str): name of binary to execute (DEFAULT: 'castep').
                Special string $seed will be parsed as the seedname,
                e.g. executable = 'pw6.x -i $seed.in > $seed.out' (requires mode='generic').
            mode (str): either 'castep' or 'generic' (DEFAULT: 'castep')
            custom_params (bool): use custom per-structure param file
                (DEFAULT: False)
            output_queue (multiprocessing.Queue): write results to queue rather than file.
            rough (int): number of small "rough" calculations (DEFAULT: 4)
            rough_iter (int): number of iterations per rough calculation
                (DEFAULT: 2)
            fine_iter (int): number of iterations per fine calculation
                (DEFAULT: 20)
            spin (bool): break spin symmetry in first calculation by amount specified
                (DEFAULT: None if not present, 5 if no argument)
            conv_cutoffs (:obj:`list` of :obj:`float`): list of cutoffs
                to use for SCF convergence test
            conv_kpts (:obj:`list` of :obj:`float`): list of kpt spacings
                to use for SCF convergence test
            kpts_1D (bool): treat z-direction as special and create
                kpt_grid [1 1 n_kz] (DEFAULT: False)
            noise (bool): add noise to the positions (DEFAULT: False)
            squeeze (bool/float): add an external pressure to the first steps (DEFAULT: False)
            archer (bool): force use of aprun over mpirun (DEFAULT: False)
            slurm (bool): force use of srun over mpirun (DEFAULT: False)
            intel (bool): force use of Intel mpirun-style calls (DEFAULT: False)
            redirect (str): file to redirect stdout to (DEFAULT: /dev/null unless debug).
            exec_test (bool): test executable with `<exec> --version`
                before progressing (DEFAULT: True)
            start (bool): begin calculation immediately or manually
                call it (DEFAULT: True)
            reopt (bool): whether to optimise one more time after
                success (DEFAULT: False)
            memcheck (bool): perform CASTEP dryrun to estimate memory
                usage, do not proceed if fails (DEFAULT: False)
            maxmem (int): maximum memory allowed in MB for memcheck
                (DEFAULT: None)
            killcheck (bool): check for file called $seed.kill during
                operation, and kill executable if present (DEFAULT: True)
            compute_dir (str): folder to run computations in; default is
                None (i.e. cwd), if not None, prepend paths with this folder
            verbosity (int): either 0, 1, 2 or >3, corresponding to ERROR, WARNING
                INFO and DEBUG logging levels.
            timings (`obj`:tuple: of `obj`:int:): tuple containing max and elapsed time in seconds

        Raises:
            WalltimeError: if desired/alotted walltime is exceeded, current run will be
                tidied up, ready to be restarted from intermediate state.
            CriticalError: if a fatal error occurs, failed run will be moved to bad_castep and
                no further calculations will be attempted.
            CalculationError: if a structure-level error occurs, causing the seed files to be moved
                to bad_castep.

        """
        # set defaults and update class with desired values
        prop_defaults = {'paths': None, 'param_dict': None, 'cell_dict': None, 'mode': 'castep', 'executable': 'castep',
                         'memcheck': False, 'rough': 4, 'rough_iter': 2, 'fine_iter': 20, 'spin': None, 'squeeze': False,
                         'output_queue': None, 'redirect': None, 'reopt': False, 'compute_dir': None, 'noise': False,
                         'custom_params': False, 'archer': False, 'maxmem': None, 'killcheck': True, 'kpts_1D': False,
                         'conv_cutoff': False, 'conv_kpt': False, 'slurm': False, 'intel': False,
                         'exec_test': True, 'timings': (None, None), 'start': True, 'verbosity': 1, 'polltime': 60,
                         'optados_executable': 'optados', 'run3_settings': dict(), 'workflow_kwargs': dict()}

        self.paths = None
        self.output_queue = None
        self.final_result = None
        self.executable = None

        self._first_run = True
        self._geom_max_iter_list = None
        self._squeeze_list = None
        self._max_iter = None
        self._num_rough_iter = None
        self._target_spacing = None
        self._target_pressure = None
        self._num_retries = 0
        self._max_num_retries = 2
        self.maxmem = None
        self.cell_dict = None
        self.param_dict = None
        self.res_dict = None
        self.calc_doc = None

        # save all keyword arguments as attributes
        self.__dict__.update(prop_defaults)
        for arg in kwargs:
            if arg in prop_defaults:
                self.__dict__.update({arg: kwargs[arg]})

        if self.maxmem is None and self.memcheck:
            self.maxmem = float(virtual_memory().available) / 1024**2

        # scrape and save seed name
        if isinstance(res, str):
            self.seed = res
            self.seed = self.seed.split('/')[-1].replace('.res', '')
        else:
            self.seed = get_root_source(res)
        self._input_res = res

        # set up logging; by default, write DEBUG level into a file called logs/$seed.log
        # and write desired verbosity to stdout (default WARN).
        if self.verbosity == 0:
            loglevel = logging.ERROR
        elif self.verbosity == 1:
            loglevel = logging.WARN
        elif self.verbosity == 2:
            loglevel = logging.INFO
        else:
            loglevel = logging.DEBUG

        if not os.path.isdir('logs'):
            os.mkdir('logs')

        LOG.handlers = []

        if self.verbosity > 1:
            stdout_handler = logging.StreamHandler(sys.stdout)
            stdout_handler.setLevel(loglevel)
            stdout_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s | %(levelname)8s: %(message)s'))
            LOG.addHandler(stdout_handler)

        logname = os.path.abspath('logs/{}.log'.format(self.seed))
        file_handler = logging.FileHandler(logname, mode='a')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s | %(levelname)8s: %(message)s'))
        LOG.addHandler(file_handler)

        splash_screen = r"""
                           _____
          _ __ _   _ _ __ |___ /
         | '__| | | | '_ \  |_ \
         | |  | |_| | | | |___) |
         |_|   \__,_|_| |_|____/

        """

        LOG.info(splash_screen)
        LOG.info('matador version: {}'.format(__version__))
        LOG.info('Host info: {}.'.format('-'.join(os.uname())))
        LOG.info('Run started for seed {}'.format(self.seed))

        # set up compute parameters
        self.ncores = ncores
        self.nnodes = nnodes
        self.node = node
        self.conv_cutoff_bool = isinstance(self.conv_cutoff, list)
        self.conv_kpt_bool = isinstance(self.conv_kpt, list)
        self.enough_memory = True
        self.success = None
        self._mpi_library = None
        self.root_folder = os.getcwd()

        if self.timings is not None:
            self.max_walltime = self.timings[0]
            self.start_time = self.timings[1]

        if self.max_walltime is not None:
            LOG.debug('Starting at {start}, max walltime allowed is {walltime}'.format(
                start=self.start_time, walltime=self.max_walltime))
            LOG.debug('{} s remaining.'.format(self.max_walltime - (time.time() - self.start_time)))

        self._redirect_filename = None

        if self.paths is None:
            self.paths = {}
            self.paths['completed_dir'] = 'completed'
            self.paths['failed_dir'] = 'bad_castep'
            self.paths['memory_fname'] = 'memory_exceeded.txt'
            self.paths['jobs_fname'] = 'jobs.txt'
            self.paths['completed_fname'] = 'finished_cleanly.txt'
            self.paths['failures_fname'] = 'failures.txt'

        elif 'completed_dir' not in self.paths:
            raise RuntimeError('Invalid paths: {}'.format(self.paths))

        self._set_input_structure(res)

        if self.start:
            self.begin()
        else:
            LOG.info('Waiting for `begin()` method to be called...')

    def _set_input_structure(self, res):
        """ Set the input structure to the given dictionary, or read it from
        the given file.

        Parameters:
            res (str/dict): either the scraped structure or the filename of a res file.

        """

        # read in initial structure and skip if failed
        if isinstance(res, str):
            self.res_dict, success = res2dict(res, db=False, verbosity=self.verbosity)
            if not success:
                msg = 'Unable to parse initial res file, error: {res_dict}'.format(res_dict=self.res_dict)
                raise CalculationError(msg)

        elif isinstance(res, dict):
            # check res is a valid crystal
            try:
                Crystal(res)
            except Exception as exc:
                msg = 'Unable to convert res to Crystal: {}'.format(exc)
                LOG.error(msg)
                LOG.error('Traceback:\n{}'.format(tb.format_exc()))
                raise CalculationError(msg)

            self.res_dict = res

    def begin(self):
        """ Run the prepared ComputeTask. Catches CalculationError objects
        and cleans up, passing all other errors upwards.

        """

        try:
            try:
                # run through CASTEP specific features
                if self.mode == 'castep':
                    calculation_parameters = {**self.cell_dict, **self.param_dict}
                    self.calculator = CastepCalculator(calculation_parameters)
                    self.success = self.run_castep()
                # otherwise run generic script
                else:
                    self.success = self.run_generic()

            except CalculationError as exc:
                self.success = False
                self._finalise_result()
                raise exc

        except Exception as exc:
            LOG.error('Process raised {} with message {}.'.format(type(exc), exc))
            if isinstance(exc, CalculationError):
                LOG.error('Full traceback:\n{}'.format(tb.format_exc()))
            for handler in LOG.handlers[:]:
                handler.close()
            raise exc

        if self.success:
            LOG.info('ComputeTask finished successfully for {seed}'.format(seed=self.seed))
        else:
            LOG.info('ComputeTask failed cleanly for {seed}'.format(seed=self.seed))
        for handler in LOG.handlers[:]:
            handler.close()

    def run_castep(self):
        """ Set up and run CASTEP calculation on the prepared structure,
        `self.res_dict`, using the parameters in `self.cell_dict` and
        `self.param_dict`.

        Raises:
            WalltimeError: if max_walltime is exceeded.
            CriticalError: if no further calculations should be performed
                on this thread.
            CalculationError: if this structure errored in some way, but
                others will hopefully be okay.

        Returns:
            bool: True if calculations were successful, False otherwise.
                In the case of convergence tests, this is always True
                unless every calculation fails.

        """
        success = False

        if self.exec_test:
            self.test_exec()

        # pre-processing
        if self.kpts_1D:
            LOG.debug('1D kpoint grid requested.')
            if 'kpoints_mp_spacing' not in self.cell_dict:
                raise CriticalError('kpoints_mp_spacing not found, but kpts_1D requested...')
            self._target_spacing = deepcopy(self.cell_dict['kpoints_mp_spacing'])

        if self.squeeze:
            self._target_pressure = deepcopy(self.cell_dict['external_pressure'])

        if self.noise:
            from matador.utils.cell_utils import add_noise
            LOG.info('Adding noise to the positions in the cell')
            self.res_dict = add_noise(self.res_dict, amplitude=0.1)

        self.calc_doc = deepcopy(self.res_dict)

        # update global doc with cell and param dicts for folder, and verify
        bad_keys = ['atom_types', 'positions_frac', 'positions_abs', 'lattice_cart', 'lattice_abc']
        self.cell_dict = {key: self.cell_dict[key] for key in self.cell_dict if key not in bad_keys}
        self.calc_doc.update(self.cell_dict)
        self.calc_doc.update(self.param_dict)
        self.calculator.verify_calculation_parameters(self.calc_doc, self.res_dict)

        try:
            LOG.info('Struture: {}'.format(Crystal(self.calc_doc)))
        except Exception as exc:
            LOG.warning('Unable to convert structure to Crystal... {}'.format(exc))

        # now verify the structure itself
        self.calculator.verify_simulation_cell(self.res_dict)

        # do memcheck, if desired, and only continue if enough memory is free
        if self.memcheck:
            memory_usage_estimate = self.do_memcheck(self.calc_doc, self.seed)
            mem_string = ('Memory estimate / Available memory (MB): {:8.0f} / {:8.0f}'
                          .format(memory_usage_estimate, self.maxmem))
            LOG.info(mem_string)

            if memory_usage_estimate > 0.9 * self.maxmem:
                msg = 'Structure {} failed memcheck, skipping... '.format(self.seed) + mem_string
                with open(self.paths['memory_fname'], 'a') as f:
                    f.write("{} {:.2f}GB/{:.2f}GB\n".format(self.seed, memory_usage_estimate / 1024, self.maxmem / 1024))

                raise MaxMemoryEstimateExceeded(msg)

        try:
            # run convergence tests
            if any([self.conv_cutoff_bool, self.conv_kpt_bool]):
                success = self.run_convergence_tests(self.calc_doc)

            # perform relaxation
            elif self.calc_doc['task'].upper() in ['GEOMETRYOPTIMISATION', 'GEOMETRYOPTIMIZATION']:
                success = self.run_castep_relaxation()

            elif self.calc_doc['task'].upper() in ['PHONON', 'THERMODYNAMICS']:
                from matador.workflows.castep import castep_full_phonon
                success = castep_full_phonon(self, self.calc_doc, self.seed, **self.workflow_kwargs)

            elif self.calc_doc['task'].upper() in ['SPECTRAL']:
                from matador.workflows.castep import castep_full_spectral
                success = castep_full_spectral(self, self.calc_doc, self.seed, **self.workflow_kwargs)

            elif self.calc_doc['task'].upper() in ['BULK_MODULUS']:
                from matador.workflows.castep import castep_elastic
                success = castep_elastic(self, self.calc_doc, self.seed, **self.workflow_kwargs)

            # run in singleshot mode, i.e. just call CASTEP on the seeds
            else:
                success = self.run_castep_singleshot(self.calc_doc, self.seed, keep=True)

            self._first_run = False

            return success

        finally:
            # always cd back to root folder, in case we have ended up on the wrong place
            if self.compute_dir is not None:
                os.chdir(self.root_folder)
                self.remove_compute_dir_if_finished(self.compute_dir)

    def run_generic(self, intermediate=False, mv_bad_on_failure=True):
        """ Run a generic mpi program on the given seed. Files from
        completed runs are moved to "completed" (unless intermediate
        is True) and failed runs to "bad_castep".

        Keyword arguments:
            intermediate (bool): whether we want to run more calculations
                on the output of this, i.e. whether to move to completed
                or not.
            mv_bad_on_failure (bool): whether to move files to bad_castep on
                failure, or leave them in place.

        Returns:
            bool: True if calculations progressed without error.

        """
        LOG.info('Calling executable {exe} MPI program on {seed}'.format(exe=self.executable, seed=self.seed))
        try:
            seed = self._input_res
            if '.' in seed:
                seed = seed.split('.')[-2]
                input_ext = seed.split('.')[-1]
            else:
                input_ext = ''

            assert isinstance(seed, str)
            self.cp_to_input(seed, ext=input_ext, glob_files=True)

            self._setup_compute_dir(self.seed, self.compute_dir, generic=True)
            if self.compute_dir is not None:
                os.chdir(self.compute_dir)

            _process = self.run_command(seed)
            self._handle_process(_process, check_walltime=True)

            if not intermediate:
                LOG.info('Writing results of generic call to res file and tidying up.')
                self.mv_to_completed(seed, keep=True, completed_dir=self.paths['completed_dir'])

            LOG.info('Executable {exe} finished cleanly.'.format(exe=self.executable))
            self._first_run = False

            return True

        except Exception as err:
            self._first_run = False
            LOG.error('Caught error inside run_generic: {error}.'.format(error=err))
            if mv_bad_on_failure:
                self.mv_to_bad(seed)
            raise err

        finally:
            if self.compute_dir is not None:
                os.chdir(self.root_folder)

    def run_castep_relaxation(self, intermediate=False):
        """ Set up a structural relaxation that is restarted intermittently
        in order to re-mesh the kpoint grid. Completed calculations are moved
        to the "completed" folder, and failures to "bad_castep".

        Keyword arguments:
            intermediate (bool): whether we want to run more calculations
                on the output of this, i.e. whether to move to completed
                or not.

        Returns:
            bool: True iff structure was optimised, False otherwise.

        Raises:
            CalculationError: if structure-level error occured.
            CriticalError: if fatal global error occured.
            WalltimeError: if walltime was reached, and jobs need to stop.

        """
        try:
            return self._relax(intermediate=intermediate)

        except WalltimeError as err:
            raise err

        except Exception as err:
            self._finalise_result(intermediate)
            raise err

        finally:
            if self.compute_dir is not None:
                os.chdir(self.root_folder)

    def _relax(self, intermediate=False):
        """ Set up a structural relaxation that is restarted intermittently
        in order to re-mesh the kpoint grid. Completed calculations are moved
        to the "completed" folder, and failures to "bad_castep".

        Keyword arguments:
            intermediate (bool): whether we want to run more calculations
                on the output of this, i.e. whether to move to completed
                or not.

        Returns:
            bool: True iff structure was optimised, False otherwise.

        Raises:
            CalculationError: if structure-level error occured.
            CriticalError: if fatal global error occured.
            WalltimeError: if walltime was reached, and jobs need to stop.

        """
        self._setup_relaxation()
        seed = self.seed
        rerun = False
        # iterate over geom iter blocks
        for ind, num_iter in enumerate(self._geom_max_iter_list):

            # preprocess last step that did not finish geom opt
            if self.reopt and rerun:
                # if we're now reoptimising after a success in last step
                # use the fine iter value
                num_iter = self.fine_iter
                LOG.info('Last step was successful, performing one last relaxation...')

            # update the geom_max_iter to use with either the number in iter_list, or the overriden value
            self.calc_doc['geom_max_iter'] = num_iter

            # delete any existing files and write new ones
            if self.squeeze:
                squeeze = int(self._squeeze_list[ind]) * float(self.squeeze)
            else:
                squeeze = None
            self._update_input_files(self.seed, self.calc_doc, squeeze=squeeze)

            if self.max_walltime is not None and self.start_time is None:
                msg = 'Somehow initial start time was not found'
                LOG.critical(msg)
                raise CriticalError(msg)

            # run CASTEP
            _process = self.run_command(seed)

            # will throw errors if the process fails
            output_filename = "{}.{}".format(seed, 'castep')

            self._handle_process(
                _process, expected_fname=output_filename, check_walltime=True
            )

            # check for errors and try to correct for them
            errors_present, errors, remedy = self._catch_castep_errors(_process)
            skip_postprocess = remedy is not None

            # try to read the CASTEP file
            opti_dict, success = castep2dict(output_filename, db=False, verbosity=self.verbosity)

            if errors_present:
                msg = 'Failed to optimise {} as CASTEP crashed with error:'.format(seed)
                msg += errors
                LOG.warning(msg)
                if isinstance(opti_dict, dict):
                    self._update_castep_output_files(seed, opti_dict)
                if remedy is not None and self._num_retries <= self._max_num_retries:
                    LOG.warning('Attempting to recover using {}'.format(remedy))
                else:
                    raise CalculationError(msg)

            if not success and remedy is None and isinstance(opti_dict, Exception):
                msg = 'Failed to parse CASTEP file... {}'.format(opti_dict)
                LOG.warning(msg)
                raise CalculationError(msg)
            if isinstance(opti_dict, Exception) and remedy is not None:
                opti_dict = {'optimised': False}

            LOG.info('Intermediate calculation completed successfully: total relaxation steps now = {} '
                     .format(opti_dict.get('geom_iter')))

            # scrub keys that need to be rescraped
            keys_to_remove = ['kpoints_mp_spacing', 'kpoints_mp_grid', 'species_pot',
                              'sedc_apply', 'sedc_scheme', 'cell_constraints']
            for key in keys_to_remove:
                if key in opti_dict:
                    del opti_dict[key]

            if not skip_postprocess:
                # now post-process the last step
                # did we try to rerun, and are now no longer optimised?
                # then reset, and go again...
                if self.reopt and rerun and not opti_dict['optimised']:
                    rerun = False
                    self._update_castep_output_files(seed, opti_dict)

                # are we optimised, but haven't yet rerun?
                # then set prepare to do one more full relaxation
                if self.reopt and not rerun and opti_dict['optimised']:
                    rerun = True

                    # disable squeezing if we've already reached optimisation
                    if self.squeeze:
                        for jnd in range(ind, len(self._squeeze_list)):
                            self._squeeze_list[jnd] = False
                    self._update_castep_output_files(seed, opti_dict)

                # or did the relaxation complete successfuly, including rerun?
                elif (not self.reopt or rerun) and opti_dict['optimised']:
                    LOG.info('Successfully relaxed {}'.format(seed))
                    self._update_castep_output_files(seed, opti_dict)
                    break

                # reached maximum number of steps
                elif ind == len(self._geom_max_iter_list) - 1:
                    msg = 'Failed to optimise {} after {} steps'.format(seed, sum(self._geom_max_iter_list))
                    LOG.info(msg)
                    raise CalculationError(msg)

                # if we weren't successful, then preprocess the next step

                # set atomic_init_spins with value from CASTEP file, if it exists
                if 'mulliken_spins' in opti_dict:
                    self.calc_doc['atomic_init_spins'] = opti_dict['mulliken_spins']

                # CASTEP prints cell constraints in castep file when running with symmetry
                # its useful to disable these and use either the initial constraints, or none
                if 'cell_constraints' in self.cell_dict:
                    self.calc_doc['cell_constraints'] = self.cell_dict['cell_constraints']

                # if writing out cell, use it for higher precision lattice_cart
                if self.calc_doc.get('write_cell_structure') and os.path.isfile('{}-out.cell'.format(seed)):
                    cell_dict, success = cell2dict('{}-out.cell'.format(seed),
                                                   verbosity=self.verbosity,
                                                   db=False, positions=True, lattice=True)
                    if success:
                        opti_dict['lattice_cart'] = list(cell_dict['lattice_cart'])

                LOG.info('N = {iters:03d} | |F| = {d[max_force_on_atom]:5.5f} eV/A | '
                         'S = {pressure:5.5f} GPa | H = {enthalpy_per_atom:5.5f} eV/atom'
                         .format(d=opti_dict,
                                 pressure=opti_dict.get('pressure', 0.0),
                                 enthalpy_per_atom=opti_dict.get('enthalpy_per_atom', 0.0),
                                 iters=opti_dict.get('geom_iter', 0)))

            # if there were errors that can be remedied, now is the time to do it
            # this will normally involve changing a parameter to avoid future failures
            self.calc_doc.update(opti_dict)

            if remedy is not None:
                LOG.info('Trying to remedy error with {}'.format(remedy))
                remedy(self.calc_doc)
                self._num_retries += 1

        return self._finalise_result(intermediate=intermediate)

    def run_castep_singleshot(self, calc_doc, seed, keep=True, intermediate=False):
        """ Perform a singleshot calculation with CASTEP. Singleshot runs do not
        attempt to remedy any errors raised.

        Files from completed runs are moved to `completed`, if not
        in intermediate mode, and failed runs to `bad_castep`.

        Parameters:
            calc_doc (dict): dictionary containing parameters and structure
            seed (str): structure filename

        Keyword arguments:
            intermediate (bool): whether we want to run more calculations
                on the output of this, i.e. whether to move to completed
                or not.
            keep (bool): whether to keep intermediate files e.g. .bands

        Returns:
            bool: True iff SCF completed successfully, False otherwise.

        """
        LOG.info('Performing single-shot CASTEP run on {}, with task: {}'.format(seed, calc_doc['task']))
        try:
            self._singleshot(calc_doc, seed, keep=keep, intermediate=intermediate)

        except WalltimeError as err:
            raise err

        except Exception as err:
            self.mv_to_bad(seed)
            if not keep:
                self.tidy_up(seed)
            raise err

        finally:
            if self.compute_dir is not None:
                os.chdir(self.root_folder)

    def _singleshot(self, calc_doc, seed, keep=True, intermediate=False):
        """ Perform a singleshot calculation with CASTEP. Singleshot runs do not
        attempt to remedy any errors raised.

        Files from completed runs are moved to `completed`, if not
        in intermediate mode, and failed runs to `bad_castep`.

        Parameters:
            calc_doc (dict): dictionary containing parameters and structure
            seed (str): structure filename

        Keyword arguments:
            intermediate (bool): whether we want to run more calculations
                on the output of this, i.e. whether to move to completed
                or not.
            keep (bool): whether to keep intermediate files e.g. .bands

        Returns:
            bool: True iff SCF completed successfully, False otherwise.

        """

        self.cp_to_input(seed)
        self._setup_compute_dir(seed, self.compute_dir, custom_params=self.custom_params)
        if self.compute_dir is not None:
            os.chdir(self.compute_dir)

        self._update_input_files(seed, calc_doc)
        doc2res(calc_doc, self.seed, info=False, hash_dupe=False, overwrite=True)

        # run CASTEP
        _process = self.run_command(seed)
        output_filename = "{}.{}".format(seed, 'castep')
        self._handle_process(
            _process, expected_fname=output_filename, check_walltime=True
        )

        # check for errors and report them
        errors_present, errors, _ = self._catch_castep_errors(_process)
        if errors_present:
            raise CalculationError(
                'CASTEP run on {} failed with errors: {}'.format(seed, errors)
            )

        # scrape dict but ignore the results
        results_dict, success = castep2dict(output_filename, db=False, verbosity=self.verbosity)
        if not success:
            raise CalculationError(
                'Error scraping CASTEP file {}: {}'.format(seed, results_dict)
            )

        self._update_castep_output_files(seed)

        if not intermediate:
            LOG.info('Writing results of singleshot CASTEP run to res file and tidying up.')
            doc2res(results_dict, seed, hash_dupe=False, overwrite=True)
            self.mv_to_completed(seed, keep=keep, completed_dir=self.paths['completed_dir'])
            self.tidy_up(seed)

        return success

    @staticmethod
    def validate_calc_doc(calc_doc, required, forbidden):
        """ Remove keys inside forbidden from calc_doc, and error
        if a required key is missing.

        Parameters:
            calc_doc (dict): dictionary of structure and parameters.
            required (list): list of required key strings.
            forbidden (list): list of forbidden keys.

        Raises:
            AssertionError: if required key is missing.

        """
        for keyword in forbidden:
            if keyword in calc_doc:
                del calc_doc[keyword]

        failures = []
        for keyword in required:
            if keyword not in calc_doc:
                failures.append(keyword)
        if failures:
            raise InputError('The following keywords are required for workflow: {}'.format(', '.join(failures)))

    @staticmethod
    def get_seekpath_compliant_input(calc_doc, spacing, debug=False):
        """ Return seekpath cell/kpoint path for the given cell and spacing.

        Parameters:
            calc_doc (dict): structural and calculation parameters.
            spacing (float): desired kpoint path spacing.

        Returns:
            (dict, list): dictionary containing the standardised unit cell
                and list containing the kpoints.

        """
        LOG.info('Getting seekpath cell and kpoint path.')
        from matador.utils.cell_utils import get_seekpath_kpoint_path
        LOG.debug('Old lattice: {}'.format(Crystal(calc_doc)))
        prim_doc, kpt_path, _ = get_seekpath_kpoint_path(calc_doc,
                                                         spacing=spacing,
                                                         symmetry_tol=calc_doc.get('symmetry_tol'),
                                                         debug=debug)
        LOG.debug('New lattice: {}'.format(Crystal(calc_doc)))
        return prim_doc, kpt_path

    def run_convergence_tests(self, calc_doc):
        """ Run kpoint and cutoff_energy convergence tests based on
        options passed to ComputeTask.

        Parameters:
            calc_doc (dict): the structure to converge.

        Returns:
            bool: True unless every single calculation failed.

        """
        LOG.info('Performing convergence tests...')
        from matador.utils.cell_utils import get_best_mp_offset_for_cell
        successes = []
        cached_cutoff = calc_doc['cut_off_energy']
        if self.conv_cutoff_bool:
            # run series of singlepoints for various cutoffs
            LOG.info('Running cutoff convergence...')
            for cutoff in self.conv_cutoff:
                LOG.info('{} eV... '.format(cutoff))
                calc_doc.update({'cut_off_energy': cutoff})
                self.paths['completed_dir'] = 'completed_cutoff'
                seed = self.seed + '_' + str(cutoff) + 'eV'
                success = self.run_castep_singleshot(calc_doc, seed, keep=False)
                successes.append(success)
        if self.conv_kpt_bool:
            # run series of singlepoints for various cutoffs
            LOG.info('Running kpt convergence tests...')
            calc_doc['cut_off_energy'] = cached_cutoff
            for kpt in self.conv_kpt:
                LOG.info('{} 1/A... '.format(kpt))
                calc_doc.update({'kpoints_mp_spacing': kpt})
                calc_doc['kpoints_mp_offset'] = get_best_mp_offset_for_cell(calc_doc)
                LOG.debug('Using offset {}'.format(calc_doc['kpoints_mp_offset']))
                self.paths['completed_dir'] = 'completed_kpts'
                seed = self.seed + '_' + str(kpt) + 'A'
                success = self.run_castep_singleshot(calc_doc, seed, keep=False)
                successes.append(success)
        return any(successes)

    def parse_executable(self, seed):
        """ Turn executable string into list with arguments to be executed.

        Example:
            With `self.executable='castep17'` and `seed='test'`, `['castep17', 'test']`
            will be returned.

        Example:
            With `self.executable='pw6.x -i $seed.in > $seed.out'` and `seed='test'`,
            `['pw6.x', '-i', 'test.in', '>' 'test.out']` will be returned.

        Parameters:
            seed (str): filename to replace $seed with in command.

        Returns:
            :obj:`list` of :obj:`str`: list called by subprocess.POpen.

        """
        if isinstance(self.executable, str):
            executable = self.executable.split()
        command = []
        found_seed = False
        for item in executable:
            if '$seed' in item:
                item = item.replace('$seed', seed)
                found_seed = True
            command.append(item)
        if not found_seed:
            command.append(seed)

        if self.redirect is not None:
            self._redirect_filename = self.redirect.replace('$seed', seed)
        else:
            self._redirect_filename = None

        LOG.debug('Executable string parsed as {}'.format(command))

        return command

    def test_exec(self):
        """ Test if <executable> --version returns a valid string.

        Raises:
            CriticalError: if executable not found.

        """
        LOG.info('Testing executable "{executable}".'.format(executable=self.executable))

        try:
            proc = self.run_command('--version')

        except FileNotFoundError:
            LOG.critical('Unable to call mpirun/aprun/srun, currently selected: {}'.format(self.mpi_library))
            message = 'Please check initialistion of ComputeTask object/CLI args.'
            LOG.debug('Raising CriticalError with message: {message}'.format(message=message))
            raise CriticalError(message)

        out, errs = proc.communicate()
        out = out.decode('utf-8')
        errs = errs.decode('utf-8')
        if 'CASTEP version' not in out:
            # this is an OpenMPI error that occurs when hyperthreading is enabled
            # best way to handle is to half the number of procs available
            if 'not enough slots' in errs:
                err_string = ('MPI library tried to use too many cores and failed, '
                              'rescaling core count and re-running with {} cores...'.format(int(self.ncores/2)))
                LOG.warning(err_string)
                LOG.warning('stdout: {stdout}'.format(stdout=out))
                LOG.warning('sterr: {stderr}'.format(stderr=errs))
                if self.ncores >= 2:
                    self.ncores = int(self.ncores/2)
                else:
                    raise CriticalError(err_string)
                self.test_exec()

            else:
                err_string = 'Executable `{}` failed testing: does it support --version flag?'.format(self.executable)
                LOG.critical(err_string)
                LOG.critical('stdout: {stdout}'.format(stdout=out))
                LOG.critical('sterr: {stderr}'.format(stderr=errs))
                LOG.debug('Raising CriticalError.')
                raise CriticalError(err_string)

        if errs:
            LOG.info('Executable {} passed test, but stderr contains the following:'.format(self.executable))
            LOG.info(errs)

        try:
            version_string = out.split('\n')
            for line in version_string:
                if 'CASTEP version' in line:
                    version = line.strip().split(' ')[-1]
                    version = float(version)

            LOG.info('CASTEP version: {}'.format(version))
            if version < 18:
                LOG.warning('Using CASTEP version {}: '
                            'some features may not work, consider updating to at least CASTEP 18'.format(version))
        except Exception:
            LOG.warning('Unable to detect CASTEP version: '
                        'some features may not work unless you are using at least version 18.')

    @property
    def mpi_library(self):
        """ Property to store/compute desired MPI library. """
        if self._mpi_library is None:
            self._mpi_library = self.set_mpi_library()
        return self._mpi_library

    def set_mpi_library(self):
        """ Combines command-line MPI arguments into string and calls
        MPI library detection is no args are present.
        """
        guessed_version = self.detect_mpi()
        LOG.info('Detected {} MPI.'.format(guessed_version))
        if sum([self.archer, self.intel, self.slurm]) > 1:
            message = 'Conflicting command-line arguments for MPI library have been supplied, exiting.'
            LOG.critical(message)
            raise CriticalError(message)

        if self.archer:
            if guessed_version != 'archer':
                message = ('Detected {} MPI, but user asked to use aprun... please check your environment.'
                           .format(self._mpi_library))
                LOG.critical(message)
                raise CriticalError(message)
            return 'archer'

        if self.intel:
            if guessed_version != 'intel':
                message = ('Detected {} MPI, but user asked to use Intel MPI... please check your environment.'
                           .format(self._mpi_library))
                LOG.critical(message)
                raise CriticalError(message)
            return 'intel'

        if self.slurm:
            if guessed_version != 'slurm':
                message = 'Detected {} MPI, but user asked to use srun... continuing with srun.'.format(self._mpi_library)
                LOG.warning(message)
            return 'slurm'

        return guessed_version

    @staticmethod
    def detect_mpi():
        """ Test which mpi library is being used when `mpirun`.

        Returns:
            mpi_library (str): 'intel', 'archer', or 'default'.

        """
        # check first for existence of mpirun command, then aprun if that fails
        try:
            try:
                mpi_version_string = str(sp.check_output('mpirun --version', shell=True))
            except sp.CalledProcessError:
                LOG.info('Failed to find mpirun, checking aprun...')
                mpi_version_string = str(sp.check_output('aprun --version', shell=True))
        except Exception as exc:
            msg = 'Failed to find mpirun or aprun.'
            LOG.critical(msg)
            LOG.debug('Error message: {exc}'.format(exc=exc))
            raise CriticalError(msg)
        if 'Intel' in mpi_version_string:
            mpi_version = 'intel'
        elif 'aprun' in mpi_version_string:
            mpi_version = 'archer'
        elif 'Open MPI' in mpi_version_string:
            mpi_version = 'default'
        else:
            LOG.debug('Could not detect MPI library so using default (OpenMPI), version string was: {response}'
                      .format(response=mpi_version_string))
            mpi_version = 'default'

        LOG.info('Using {version} MPI library.'.format(version=mpi_version))
        return mpi_version

    def do_memcheck(self, calc_doc, seed):
        """ Perform a CASTEP dryrun to estimate memory usage.

        Parameters:
            calc_doc (dict): dictionary of structure and CASTEP parameters
            seed (str): filename for structure

        Returns:
            bool: True if the memory estimate is <90% of node RAM or
                `self.maxmem`, if set

        """
        LOG.info('Performing memory check for {seed}'.format(seed=seed))
        memcheck_seed = seed + '_memcheck'

        memcheck_doc = deepcopy(calc_doc)

        if memcheck_doc['task'] in MATADOR_CUSTOM_TASKS:
            memcheck_doc['task'] = 'singlepoint'

        doc2param(memcheck_doc, memcheck_seed, hash_dupe=False)
        doc2cell(memcheck_doc, memcheck_seed, hash_dupe=False)

        LOG.debug('Running CASTEP dryrun.')
        self.executable += ' --dryrun'
        process = self.run_command(memcheck_seed)
        process.communicate()
        self.executable = self.executable.replace(' --dryrun', '')

        results, success = castep2dict(memcheck_seed + '.castep', db=False, verbosity=self.verbosity)

        for _file in glob.glob(memcheck_seed + '*'):
            if _file.endswith('.res'):
                continue
            else:
                os.remove(_file)

        if not success or 'estimated_mem_MB' not in results:
            msg = 'CASTEP dryrun failed with output {results}'.format(results=results)
            raise MaxMemoryEstimateExceeded(msg)

        estimate = results['estimated_mem_per_process_MB']
        if self.ncores is not None:
            estimate *= self.ncores
        if self.nnodes is not None:
            estimate *= self.nnodes

        return estimate

    def run_command(self, seed):
        """ Calls executable on seed with desired number of cores.

        Parameters:
            seed (str): seedname to pass append to CASTEP command,
                e.g. <seed> or --version.

        Returns:
            subprocess.Popen: process to run.

        """
        command = self.parse_executable(seed)
        if self.nnodes is None or self.nnodes == 1:
            if self.node is None:
                if self.ncores == 1 and self.node is None:
                    command = ['nice', '-n', '15'] + command
                elif self.mpi_library == 'archer':
                    command = ['aprun', '-n', str(self.ncores)] + command
                elif self.mpi_library == 'slurm':
                    command = ['srun'] + command
                elif self.mpi_library in ['intel', 'default']:
                    command = ['mpirun', '-n', str(self.ncores)] + command
                else:
                    command = ['mpirun', '-n', str(self.ncores)] + command
            elif self.node is not None:
                cwd = os.getcwd()
                command = ['ssh', '{}'.format(self.node), 'cd', '{};'.format(cwd), 'mpirun', '-n',
                           str(self.ncores)] + command
        else:
            if self.mpi_library == 'archer':
                command = [
                    'aprun', '-n',
                    str(self.ncores * self.nnodes), '-N',
                    str(self.ncores), '-S', '12', '-d', '1'
                ] + command
            elif self.mpi_library == 'slurm':
                command = ['srun'] + command
            elif self.mpi_library == 'intel':
                command = ['mpirun', '-n', str(self.ncores * self.nnodes), '-ppn', str(self.ncores)] + command
            else:
                command = ['mpirun', '-n', str(self.ncores * self.nnodes), '-npernode', str(self.ncores)] + command

        if self.ncores > 1 and 'run' not in ' '.join(command):
            raise RuntimeError('Issue running command {} in parallel, this is probably a matador bug!'.format(command))

        stdout = sp.PIPE
        stderr = sp.PIPE

        if self._redirect_filename is not None:
            LOG.info('Redirecting output to {redirect}'.format(redirect=self._redirect_filename))
            redirect_file = open(self._redirect_filename, 'w')
            stdout = redirect_file

        LOG.info('Running {}'.format(command))
        process = sp.Popen(command, shell=False, stdout=stdout, stderr=stderr)
        try:
            redirect_file.close()
        except Exception:
            pass

        return process

    def _handle_process(self, process, check_walltime=False, expected_fname=None):
        """ Poll, wait, communicate with running process, with optional
        checks on walltime and file creation.

        Parameters:
            process (subprocess.Popen): process object to handle.

        Keyword arguments:
            check_walltime (bool): timeout the process 2 x :attr:`polltime` before
                the allotted walltime, and clean up any files.
            expected_fname (str): the filename to check for, assuming filename

        Raises:
            CalculationError: if calculation reports its own failure through
                stderr or non-zero return code, or if expected_fname check fails.
            WalltimeError: if the calculation has not completed before the end
                of the allotted walltime.

        Returns:
            (str, str): any output to stdout/stderr.

        """

        proc_clock = time.time()

        def _check_file_has_been_written(fname, proc_start):
            """ Check if the file at `fname` has been written to by the current process.
            Assumes a 1 second lee-way in the IO, in case the process wrote to the file
            just before starting the timer.

            Parameters:
                fname (str): the filename to check.
                proc_start (int): Unix time in seconds at which the process started.

            Raises:
                CalculationError: if the file does not exist, or if it is too old.

            """

            if not os.path.isfile(fname):
                msg = ('File {} was not created after {} s, please check your executable: {}.'
                       .format(fname, self.polltime, self.executable))
                LOG.critical(msg)
                raise CalculationError(msg)
            if os.path.getmtime(fname) - proc_clock < -1:
                msg = ('File {} present, but too old to be made by this process. '
                       'Please check your executable: {}.'
                       .format(fname, self.executable))
                LOG.critical(msg)
                raise CalculationError(msg)

        if expected_fname is not None:
            LOG.info("Watching for file write to {} after {} s".format(expected_fname, self.polltime))
            while process.poll() is None:
                proc_elapsed = time.time() - proc_clock
                if proc_elapsed > self.polltime:
                    _check_file_has_been_written(expected_fname, proc_clock)
                    break
                time.sleep(1)

        timeout = None
        if check_walltime and self.max_walltime is not None:
            LOG.info("Will not let process walltime exceed {} seconds".format(self.max_walltime))
            proc_elapsed = time.time() - proc_clock
            timeout = max(self.max_walltime - 2 * self.polltime - proc_elapsed, 1)

        try:
            out, errs = process.communicate(timeout=timeout)

            out = out.decode("utf-8")
            errs = errs.decode("utf-8")
            if process.returncode != 0 or errs:
                # as there are several reasons why the process can return != 0, handle
                # them in turn for each case, rather than raising an error here
                LOG.warning('Process returned error code {}'.format(process.returncode))
                LOG.warning('\nstdout: {}'.format(out))
                LOG.warning('\nstderr: {}'.format(errs))

        except sp.TimeoutExpired:
            LOG.error("Process reached maximum walltime, cleaning up...")
            self._times_up(process)
            raise WalltimeError("Cleaned up process after reaching maximum walltime")

        except Exception as err:
            LOG.error('Unexpected Exception {} caught: terminating job for {}.'.format(type(err).__name__, self.seed))
            raise err

        finally:
            LOG.error("Explicitly terminating process.")
            process.terminate()

        return out, errs

    def _catch_castep_errors(self, process):
        """ Look for CASTEP error files and fallover appropriately. If
        the magic string 'Work-around was succesful' is found in error,
        no errors will be reported and the file will be deleted,
        provided no other errors exist.
            results_dict (dict/str): dictionary containing data to u

        Returns:
            (bool, str, function): True if error files were found, otherwise False,
                followed by the error messages. If the error can be remedied,
                return the function to attempt to fix it.

        """
        err_file = '{}*err'.format(self.seed)
        error_str = ''
        errors_present = False
        remedy = None

        if process.returncode != 0:
            msg = 'CASTEP returned non-zero error code {}.'.format(process.returncode)
            error_str += msg + '\n'
            errors_present = True

        for globbed in glob.glob(err_file):
            if globbed.endswith('opt_err'):
                continue
            if os.path.isfile(globbed):
                with open(globbed, 'r') as f:
                    flines = f.readlines()

                for line in flines:
                    if 'Work-around was successful, continuing with calculation.' in line:
                        LOG.info('Found LAPACK issue that was circumvented, removing error file.')
                        os.remove(globbed)
                        break
                    elif 'ERROR in cell constraints: attempt to fix' in line:
                        remedy = self._remedy_castep_symmetry_error
                        break
                    elif 'Symmetry matrix not integer' in line:
                        remedy = self._remedy_symmetry_matrix_not_integer
                        break
                else:
                    error_str += ' '.join(flines)
                    error_str += '\n'
                    errors_present = True

        return errors_present, error_str, remedy

    @staticmethod
    def _remedy_castep_symmetry_error(opti_dict):
        """ Remedy the common cell constrains error by
        disabling symmetry.

        Parameters:
            opti_dict (dict): the dictionary of parameters to change.

        """
        LOG.info('Trying to remedy CASTEP symmetry error...')
        if 'symmetry_generate' in opti_dict:
            del opti_dict['symmetry_generate']
        if 'symmetry_tol' in opti_dict:
            del opti_dict['symmetry_tol']
        if 'snap_to_symmetry' in opti_dict:
            del opti_dict['snap_to_symmetry']

    def _remedy_symmetry_matrix_not_integer(self, opti_dict):
        LOG.info('Trying to remedy symmetry matrix not integer bug: symmetry tol now = {}'
                 .format(opti_dict['symmetry_tol']))
        if 'symmetry_tol' in opti_dict:
            opti_dict['symmetry_tol'] /= 10

    def mv_to_bad(self, seed):
        """ Move all files associated with "seed" to bad_castep, from both
        the compute directory (if it exists) and the root dir..

        Parameters:
            seed (str): filename of structure.

        """
        try:
            bad_dir = self.root_folder + '/' + self.paths["failed_dir"]
            if not os.path.exists(bad_dir):
                os.makedirs(bad_dir, exist_ok=True)
            seed_files = glob.glob(seed + '.*') + glob.glob(seed + '-out.cell*')
            if seed_files:
                LOG.info('Moving files to bad_castep: {bad}.'.format(bad=bad_dir))
                LOG.debug('Files to move: {seed}'.format(seed=seed_files))
                for _file in seed_files:
                    try:
                        shutil.copy2(_file, bad_dir)
                        os.remove(_file)
                    except Exception as exc:
                        LOG.warning('Error moving files to bad: {error}'.format(error=exc))
            # check root folder for any matching files and remove them
            fname = '{}/{}'.format(self.root_folder, seed)
            for ext in ['.res', '.res.lock', '.castep', '-out.cell']:
                if os.path.isfile('{}{}'.format(fname, ext)):
                    os.remove('{}{}'.format(fname, ext))
        except Exception as exc:
            LOG.warning('Error moving files to bad: {error}'.format(error=exc))

    def mv_to_completed(self, seed, completed_dir='completed', keep=False, skip_existing=False):
        """ Move all associated files to completed, removing any
        remaining files in the root_folder and compute_dir.

        Parameters:
            seed (str): filename for structure.

        Keyword arguments:
            completed_dir (str): folder for completed jobs.
            keep (bool): whether to also move intermediate files.
            skip_existing (bool): if True, skip files that already exist,
                otherwise throw an error.

        """
        completed_dir = self.root_folder + '/' + completed_dir
        LOG.info('Moving files to completed: {completed}'.format(completed=completed_dir))

        if seed.endswith('.res'):
            seed = str(seed.replace('.res', ''))

        if not os.path.exists(completed_dir):
            os.makedirs(completed_dir, exist_ok=True)

        for _file in glob.glob(seed + '*_bak') + glob.glob(seed + '*.lock'):
            os.remove(_file)
        if keep:
            seed_files = glob.glob(seed + '.*') + glob.glob(seed + '-out.cell*')
            if seed_files:
                LOG.debug('Files to move: {files}.'.format(files=seed_files))
                for _file in seed_files:
                    if skip_existing and os.path.isfile(completed_dir + '/' + _file):
                        LOG.warning('File already found {}...'.format(_file))
                    else:
                        shutil.move(_file, completed_dir)
        else:
            # move castep/param/res/out_cell files to completed
            file_exts = ['.castep']
            if self.kpts_1D:
                file_exts.append('.param')
            if not self.conv_kpt_bool and not self.conv_cutoff_bool:
                file_exts.append('.res')
            if os.path.isfile(seed + '-out.cell'):
                file_exts.append('-out.cell')
            if self.calc_doc.get('write_formatted_density'):
                file_exts.append('.den_fmt')
            for ext in file_exts:
                try:
                    shutil.move('{}{}'.format(seed, ext), completed_dir)
                except Exception as exc:
                    LOG.warning('Error moving files to completed: {error}'.format(error=exc))

        # delete whatever is left
        wildcard_fnames = glob.glob(seed + '.*')
        wildcard_fnames += glob.glob(seed + '-out.*')
        for fname in wildcard_fnames:
            os.remove(fname)

    def cp_to_input(self, seed, ext='res', glob_files=False):
        """ Copy initial cell and res to input folder.

        Parameters:
            seed (str): filename of structure.

        Keyword arguments:
            ext (str): file extension for structure.
            glob_files (bool): whether to glob all related seed files.

        """
        if not self._first_run:
            return
        input_dir = self.root_folder + '/input'
        LOG.debug('Copying file to input_dir: {input}'.format(input=input_dir))
        if not os.path.exists(input_dir):
            os.makedirs(input_dir, exist_ok=True)
        if glob_files:
            files = glob.glob('{}.*'.format(seed))
            LOG.debug('Files to copy: {files}'.format(files=files))
            for f in files:
                if f.endswith('.lock'):
                    continue
                if not os.path.isfile(input_dir + '/' + f):
                    shutil.copy2('{}'.format(f), input_dir)
        else:
            LOG.debug('File to copy: {file}'.format(file='{}.{}'.format(seed, ext)))
            if os.path.isfile('{}.{}'.format(seed, ext)):
                if not os.path.isfile('{}/{}.{}'.format(input_dir, seed, ext)):
                    shutil.copy2('{}.{}'.format(seed, ext), input_dir)

    def _setup_relaxation(self):
        """ Set up directories and files for relaxation.

        Raises:
            CalculationError: if structure has already exceeded
                geom_max_iter.
            CriticalError: if unable to split up calculation, normally
                indicating geom_max_iter is too small.

        """

        LOG.info('Preparing to relax {seed}'.format(seed=self.seed))
        self._setup_compute_dir(self.seed, self.compute_dir, custom_params=self.custom_params)

        # update res file with intermediate calculation if castep file is newer than res
        if os.path.isfile(self.seed + '.castep') and os.path.isfile(self.seed + '.res'):
            LOG.info('Trying to update res file with result from intermediate CASTEP file found in root_dir')
            if self.compute_dir is not None:
                shutil.copy2(self.seed + '.castep', self.compute_dir)
            castep_dict, success = castep2dict(self.seed + '.castep', db=False, verbosity=self.verbosity)
            if success:
                self.res_dict['geom_iter'] = castep_dict.get('geom_iter', 0)

            if os.path.getmtime(self.seed + '.res') < os.path.getmtime(self.seed + '.castep'):
                LOG.info('CASTEP file was updated more recently than res file, using intermediate structure...')
                self.res_dict.update(castep_dict)

        if self.compute_dir is not None:
            os.chdir(self.compute_dir)

        # copy initial res file to seed
        LOG.info('Writing fresh res file to start calculation from.')
        doc2res(self.res_dict, self.seed, info=False, hash_dupe=False, overwrite=True)
        self.cp_to_input(self.seed)

        # set up geom opt iteration options based on input/scraped parameters
        self._max_iter = self.calc_doc['geom_max_iter']
        if self.res_dict.get('geom_iter'):
            self._max_iter -= self.res_dict['geom_iter']
        else:
            self.res_dict['geom_iter'] = 0

        if self.res_dict['geom_iter'] > self.calc_doc['geom_max_iter']:
            msg = '{} iterations already performed on structure, exiting...'.format(self.res_dict['geom_iter'])
            LOG.critical(msg)
            raise CalculationError(msg)

        # number of steps in fine and rough calcs respectively
        fine_iter = self.fine_iter
        rough_iter = self.rough_iter
        # number of calculations to run with rough_iter steps
        num_rough_iter = self.rough
        if 'geom_method' in self.calc_doc:
            if self.calc_doc['geom_method'].lower() == 'tpsd' and rough_iter < 3:
                rough_iter = 3
        self._geom_max_iter_list = (num_rough_iter * [rough_iter])
        self._max_iter -= num_rough_iter * rough_iter
        if self.squeeze:
            self._squeeze_list = [True for val in self._geom_max_iter_list]
        else:
            self._squeeze_list = [False for val in self._geom_max_iter_list]

        num_fine_iter = ceil(int(self._max_iter) / fine_iter)
        if self._max_iter > 0:
            if self._max_iter < fine_iter:
                fine_iter = self._max_iter
                num_fine_iter = 1
            self._geom_max_iter_list.extend(num_fine_iter * [fine_iter])
            self._squeeze_list.extend(num_fine_iter * [False])

        LOG.info('Geometry optimisation iteration scheme set to {}'.format(self._geom_max_iter_list))
        if self.squeeze:
            LOG.info('Squeeze scheme set to {}'.format(self._geom_max_iter_list))

        if not self._geom_max_iter_list:
            msg = 'Could not divide up relaxation; consider increasing geom_max_iter'
            LOG.critical(msg)
            raise CriticalError(msg)

    def _update_input_files(self, seed, calc_doc, squeeze=None):
        """ Update the cell and param files for the next relaxation.

        Parameters:
            seed (str): the seedname to update.
            calc_doc (dict): the calculation dictionary to write to file.

        Keyword arguments:
            squeeze (float): external pressure to add this step

        """
        if seed is None:
            seed = self.seed

        this_calc_doc = deepcopy(calc_doc)

        # update cell
        if os.path.isfile(seed + '.cell'):
            os.remove(seed + '.cell')
        if squeeze is not None:
            if squeeze:
                LOG.info('Applying pressure of {} GPa to this calculation'.format(squeeze))
                this_calc_doc['external_pressure'] = [[squeeze, 0, 0], [0, squeeze, 0], [0, 0, squeeze]]
            else:
                LOG.info('Pressure reset to {}'.format(self._target_pressure))
                this_calc_doc['external_pressure'] = self._target_pressure

        if self.kpts_1D:
            n_kz = ceil(1 / (this_calc_doc['lattice_abc'][0][2] * self._target_spacing))
            if n_kz % 2 == 1:
                n_kz += 1
            this_calc_doc['kpoints_mp_grid'] = [1, 1, n_kz]
            if 'kpoints_mp_spacing' in calc_doc:
                del calc_doc['kpoints_mp_spacing']
        doc2cell(this_calc_doc, seed, hash_dupe=False, spin=self.spin)

        # update param
        if not self.custom_params:
            if os.path.isfile(seed + '.param'):
                os.remove(seed + '.param')
            doc2param(this_calc_doc, seed, hash_dupe=False, spin=self.spin)

        LOG.debug('Calculation dictionary: {}'
                  .format(matador.utils.print_utils.dumps(this_calc_doc,
                                                          indent=None)))

    def tidy_up(self, seed):
        """ Delete all created files before quitting.

        Parameters:
            seed (str): filename for structure.

        """
        files = glob.glob(seed + '.*')
        if self.compute_dir is not None:
            files += glob.glob(self.root_folder + '/' + seed + '.*')
        if files:
            LOG.info('Tidying up remaining files: {files}'.format(files=files))
            for f in files:
                # if we're working in a compute dir, then delete any remaining files in base dir
                # otherwise, only delete things that we dont care about, i.e. non-res/castep in
                # case they were not correctly moved to completed/bad_castep by the other routines
                if self.compute_dir is not None or (not (f.endswith('.res') or f.endswith('.castep'))):
                    os.remove(f)

    def _update_castep_output_files(self, seed, opti_dict=None):
        """ Copy new data to root CASTEP output files and update
         the results dict.

        Keyword arguments:
            opti_dict (dict): intermediate calculation results.

       """
        LOG.info('Updating .res and .castep files in root_dir with new results')
        if opti_dict is not None:
            if os.path.isfile(seed + '.res'):
                os.rename('{}.res'.format(seed), '{}.res_bak'.format(seed))
            try:
                doc2res(opti_dict, seed, hash_dupe=False)
            except CalculationError:
                doc2res(opti_dict, seed, hash_dupe=False, info=False)
            if os.path.isfile(seed + '.res_bak'):
                os.remove(seed + '.res_bak')
            self.res_dict.update(opti_dict)

        if self.compute_dir is not None:
            if os.path.isfile(seed + '.res'):
                shutil.copy2(seed + '.res', self.root_folder)
            if os.path.isfile(seed + '.castep'):
                shutil.copy2(seed + '.castep', self.root_folder)

    def _finalise_result(self, intermediate=False):
        """ Push to queue if necessary and return status.

        Keyword arguments:
            intermediate (bool): whether we want to run more calculations
                on the output of this, i.e. whether to move to completed
                or not.

        Returns:
            bool: True is relaxation was successful, False otherwise.

        """
        LOG.info('Finalising calculation...')
        try:
            success = self.res_dict.get('optimised', False)
        except AttributeError:
            success = False

        LOG.info('Was calculation successful? {success}'.format(success=success))
        if self.output_queue is not None:
            LOG.info('Pushing results to output queue')
            self.output_queue.put(self.res_dict)
        if success:
            if not intermediate:
                self.mv_to_completed(self.seed, completed_dir=self.paths['completed_dir'])
        else:
            self.mv_to_bad(self.seed)

        if success:
            self.final_result = self.res_dict

        if not intermediate:
            # clean up rest of files
            self.tidy_up(self.seed)

        return success

    def _times_up(self, process):
        """ If walltime has nearly expired, run this function
        to kill the process and unlock it for restarted calculations.

        Parameters:
            subprocess.Popen: running process to be killed.

        """
        if self.compute_dir is not None:
            LOG.info('Cleaning up compute_dir: {dir}'.format(dir=self.compute_dir))
            for f in glob.glob('{}.*'.format(self.seed)):
                shutil.copy2(f, self.root_folder)
                os.remove(f)
        LOG.info('Removing lock file so calculation can be continued.')
        if os.path.isfile('{}/{}{}'.format(self.root_folder, self.seed, '.res.lock')):
            os.remove('{}/{}{}'.format(self.root_folder, self.seed, '.res.lock'))

    @staticmethod
    def remove_compute_dir_if_finished(compute_dir):
        """ Delete the compute directory, provided it contains no
        calculation data.

        Parameters:
            compute_dir (str): path to compute directory.

        Returns:
            bool: True if folder was deleted as no res/castep files
                were found, otherwise False.

        """
        LOG.info('Checking if compute_dir still contains calculations...')
        if not os.path.isdir(compute_dir):
            return False

        files = glob.glob(compute_dir + '/*')
        LOG.debug('Found {files} in {dir}'.format(files=files, dir=compute_dir))

        for fname in files:
            if fname.endswith('.res') or fname.endswith('.castep'):
                LOG.debug('Not removing {dir} as it still contains calculation {fname}'.format(
                    dir=compute_dir, fname=fname))
                return False

        # remove files in directory, then delete directory
        LOG.debug('Deleting files {files} from {dir}'.format(files=files, dir=compute_dir))
        for fname in files:
            if os.path.isfile(fname):
                try:
                    os.remove(fname)
                except FileNotFoundError:
                    pass

        if os.path.isdir(compute_dir):
            LOG.debug('Deleting directory {dir}'.format(dir=compute_dir))
            try:
                os.rmdir(compute_dir)
            except OSError:
                LOG.debug('Unable to delete directory {} as it still contains files.'.format(compute_dir))

        if os.path.islink(compute_dir.split('/')[-1] + '_link'):
            os.remove(compute_dir.split('/')[-1] + '_link')

        return True

    def _setup_compute_dir(self, seed, compute_dir, custom_params=False, generic=False):
        """ Create the desired directory if it doens't exist,
        and try to link to it in the current folder.

        Parameters:
            seed (str): name of seed.
            compute_dir (str): name of directory to make.

        Keyword arguments:
            custom_params (bool): whether to try to copy custom
                param files into this directory.

        """
        if compute_dir is None:
            return

        LOG.info('Using compute_dir: {}'.format(compute_dir))
        if not os.path.isdir(compute_dir):
            try:
                os.makedirs(compute_dir)
            except PermissionError as exc:
                raise CriticalError('Invalid compute dir requested: {} {}'.format(compute_dir, exc))
        # if compute_dir isn't simply inside this folder, make a symlink that is
        if compute_dir.startswith('/'):
            link_name = compute_dir.split('/')[-1] + '_link'
            if not os.path.exists(link_name):
                os.symlink(compute_dir, link_name)

        if generic:
            # if generic, copy all seed files to compute dir
            for f in glob.glob(seed + '*'):
                shutil.copy2(f, compute_dir)

        # copy pspots and any intermediate calcs to compute_dir
        LOG.info('Copying pspots into compute_dir')
        if self.cell_dict is not None:
            for pspot in self.cell_dict.get('species_pot', {}).values():
                if os.path.isfile(pspot):
                    try:
                        shutil.copy2(pspot, compute_dir)
                    except shutil.SameFileError:
                        pass

        if custom_params:
            shutil.copy2(seed + '.param', compute_dir)

        # if a checkfile exists, copy it to the new dir
        # so that it can be restarted from
        if os.path.isfile(seed + '.check'):
            shutil.copy2(seed + '.check', compute_dir)

    def scf(self, *args, **kwargs):
        """ Alias for backwards-compatibility. """
        import warnings
        warnings.warn(
            '`scf` method name is deprecated, please use run_castep_singleshot to avoid future issues.',
            DeprecationWarning
        )
        return self.run_castep_singleshot(*args, **kwargs)

    def relax(self, *args, **kwargs):
        """ Alias for backwards-compatibility. """
        import warnings
        warnings.warn(
            '`relax` method name is deprecated, please use run_castep_relaxation to avoid future issues.',
            DeprecationWarning
        )
        return self.run_castep_relaxation(*args, **kwargs)
