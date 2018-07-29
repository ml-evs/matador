# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the FullRelaxer class for handling
calculations on a single structure, and some useful associated
errors.

"""


import os
import sys
import shutil
import subprocess as sp
import glob
import time
import logging
from copy import deepcopy
from math import ceil
from psutil import virtual_memory

from matador.scrapers.castep_scrapers import cell2dict
from matador.scrapers.castep_scrapers import res2dict, castep2dict
from matador.export import doc2cell, doc2param, doc2res


class FullRelaxer:
    """ The main use of this class is to call an executable on a given
    structure. The various parameters are passed to this class by the
    common entrypoints, run3 and ilustrado. It is unlikely that you
    will want to use this class directly.

    Note:
        By default, calculations are run inside a folder with the same
        name as the host (e.g. node12, or whatever). This decreases the
        load on parallel file systems such as Lustre.

    """

    def __init__(self, res, ncores, nnodes, node, **kwargs):
        """ Make the files to run the calculation and call the desired program.

        Parameters:
            res (str/dict): filename or input structure dict
            param_dict (dict): dictionary of CASTEP parameters
            cell_dict (dict): dictionary of CASTEP cell options
            ncores (int): number of cores *per node* for mpirun call
            nnodes (int): number of nodes for mpirun call (if None, use 1)
            node (str): node name to run on (if None, run on localhost)

        Keyword arguments:
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
            spin (bool): break spin symmetry in first calculation
                (DEFAULT: False)
            conv_cutoffs (:obj:`list` of :obj:`float`): list of cutoffs
                to use for SCF convergence test
            conv_kpts (:obj:`list` of :obj:`float`): list of kpt spacings
                to use for SCF convergence test
            kpts_1D (bool): treat z-direction as special and create
                kpt_grid [1 1 n_kz] (DEFAULT: False)
            archer (bool): force use of aprun over mpirun (DEFAULT: False)
            slurm (bool): force use of srun over mpirun (DEFAULT: False)
            intel (bool): force use of Intel mpirun-style calls (DEFAULT: False)
            profile (bool): use cProfile to profile runtime (DEFAULT: False)
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
            SystemExit: if a fatal error occurs, failed run will be moved to bad_castep and
                no further calculations will be attempted.
            RuntimeError: if a structure-level error occurs, causing the seed files to be moved
                to bad_castep.

        """
        # set defaults and update class with desired values
        prop_defaults = {'paths': None, 'param_dict': None, 'cell_dict': None, 'mode': 'castep', 'executable': 'castep',
                         'memcheck': False, 'rough': 4, 'rough_iter': 2, 'fine_iter': 20, 'spin': False,
                         'output_queue': None, 'redirect': None, 'reopt': False, 'compute_dir': None, 'noise': False,
                         'custom_params': False, 'archer': False, 'maxmem': None, 'killcheck': True, 'kpts_1D': False,
                         'conv_cutoff': False, 'conv_kpt': False, 'profile': False, 'slurm': False, 'intel': False,
                         'exec_test': True, 'timings': (None, None), 'start': True, 'verbosity': 1, 'polltime': 30}
        self.paths = None
        self.process = None
        self.output_queue = None

        self.geom_max_iter_list = None
        self.max_iter = None
        self.num_rough_iter = None
        self.__dict__.update(prop_defaults)
        self.__dict__.update(kwargs)

        if self.profile:
            import cProfile
            import pstats
            from sys import version_info
            from matador import __version__
            profile = cProfile.Profile()
            profile.enable()

        # scrape and save seed name
        self.res = res
        if isinstance(self.res, str):
            self.seed = self.res
            if '/' in self.seed:
                shutil.copy(self.seed, self.seed.split('/')[-1])
        else:
            assert isinstance(self.res['source'], list)
            assert len(self.res['source']) == 1
            self.seed = self.res['source'][0]

        self.seed = self.seed.split('/')[-1].replace('.res', '')

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

        logging.basicConfig(level=logging.DEBUG)
        logging.getLogger().handlers = []

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(loglevel)
        stdout_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)8s: %(message)s'))
        logging.getLogger().addHandler(stdout_handler)

        logname = os.path.abspath('logs/{}.log'.format(self.seed))
        file_handler = logging.FileHandler(logname, mode='a')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)8s: %(message)s'))
        logging.getLogger().addHandler(file_handler)

        logging.info('Initialising FullRelaxer object for {seed}'.format(seed=self.seed))

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

        logging.debug('Starting at {start}, max walltime allowed is {walltime}'.format(
            start=self.start_time, walltime=self.max_walltime))

        self._redirect_filename = None

        if self.paths is None:
            self.paths = {}
            self.paths['completed_dir'] = 'completed'
        else:
            assert 'completed_dir' in self.paths

        # actually run the calculations, catching only RuntimeErrors
        try:
            # run through CASTEP specific features
            if self.mode == 'castep':
                self.success = self.run_castep(res)
            # otherwise run generic script
            else:
                self.success = self.run_generic(res)

        except RuntimeError as exc:
            logging.error('Process raised RuntimeError with message {error}.'.format(error=exc))
            self.success = False
            self._finalise_result()

        if self.profile:
            profile.disable()
            fname = 'relaxer-{}-{}-{}.{}.{}'.format(__version__, os.hostname()[1], version_info.major,
                                                    version_info.minor, version_info.micro)
            profile.dump_stats(fname + '.prof')
            with open(fname + '.pstats', 'w') as fp:
                stats = pstats.Stats(profile, stream=fp).sort_stats('cumulative')
                stats.print_stats()

        if self.success:
            logging.info('FullRelaxer finished successfully for {seed}'.format(seed=self.seed))
        else:
            logging.info('FullRelaxer failed cleanly for {seed}'.format(seed=self.seed))

    def run_castep(self, res):
        """ Set up and run CASTEP calculation.

        Parameters:
            res (str/dict): either the filename or a dict
                containing the structure.

        Raises:
            WalltimeError: if max_walltime is exceeded.
            SystemExit: if no further calculations should be performed
                on this thread.
            RuntimeError: if this structure errored in some way, but
                others will hopefully be okay.

        Returns:
            bool: True if calculations were successful, False otherwise.
                In the case of convergence tests, this is always True
                unless every calculation fails.

        """
        logging.info('Calling CASTEP on {seed}'.format(seed=self.seed))

        success = False

        if self.exec_test:
            self.test_exec()

        # pre-processing
        if self.kpts_1D:
            logging.debug('1D kpoint grid requested.')
            if 'kpoints_mp_spacing' not in self.cell_dict:
                raise SystemExit('kpoints_mp_spacing not found, but kpts_1D requested...')
            self._target_spacing = deepcopy(self.cell_dict['kpoints_mp_spacing'])

        # read in initial structure and skip if failed
        if isinstance(res, str):
            self.res_dict, success = res2dict(res, db=False)
            if not success:
                msg = 'Unable to parse initial res file, error: {res_dict}'.format(res_dict=self.res_dict)
                logging.error(msg)
                raise RuntimeError(msg)
        elif isinstance(res, dict):
            self.res_dict = res

        if self.noise:
            from matador.utils.cell_utils import add_noise
            logging.info('Adding noise to the positions in the cell')
            self.res_dict = add_noise(self.res_dict, amplitude=0.1)

        calc_doc = deepcopy(self.res_dict)

        # update global doc with cell and param dicts for folder
        calc_doc.update(self.cell_dict)
        calc_doc.update(self.param_dict)

        # check for pseudos
        if 'library' not in calc_doc['species_pot']:
            for elem in self.res_dict['stoichiometry']:
                if ('|' not in calc_doc['species_pot'][elem[0]] and
                        not os.path.isfile(os.path.expanduser(calc_doc['species_pot'][elem[0]]))):
                    msg = 'Unable to find pseudopotential file/string: {}'.format(calc_doc['species_pot'][elem[0]])
                    logging.critical(msg)
                    raise SystemExit(msg)

        # this is now a dict containing the exact calculation we are going to run
        self.calc_doc = calc_doc

        # do memcheck, if desired, and only continue if enough memory is free
        if self.memcheck:
            self.enough_memory = self.do_memcheck(calc_doc, self.seed)
        else:
            self.enough_memory = True

        if not self.enough_memory:
            return False

        # run convergence tests
        if any([self.conv_cutoff_bool, self.conv_kpt_bool]):
            success = self.run_convergence_tests(calc_doc)

        # perform relaxation
        elif calc_doc['task'].upper() in ['GEOMETRYOPTIMISATION', 'GEOMETRYOPTIMIZATION']:

            # begin relaxation
            if self.start:
                try:
                    success = self.relax()
                except Exception as err:
                    logging.debug('Cleaning up after catching error from relax.')
                    if self.compute_dir is not None:
                        # always cd back to root folder
                        os.chdir(self.root_folder)
                        self.remove_compute_dir_if_finished(self.compute_dir)
                    raise err

                if self.compute_dir is not None:
                    logging.debug('Cleaning up after relaxation.')
                    # always cd back to root folder
                    os.chdir(self.root_folder)
                    self.remove_compute_dir_if_finished(self.compute_dir)

        elif calc_doc['task'].upper() in ['PHONON', 'THERMODYNAMICS']:
            from matador.workflows.castep import castep_full_phonon
            success = castep_full_phonon(self, calc_doc, self.seed)

        elif calc_doc['task'].upper() in ['SPECTRAL']:
            from matador.workflows.castep import castep_full_spectral
            success = castep_full_spectral(self, calc_doc, self.seed)

        # run in SCF mode, i.e. just call CASTEP on the seeds
        else:
            success = self.scf(calc_doc, self.seed, keep=True)

        return success

    def run_generic(self, seed, intermediate=False):
        """ Run a generic mpi program on the given seed. Files from
        completed runs are moved to "completed" (unless intermediate
        is True) and failed runs to "bad_castep".

        Parameters:
            seed (str): filename of structure

        Keyword arguments:
            intermediate (bool): whether we want to run more calculations
                on the output of this, i.e. whether to move to completed
                or not.

        Returns:
            bool: True if calculations progressed without error.

        """
        logging.info('Calling "generic" MPI program on {seed}'.format(seed=self.seed))
        try:
            self.seed = seed
            if '.' in self.seed:
                self.seed = self.seed.split('.')[-2]
                self.input_ext = self.seed.split('.')[-1]
            else:
                self.input_ext = ''
            assert isinstance(self.seed, str)
            self.cp_to_input(seed, ext=self.input_ext, glob_files=True)
            self.process = self.run_command(seed)
            self.process.communicate()
            if self.process.returncode != 0:
                raise RuntimeError('Process returned code {}'.format(self.process.returncode))

            if not intermediate:
                logging.info('Writing results of generic call to res file and tidying up.')
                self.mv_to_completed(seed, keep=True, completed_dir=self.paths['completed_dir'])

            return True

        except Exception as err:
            logging.error('Caught error inside run_generic: {error}.'.format(error=err))
            self.mv_to_bad(seed)
            raise err

    def relax(self):
        """ Set up a structural relaxation that is restarted intermittently
        in order to re-mesh the kpoint grid. Completed calculations are moved
        to the "completed" folder, and failures to "bad_castep".

        Returns:
            bool: True iff structure was optimised, False otherwise.

        Raises:
            RuntimeError: if structure-level error occured.
            SystemExit: if fatal global error occured.
            WalltimeError: if walltime was reached, and jobs need to stop.

        """
        logging.info('Attempting to relax {}'.format(self.seed))

        self._setup_relaxation()
        seed = self.seed
        rerun = False
        try:
            # iterate over geom iter blocks
            for ind, num_iter in enumerate(self.geom_max_iter_list):

                # preprocess last step that did not finish geom opt
                if self.reopt and rerun:
                    # if we're now reoptimising after a success in last step
                    # use the fine iter value
                    num_iter = self.fine_iter
                    logging.info('Last step was successful, performing one last relaxation...')

                # update the geom_max_iter to use with either the number in iter_list, or the overriden value
                self.calc_doc['geom_max_iter'] = num_iter

                # delete any existing files and write new ones
                self._update_input_files()

                # run CASTEP
                self.process = self.run_command(seed)

                # if specified max_walltime (or found through SLURM), then monitor job
                if self.max_walltime is not None:
                    if self.start_time is None:
                        msg = 'Somehow initial start time was not found'
                        logging.critical(msg)
                        raise SystemExit(msg)

                    logging.info('Polling process every {} s'.format(self.polltime))

                    while self.process.poll() is None:
                        elapsed = time.time() - self.start_time
                        logging.debug('Elapsed time: {:.0f} s / {:.0f} s'
                                      .format(elapsed, self.max_walltime))

                        # leave 1 minute to clean up
                        if elapsed > abs(self.max_walltime - 3*self.polltime):
                            msg = 'About to run out of time on seed {}, killing early...'.format(self.seed)
                            logging.info(msg)
                            raise WalltimeError(msg)

                        time.sleep(self.polltime)

                self.process.communicate()
                logging.debug('Process returned {}'.format(self.process.returncode))
                if self.process.returncode != 0:
                    msg = 'CASTEP returned non-zero error code.'
                    logging.warning(msg)
                    raise RuntimeError(msg)

                # scrape new structure from castep file
                if not os.path.isfile(seed + '.castep'):
                    msg = ('CASTEP file was not created, please check your executable: {}.'
                           .format(self.executable))
                    logging.critical(msg)
                    raise SystemExit(msg)

                opti_dict, success = castep2dict(seed + '.castep', db=False, verbosity=self.verbosity)
                if not success and isinstance(opti_dict, str):
                    msg = 'Failed to parse CASTEP file...'
                    logging.warning(msg)
                    raise RuntimeError(msg)

                logging.debug('Intermediate calculation completed successfully: num_iter = {} '.format(num_iter))

                # scrub keys that need to be rescraped
                keys_to_remove = ['kpoints_mp_spacing', 'kpoints_mp_grid', 'species_pot', 'sedc_apply', 'sedc_scheme']
                for key in keys_to_remove:
                    if key in opti_dict:
                        del opti_dict[key]

                # now post-process the last step
                # did we try to rerun, and are now no longer optimised?
                # then reset, and go again...
                if self.reopt and rerun and not opti_dict['optimised']:
                    rerun = False

                # are we optimised, but haven't yet rerun?
                # then set prepare to do one more full relaxation
                if self.reopt and not rerun and opti_dict['optimised']:
                    rerun = True
                    self._update_output_files(opti_dict)

                # or did the relaxation complete successfuly, including rerun?
                elif (not self.reopt or rerun) and opti_dict['optimised']:
                    logging.info('Successfully relaxed {}'.format(seed))
                    self._update_output_files(opti_dict)
                    break

                # reached maximum number of steps
                elif ind == len(self.geom_max_iter_list) - 1:
                    logging.info('Failed to optimise {}'.format(seed))
                    break

                errors_present, errors = self._catch_castep_errors()
                if errors_present:
                    msg = 'Failed to optimise {} as CASTEP crashed with error:'.format(seed)
                    msg += errors
                    logging.warning(msg)
                    self._update_output_files(opti_dict)
                    break

                # if we weren't successful, then preprocess the next step

                # set atomic_init_spins with value from CASTEP file, if it exists
                if 'mulliken_spins' in opti_dict:
                    self.calc_doc['atomic_init_spins'] = opti_dict['mulliken_spins']

                # if writing out cell, use it for higher precision lattice_cart
                if self.calc_doc.get('write_cell_structure'):
                    cell_dict, success = cell2dict(seed + '-out.cell',
                                                   verbosity=self.verbosity,
                                                   db=False, outcell=True)
                    if success:
                        opti_dict['lattice_cart'] = list(cell_dict['lattice_cart'])

                logging.debug('N = {iters:03d} | |F| = {d[max_force_on_atom]:5.5f} eV/A | '
                              'S = {d[pressure]:5.5f} GPa | H = {d[enthalpy_per_atom]:5.5f} eV/atom'
                              .format(d=opti_dict, iters=sum(self.geom_max_iter_list[:ind+1])))

                self.calc_doc.update(opti_dict)

            return self._finalise_result()

        # catch WalltimeErrors and reset the job folder ready for continuation
        except WalltimeError as err:
            logging.error('WalltimeError thrown; calling times_up')
            self.times_up(self.process)
            raise err

        # All other errors mean something bad has happened, so we should clean up this job
        # more jobs will run unless this exception is either SystemExit or KeyboardInterrupt
        except Exception as err:
            logging.error('Error caught: terminating job for {}. Error = {}'.format(self.seed, err))
            self.process.terminate()
            self._finalise_result()
            raise err

    def scf(self, calc_doc, seed, keep=True, intermediate=False):
        """ Perform a single-shot calculation with CASTEP.

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
        logging.info('Performing single-shot CASTEP run on {}'.format(seed))
        try:
            self.cp_to_input(self.seed)

            if not self.custom_params:
                doc2param(calc_doc, seed, hash_dupe=False, overwrite=True)
            doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False, overwrite=True)

            # run CASTEP
            self.process = self.run_command(seed)
            self.process.communicate()

            # scrape dict but ignore the results
            results_dict, success = castep2dict(seed + '.castep', db=False)
            # check for errors
            errors_present, errors = self._catch_castep_errors()
            if errors_present:
                msg = 'CASTEP run on {} failed with errors: {}'.format(seed, errors)
                logging.error(msg)
                raise RuntimeError(msg)

            if not success:
                msg = 'Error scraping CASTEP file {}: {}'.format(seed, results_dict)
                logging.error(msg)
                raise RuntimeError(msg)

            success = True

            if not intermediate:
                logging.info('Writing results of singleshot CASTEP run to res file and tidying up.')
                doc2res(results_dict, seed, hash_dupe=False, overwrite=True)
                self.mv_to_completed(seed, keep=keep, completed_dir=self.paths['completed_dir'])
                if not keep:
                    self.tidy_up(seed)

            return success

        except Exception as err:
            logging.error('Caught error: {}'.format(err))
            self.mv_to_bad(seed)
            if not keep:
                self.tidy_up(seed)
            raise err

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

        for keyword in required:
            assert keyword in calc_doc

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
        logging.info('Getting seekpath cell and kpoint path.')
        from matador.utils.cell_utils import get_seekpath_kpoint_path
        from matador.crystal import Crystal
        logging.debug('Old lattice: {}'.format(Crystal(calc_doc)))
        prim_doc, kpt_path, _ = get_seekpath_kpoint_path(calc_doc,
                                                         spacing=spacing,
                                                         debug=debug)
        logging.debug('New lattice: {}'.format(Crystal(calc_doc)))
        return prim_doc, kpt_path

    def run_convergence_tests(self, calc_doc):
        """ Run kpoint and cutoff_energy convergence tests based on
        options passed to FullRelaxer.

        Parameters:
            calc_doc (dict): the structure to converge.

        Returns:
            bool: True unless every single calculation failed.

        """
        logging.info('Performing convergence tests...')
        from matador.utils.cell_utils import get_best_mp_offset_for_cell
        successes = []
        cached_cutoff = calc_doc['cut_off_energy']
        if self.conv_cutoff_bool:
            # run series of singlepoints for various cutoffs
            logging.info('Running cutoff convergence...')
            for cutoff in self.conv_cutoff:
                logging.info('{} eV... '.format(cutoff))
                calc_doc.update({'cut_off_energy': cutoff})
                self.paths['completed_dir'] = 'completed_cutoff'
                seed = self.seed + '_' + str(cutoff) + 'eV'
                success = self.scf(calc_doc, seed, keep=False)
                successes.append(success)
        if self.conv_kpt_bool:
            # run series of singlepoints for various cutoffs
            logging.info('Running kpt convergence tests...')
            calc_doc['cut_off_energy'] = cached_cutoff
            for kpt in self.conv_kpt:
                logging.info('{} 1/A... '.format(kpt))
                calc_doc.update({'kpoints_mp_spacing': kpt})
                calc_doc['kpoints_mp_offset'] = get_best_mp_offset_for_cell(calc_doc)
                logging.debug('Using offset {}'.format(calc_doc['kpoints_mp_offset']))
                self.paths['completed_dir'] = 'completed_kpts'
                seed = self.seed + '_' + str(kpt) + 'A'
                success = self.scf(calc_doc, seed, keep=False)
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

        logging.debug('Executable string parsed as {}'.format(command))

        return command

    def test_exec(self):
        """ Test if <executable> --version returns a valid string.

        Raises:
            SystemExit: if executable not found.

        """
        logging.info('Testing executable {executable}.'.format(executable=self.executable))
        try:
            proc = self.run_command('--version', exec_test=True)
        except FileNotFoundError:
            logging.critical('Unable to call mpirun/aprun/srun, currently selected: {}'.format(self.mpi_library))
            message = 'Please check initialistion of FullRelaxer object/CLI args.'
            logging.debug('Raising SystemExit with message: {message}'.format(message=message))
            raise SystemExit(message)

        out, errs = proc.communicate()
        if 'version' not in out.decode('utf-8') and errs is not None:
            err_string = 'Executable {exc} failed testing. Is it on your PATH?'.format(exc=self.executable)
            logging.critical(err_string)
            logging.critical('stdout: {stdout}'.format(stdout=out.decode('utf-8')))
            logging.critical('sterr: {stderr}'.format(stderr=errs.decode('utf-8')))
            logging.debug('Raising SystemExit.')
            raise SystemExit(err_string)

    @property
    def mpi_library(self):
        """ Property to store/compute desired MPI library. """
        if self._mpi_library is None:
            self._mpi_library = self.set_mpi_library()
            logging.info('Detected {mpi} MPI.'.format(mpi=self._mpi_library))
        return self._mpi_library

    def set_mpi_library(self):
        """ Combines command-line MPI arguments into string and calls
        MPI library detection is no args are present.
        """
        if sum([self.archer, self.intel, self.slurm]) > 1:
            message = 'Conflicting command-line arguments for MPI library have been supplied, exiting.'
            logging.critical(message)
            raise SystemExit(message)
        elif self.archer:
            return 'archer'
        elif self.intel:
            return 'intel'
        elif self.slurm:
            return 'slurm'
        else:
            return self.detect_mpi()

    @staticmethod
    def detect_mpi():
        """ Test which mpi library is being used when `mpirun`.

        Returns:
            mpi_library (str): 'intel', 'archer', or 'default'.

        """
        # check first for existence of mpirun command, then aprun if that fails
        try:
            try:
                logging.info('Failed to find mpirun, checking aprun...')
                mpi_version_string = str(sp.check_output('mpirun --version', shell=True))
            except sp.CalledProcessError:
                mpi_version_string = str(sp.check_output('aprun --version', shell=True))
        except Exception as exc:
            msg = 'Failed to find mpirun or aprun.'
            logging.critical(msg)
            logging.debug('Error message: {exc}'.format(exc=exc))
            raise SystemExit(msg)
        if 'Intel' in mpi_version_string:
            mpi_version = 'intel'
        elif 'aprun' in mpi_version_string:
            mpi_version = 'archer'
        elif 'Open MPI' in mpi_version_string:
            mpi_version = 'default'
        else:
            logging.debug('Could not detect MPI library so using default (OpenMPI), version string was: {response}'
                          .format(response=mpi_version_string))
            mpi_version = 'default'

        logging.info('Using {version} MPI library.'.format(version=mpi_version))
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
        logging.info('Performing memory check for {seed}'.format(seed=seed))
        memcheck_seed = seed + '_memcheck'
        doc2param(calc_doc, memcheck_seed, hash_dupe=False)
        doc2cell(calc_doc, memcheck_seed, hash_dupe=False, copy_pspots=False)
        free_memory = float(virtual_memory().available) / 1024**2
        if self.maxmem is None:
            maxmem = 0.9 * free_memory
        else:
            maxmem = self.maxmem

        # check if cell is totally pathological, as CASTEP dryrun will massively underestimate mem
        if all([angle < 30 for angle in calc_doc['lattice_abc'][1]]):
            logging.error('Cell is pathological (at least one angle < 30), failing memory check.')
            return False

        logging.debug('Running CASTEP dryrun.')
        process = sp.Popen(['nice', '-n', '15', self.executable, '-d', memcheck_seed])
        process.communicate()

        skip = False
        results, success = castep2dict(memcheck_seed + '.castep', db=False)
        if not success or 'estimated_mem_MB' not in results:
            skip = True
            logging.warning('CASTEP dryrun failed with output {results}'.format(results=results))

        for _file in glob.glob(memcheck_seed + '*'):
            if _file.endswith('.res'):
                continue
            else:
                os.remove(_file)

        if 'estimated_mem_MB' in results:
            mem_string = ('Memory estimate / Available memory (MB): {:8.0f} / {:8.0f}'
                          .format(results['estimated_mem_MB'], maxmem))
            logging.info(mem_string)

        if skip or results['estimated_mem_MB'] > maxmem:
            return False

        return True

    def run_command(self, seed, exec_test=False):
        """ Calls executable on seed with desired number of cores.

        Parameters:
            seed (str): seedname to pass append to CASTEP command,
                e.g. <seed> or --version.

        Keyword arguments:
            exec_test (bool): run executable in test mode, with output
                piped to stdout.

        Returns:
            subprocess.Popen: process to run.

        """
        command = self.parse_executable(seed)
        if self.nnodes is None or self.nnodes == 1:
            if self.ncores == 1 and self.node is None:
                command = ['nice', '-n', '15'] + command
            elif self.mpi_library == 'archer':
                command = ['aprun', '-n', str(self.ncores)] + command
            elif self.mpi_library == 'slurm':
                command = ['srun', '--exclusive', '-N', '1', '-n', str(self.ncores)] + command
            elif self.mpi_library == 'intel':
                command = ['mpirun', '-n', str(self.ncores), '-ppn', str(self.ncores)] + command
            elif self.node is not None:
                cwd = os.getcwd()
                command = ['ssh', '{}'.format(self.node), 'cd', '{};'.format(cwd), 'mpirun', '-n',
                           str(self.ncores)] + command
            else:
                command = ['nice', '-n', '15', 'mpirun', '-n', str(self.ncores)] + command
        else:
            if self.mpi_library == 'archer':
                command = [
                    'aprun', '-n',
                    str(self.ncores * self.nnodes), '-N',
                    str(self.ncores), '-S', '12', '-d', '1'
                ] + command
            elif self.mpi_library == 'slurm':
                command = ['srun', '--exclusive', '-N',
                           str(self.nnodes), '-n',
                           str(self.ncores * self.nnodes)] + command
            elif self.mpi_library == 'intel':
                command = ['mpirun', '-n', str(self.ncores * self.nnodes), '-ppn', str(self.ncores)] + command
            else:
                command = ['mpirun', '-n', str(self.ncores * self.nnodes), '-npernode', str(self.ncores)] + command

        # ensure default stdout is None for generic mpi calls
        stdout = None
        stderr = None

        if exec_test:
            stdout = sp.PIPE
            stderr = sp.PIPE
        elif self.debug:
            stdout = None
            stderr = None
        else:
            dev_null = open(os.devnull, 'w')
            stdout = dev_null
            stderr = dev_null

        if self._redirect_filename is not None:
            logging.info('Redirecting output to {redirect}'.format(redirect=self._redirect_filename))
            redirect_file = open(self._redirect_filename, 'w')
            stdout = redirect_file

        process = sp.Popen(command, shell=False, stdout=stdout, stderr=stderr)
        try:
            redirect_file.close()
        except Exception:
            pass
        try:
            dev_null.close()
        except Exception:
            pass

        if not exec_test:
            logging.debug('Running {}'.format(command))
        return process

    def _catch_castep_errors(self):
        """ Look for CASTEP error files and fallover appropriately. If
        the magic string 'Work-around was succesful' is found in error,
        no errors will be reported and the file will be deleted,
        provided no other errors exist.

        Returns:
            bool, str: True if error files were found, otherwise False,
                followed by the error messages.

        """
        err_file = '{}*err'.format(self.seed)
        error_str = ''
        errors_present = False

        for globbed in glob.glob(err_file):
            if os.path.isfile(globbed):
                with open(globbed, 'r') as f:
                    flines = f.readlines()

                for line in flines:
                    if 'Work-around was successful, continuing with calculation.' in line:
                        logging.info('Found LAPACK issue that was circumvented, removing error file.')
                        os.remove(globbed)
                        break
                else:
                    error_str += ' '.join(flines)
                    error_str += '\n'
                    errors_present = True

        return errors_present, error_str

    def mv_to_bad(self, seed):
        """ Move all files associated with "seed" to bad_castep, from both
        the compute directory (if it exists) and the root dir..

        Parameters:
            seed (str): filename of structure.

        """
        try:
            bad_dir = self.root_folder + '/bad_castep'
            logging.info('Moving files to bad_castep: {bad}.'.format(bad=bad_dir))
            if not os.path.exists(bad_dir):
                os.makedirs(bad_dir, exist_ok=True)
            seed_files = glob.glob(seed + '.*')
            if seed_files:
                logging.debug('Files to move: {seed}'.format(seed=seed_files))
                for _file in seed_files:
                    try:
                        shutil.copy(_file, bad_dir)
                        os.remove(_file)
                    except Exception as exc:
                        logging.warning('Error moving files to bad: {error}'.format(error=exc))
            # check root folder for any matching files and remove them
            fname = '{}/{}'.format(self.root_folder, seed)
            for ext in ['.res', '.res.lock', '.castep']:
                if os.path.isfile('{}{}'.format(fname, ext)):
                    os.remove('{}{}'.format(fname, ext))
        except Exception as exc:
            logging.warning('Error moving files to bad: {error}'.format(error=exc))

    def mv_to_completed(self, seed, completed_dir='completed', keep=False):
        """ Move all associated files to completed, removing any
        remaining files in the root_folder and compute_dir.

        Parameters:
            seed (str): filename for structure.

        Keyword arguments:
            completed_dir (str): folder for completed jobs.
            keep (bool): whether to also move intermediate files.


        """
        completed_dir = self.root_folder + '/' + completed_dir
        logging.info('Moving files to completed: {completed}.'.format(completed=completed_dir))
        if not os.path.exists(completed_dir):
            os.makedirs(completed_dir, exist_ok=True)
        if keep:
            seed_files = glob.glob(seed + '.*') + [seed + '-out.cell']
            if seed_files:
                logging.debug('Files to move: {files}.'.format(files=seed_files))
                for _file in seed_files:
                    shutil.move(_file, completed_dir)
        else:
            # move castep/param/res/out_cell files to completed
            file_exts = ['.castep']
            if self.kpts_1D:
                file_exts.append('.param')
            if not self.conv_kpt_bool and not self.conv_cutoff_bool:
                file_exts.append('.res')
            if self.calc_doc.get('write_cell_structure'):
                file_exts.append('-out.cell')
            if self.calc_doc.get('write_formatted_density'):
                file_exts.append('.den_fmt')
            for ext in file_exts:
                try:
                    shutil.move('{}{}'.format(seed, ext), completed_dir)
                except Exception as exc:
                    logging.warning('Error moving files to completed: {error}'.format(error=exc))
            # check root folder for any matching files and remove them
            fname = '{}/{}'.format(self.root_folder, seed)
            for ext in file_exts + ['.res.lock']:
                if os.path.isfile('{}{}'.format(fname, ext)):
                    os.remove('{}{}'.format(fname, ext))

            wildcard_fnames = glob.glob('{}/{}.*'.format(self.root_folder, seed))
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
        input_dir = self.root_folder + '/input'
        logging.debug('Copying file to input_dir: {input}'.format(input=input_dir))
        if not os.path.exists(input_dir):
            os.makedirs(input_dir, exist_ok=True)
        if glob_files:
            files = glob.glob('{}*'.format(seed))
            logging.debug('Files to copy: {files}'.format(files=files))
            for f in files:
                if f.endswith('.lock'):
                    continue
                if not os.path.isfile(f):
                    shutil.copy('{}'.format(f), input_dir)
        else:
            logging.debug('File to copy: {file}'.format(file='{}.{}'.format(seed, ext)))
            if os.path.isfile('{}.{}'.format(seed, ext)):
                if not os.path.isfile('{}/{}.{}'.format(input_dir, seed, ext)):
                    shutil.copy('{}.{}'.format(seed, ext), input_dir)

    def _setup_relaxation(self):
        """ Set up directories and files for relaxation. """

        logging.info('Preparing to relax {seed}'.format(seed=self.seed))
        if self.compute_dir is not None:
            logging.info('Using compute_dir: {compute}'.format(compute=self.compute_dir))
            if not os.path.isdir(self.compute_dir):
                os.makedirs(self.compute_dir)

            # copy pspots and any intermediate calcs to compute_dir
            logging.info('Copying pspots into compute_dir')
            pspots = glob.glob('*.usp')
            for pspot in pspots:
                shutil.copy(pspot, self.compute_dir)

        # update res file with intermediate calculation if castep file is newer than res
        if os.path.isfile(self.seed + '.castep') and os.path.isfile(self.seed + '.res'):
            logging.info('Trying to update res file with result from intermediate CASTEP file found in root_dir')
            if self.compute_dir is not None:
                shutil.copy(self.seed + '.castep', self.compute_dir)
            castep_dict, success = castep2dict(self.seed + '.castep', db=False)
            if success:
                self.res_dict['geom_iter'] = castep_dict.get('geom_iter')

            if os.path.getmtime(self.seed + '.res') < os.path.getmtime(self.seed + '.castep'):
                logging.info('CASTEP file was updated more recently than res file, using intermediate structure...')
                self.res_dict.update(castep_dict)

        if self.compute_dir is not None:
            os.chdir(self.compute_dir)

        # copy initial res file to seed
        logging.info('Writing fresh res file to start calculation from.')
        doc2res(self.res_dict, self.seed, info=False, hash_dupe=False, overwrite=True)
        self.cp_to_input(self.seed)

        # set up geom opt iteration options based on input/scraped parameters
        self.max_iter = self.calc_doc['geom_max_iter']
        if 'geom_iter' in self.res_dict:
            self.max_iter -= self.res_dict['geom_iter']

        if self.max_iter < 0:
            msg = '{} iterations already performed on structure, exiting...'
            logging.critical(msg)
            raise RuntimeError(msg)

        # number of steps in fine and rough calcs respectively
        fine_iter = self.fine_iter
        rough_iter = self.rough_iter
        # number of calculations to run with rough_iter steps
        num_rough_iter = self.rough
        if 'geom_method' in self.calc_doc:
            if self.calc_doc['geom_method'].lower() == 'tpsd' and rough_iter < 3:
                rough_iter = 3
        self.geom_max_iter_list = (num_rough_iter * [rough_iter])
        self.max_iter -= num_rough_iter * rough_iter

        num_fine_iter = ceil(int(self.max_iter) / fine_iter)
        if self.max_iter > 0:
            if self.max_iter < fine_iter:
                fine_iter = self.max_iter
                num_fine_iter = 1
            self.geom_max_iter_list.extend(num_fine_iter * [fine_iter])

        logging.info('Geometry optimisation iteration scheme set to {}'.format(self.geom_max_iter_list))

        if not self.geom_max_iter_list:
            msg = 'Could not divide up relaxation; consider increasing geom_max_iter'
            logging.critical(msg)
            raise SystemExit(msg)

    def _update_input_files(self):
        """ Update the cell and param files for the next relaxation. """
        calc_doc = self.calc_doc
        # update cell
        if os.path.isfile(self.seed + '.cell'):
            os.remove(self.seed + '.cell')
        if self.kpts_1D:
            n_kz = ceil(1 / (calc_doc['lattice_abc'][0][2] * self._target_spacing))
            if n_kz % 2 == 1:
                n_kz += 1
            calc_doc['kpoints_mp_grid'] = [1, 1, n_kz]
            if 'kpoints_mp_spacing' in calc_doc:
                del calc_doc['kpoints_mp_spacing']
        doc2cell(calc_doc, self.seed, hash_dupe=False, copy_pspots=False, spin=self.spin)

        # update param
        if not self.custom_params:
            if os.path.isfile(self.seed + '.param'):
                os.remove(self.seed + '.param')
        doc2param(calc_doc, self.seed, hash_dupe=False, spin=self.spin)

    @staticmethod
    def tidy_up(seed):
        """ Delete all created files before quitting.

        Parameters:
            seed (str): filename for structure.

        """
        files = glob.glob(seed + '.*')
        if files:
            logging.info('Tidying up remaining files: {files}'.format(files=files))
            for f in files:
                if not (f.endswith('.res') or f.endswith('.castep')):
                    os.remove(f)

    def _update_output_files(self, opti_dict):
        """ Copy new data to output files and update
         the results dict.

        Parameter:
            opti_dict (dict): intermediate calculation results.

        """
        logging.info('Updating .res and .castep files in root_dir with new results')
        if os.path.isfile(self.seed + '.res'):
            os.remove(self.seed + '.res')
        doc2res(opti_dict, self.seed, hash_dupe=False)
        if self.compute_dir is not None:
            shutil.copy(self.seed + '.res', self.root_folder)
            if os.path.isfile(self.seed + '.castep'):
                shutil.copy(self.seed + '.castep', self.root_folder)
        self.res_dict.update(opti_dict)

    def _finalise_result(self):
        """ Push to queue if necessary and return status.

        Returns:
            bool: True is relaxation was successful, False otherwise.

        """
        logging.info('Finalising calculation...')
        try:
            success = self.res_dict.get('optimised', False)
        except AttributeError:
            success = False

        logging.info('Was calculation successful? {success}'.format(success=success))
        if self.output_queue is not None:
            logging.info('Pushing results to output queue')
            self.output_queue.put(self.res_dict)
        if success:
            self.mv_to_completed(self.seed, completed_dir=self.paths['completed_dir'])
        else:
            self.mv_to_bad(self.seed)

        # clean up rest of files
        self.tidy_up(self.seed)

        return success

    def times_up(self, process):
        """ If walltime has nearly expired, run this function
        to kill the process and unlock it for restarted calculations.

        Parameters:
            subprocess.Popen: running process to be killed.

        """
        logging.info('Ending process early for seed: {seed}'.format(seed=self.seed))
        process.terminate()
        logging.info('Ended process early for seed: {seed}'.format(seed=self.seed))
        if self.compute_dir is not None:
            logging.info('Cleaning up compute_dir: {dir}'.format(dir=self.compute_dir))
            for f in glob.glob('{}.*'.format(self.seed)):
                shutil.copy(f, self.root_folder)
                os.remove(f)

        logging.info('Removing lock file so calculation can be continued.')
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
        logging.info('Checking if compute_dir still contains calculations...')
        if not os.path.isdir(compute_dir):
            return False

        files = glob.glob(compute_dir + '/*')
        logging.debug('Found {files} in {dir}'.format(files=files, dir=compute_dir))

        for fname in files:
            if fname.endswith('.res') or fname.endswith('.castep'):
                logging.debug('Not removing {dir} as it still contains calculation {fname}'.format(
                    dir=compute_dir, fname=fname))
                return False

        # remove files in directory, then delete directory
        logging.debug('Deleting files {files} from {dir}'.format(files=files, dir=compute_dir))
        for fname in files:
            if os.path.isfile(fname):
                try:
                    os.remove(fname)
                except FileNotFoundError:
                    pass

        if os.path.isdir(compute_dir):
            logging.debug('Deleting directory {dir}'.format(dir=compute_dir))
            os.rmdir(compute_dir)
        return True


class WalltimeError(Exception):
    """ Raise this when you don't want any more jobs to run
    because they're about to exceed the max walltime.

    """
    pass
