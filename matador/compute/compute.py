# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the FullRelaxer class for handling
calculations on a single structure, and some useful associated
errors.

"""


import os
import shutil
import subprocess as sp
import glob
import time
from copy import deepcopy
from traceback import print_exc
from math import ceil
from psutil import virtual_memory

from matador.scrapers.castep_scrapers import cell2dict
from matador.scrapers.castep_scrapers import res2dict, castep2dict
from matador.utils.print_utils import print_success, print_warning, print_notify, print_failure
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
            ncores (int): number of cores for mpirun call
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
            timings (`obj`:tuple: of `obj`:int:): tuple containing max and elapsed time in seconds

        Raises:
            WalltimeError: if desired walltime is exceeded.
            SystemExit: if a fatal error occurs.
            RuntimeError: if a structure-level error occurs.

        """
        # set defaults and update class with desired values
        prop_defaults = {'paths': None, 'param_dict': None, 'cell_dict': None, 'mode': 'castep', 'executable': 'castep',
                         'memcheck': False, 'rough': 4, 'rough_iter': 2, 'fine_iter': 20, 'spin': False, 'output_queue': None,
                         'redirect': None, 'reopt': False, 'compute_dir': None, 'custom_params': False, 'archer': False,
                         'maxmem': None, 'killcheck': True, 'kpts_1D': False, 'conv_cutoff': False, 'conv_kpt': False,
                         'debug': False, 'profile': False, 'slurm': False, 'intel': False, 'exec_test': True, 'timings': (None, None),
                         'start': True, 'verbosity': 0, 'polltime': 30}

        self.paths = None
        self.process = None
        self.output_queue = None
        self.__dict__.update(prop_defaults)
        self.__dict__.update(kwargs)

        if self.profile:
            import cProfile
            import pstats
            from sys import version_info
            from matador import __version__
            profile = cProfile.Profile()
            profile.enable()

        self.ncores = ncores
        self.res = res
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

        if self.verbosity > 2:
            print('Timing data: ', self.timings)

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

        except RuntimeError:
            self.success = False

        if self.profile:
            profile.disable()
            fname = 'relaxer-{}-{}-{}.{}.{}'.format(__version__, os.hostname()[1], version_info.major,
                                                    version_info.minor, version_info.micro)
            profile.dump_stats(fname + '.prof')
            with open(fname + '.pstats', 'w') as fp:
                stats = pstats.Stats(profile, stream=fp).sort_stats('cumulative')
                stats.print_stats()

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
        if self.exec_test:
            self.test_exec()

        if self.kpts_1D:
            if 'kpoints_mp_spacing' not in self.cell_dict:
                raise SystemExit('kpoints_mp_spacing not found, but kpts_1D requested...')
            self._target_spacing = deepcopy(self.cell_dict['kpoints_mp_spacing'])

        # read in initial structure and skip if failed
        if isinstance(res, str):
            self.res_dict, success = res2dict(res, db=False)
            if not success:
                raise RuntimeError('Failed to parse initial res file {}'.format(res))

        elif isinstance(res, dict):
            self.res_dict = res

        calc_doc = deepcopy(self.res_dict)

        # set seed name
        assert isinstance(calc_doc['source'], list)
        self.seed = calc_doc['source'][0].replace('.res', '')

        # update global doc with cell and param dicts for folder
        calc_doc.update(self.cell_dict)
        calc_doc.update(self.param_dict)

        # check for pseudos
        if 'library' not in calc_doc['species_pot']:
            for elem in self.res_dict['stoichiometry']:
                if ('|' not in calc_doc['species_pot'][elem[0]] and
                        not os.path.isfile(calc_doc['species_pot'][elem[0]])):
                    raise SystemExit('You forgot your pseudos, you silly goose!')

        # this is now a dict containing the exact calculation we are going to run
        self.calc_doc = calc_doc

        # do memcheck, if desired, and only continue if enough memory is free
        if self.memcheck:
            self.enough_memory = self.do_memcheck(calc_doc, self.seed)
        else:
            self.enough_memory = True
        if not self.enough_memory:
            return self.enough_memory

        # run convergence tests
        if any([self.conv_cutoff_bool, self.conv_kpt_bool]):
            if self.verbosity > 0:
                print('Running convergence tests...')
            success = self.run_convergence_tests(calc_doc)
            return success

        # perform relaxation
        elif calc_doc['task'].upper() in ['GEOMETRYOPTIMISATION', 'GEOMETRYOPTIMIZATION']:
            # set up geom opt parameters
            self.max_iter = calc_doc['geom_max_iter']
            self.num_rough_iter = self.rough
            fine_iter = self.fine_iter
            rough_iter = self.rough_iter
            if 'geom_method' in calc_doc:
                if calc_doc['geom_method'].lower() == 'tpsd' and rough_iter < 3:
                    rough_iter = 3
            num_fine_iter = int(int(self.max_iter) / fine_iter)
            self.geom_max_iter_list = (self.num_rough_iter * [rough_iter])
            self.geom_max_iter_list.extend(num_fine_iter * [fine_iter])
            if not self.geom_max_iter_list:
                raise SystemExit('Could not divide up relaxation; consider increasing geom_max_iter')

            # begin relaxation
            if self.start:
                try:
                    success = self.relax()
                except Exception as err:
                    if self.compute_dir is not None:
                        # always cd back to root folder
                        os.chdir(self.root_folder)
                        self.remove_compute_dir_if_finished(self.compute_dir, debug=self.debug)
                    raise err

                if self.compute_dir is not None:
                    # always cd back to root folder
                    os.chdir(self.root_folder)
                    self.remove_compute_dir_if_finished(self.compute_dir, debug=self.debug)

                return success

        elif calc_doc['task'].upper() in ['PHONON', 'THERMODYNAMICS']:
            success = self.castep_full_phonon(calc_doc, self.seed)

        # run in SCF mode, i.e. just call CASTEP on the seeds
        else:
            success = self.scf(calc_doc, self.seed, keep=True)
            return success

    def run_generic(self, seed):
        """ Run a generic mpi program on the given seed. Files from
        completed runs are moved to "completed" and failed runs to "bad_castep".

        Parameters:
            seed (str): filename of structure

        Returns:
            bool: True if calculations progressed without error.

        """
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
                self.mv_to_bad(seed)
                return False

            self.mv_to_completed(seed, keep=True, completed_dir=self.paths['completed_dir'])
            return True

        except Exception as err:
            print_exc()
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
        self._setup_relaxation_dirs()
        seed = self.seed
        rerun = False
        try:
            # iterate over geom iter blocks
            for ind, num_iter in enumerate(self.geom_max_iter_list):

                # print some blurb about what we're doing
                if self.verbosity > 0:
                    if ind == 0:
                        if self.verbosity > 2:
                            if self.custom_params is not False:
                                print_notify('custom params: {}'.format(self.custom_params))
                        print_notify('Beginning geometry optimisation...')
                    elif ind == self.num_rough_iter:
                        print_notify('Beginning fine geometry optimisation...')

                # preprocess last step that did not finish geom opt
                if self.reopt and rerun:
                    # if we're now reoptimising after a success in last step
                    # use the fine iter value
                    num_iter = self.fine_iter
                    if self.verbosity > 0:
                        print_notify('Last step was successful, performing one last relaxation...')

                # update the geom_max_iter to use with either the number in iter_list, or the overriden value
                self.calc_doc['geom_max_iter'] = num_iter

                # delete any existing files and write new ones
                self._update_input_files()

                # run CASTEP
                self.process = self.run_command(seed)

                # if specified max_walltime (or found through SLURM), then monitor job
                if self.max_walltime is not None:
                    if self.start_time is None:
                        raise SystemExit('Somehow initial start time was not found')

                    if self.debug:
                        print('Polling process every {} s'.format(self.polltime))

                    while self.process.poll() is None:
                        elapsed = time.time() - self.start_time
                        if self.debug:
                            print('Elapsed time: {:>10.1f} s'.format(elapsed))

                        # leave 1 minute to clean up
                        if elapsed > abs(self.max_walltime):  # - 60):
                            raise WalltimeError('Ran out of time on seed {}'.format(self.seed))
                        time.sleep(self.polltime)

                self.process.communicate()

                # scrape new structure from castep file
                if not os.path.isfile(seed + '.castep'):
                    raise SystemExit('CASTEP file was not created, '
                                     'please check your executable: {}.'.format(self.executable))

                opti_dict, success = castep2dict(seed + '.castep', db=False, verbosity=self.verbosity)
                if not success and isinstance(opti_dict, str):
                    raise RuntimeError('Failed to parse CASTEP file...')

                if self.debug:
                    print_notify('Intermediate calculation finished')
                    print(opti_dict)

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
                    if self.verbosity > 0:
                        print_success('Successfully relaxed ' + seed)
                    self._update_output_files(opti_dict)
                    return self._finalise_result()

                # reached maximum number of steps
                elif ind == len(self.geom_max_iter_list) - 1:
                    if self.verbosity > 0:
                        print_warning('Failed to optimise ' + seed)
                    return self._finalise_result()

                errors_present, errors = self._catch_castep_errors()
                if errors_present:
                    if self.verbosity > 0:
                        print_warning('Failed to optimise {} as CASTEP crashed with error:'.format(seed))
                        print(errors)

                    self._update_output_files(opti_dict)
                    return self._finalise_result()

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

                if self.debug:
                    print_notify('Restarting calculation with current state:')
                    print(self.calc_doc)
                if self.verbosity > 1:
                    print(('num_iter: {:3d} | max F: {:5f} eV/A | stress: {: 5f} GPa | ' +
                           'cell volume: {:5f} A^3 | enthalpy per atom {:5f} eV')
                          .format(sum(self.geom_max_iter_list[:ind+1]),
                                  opti_dict['max_force_on_atom'],
                                  opti_dict['pressure'],
                                  opti_dict['cell_volume'],
                                  opti_dict['enthalpy_per_atom']))

                self.calc_doc.update(opti_dict)

        # catch WalltimeErrors and reset the job folder ready for continuation
        except WalltimeError as err:
            self.times_up(self.process)
            raise err

        # All other errors mean something bad has happened, so we should clean up this job
        # more jobs will run unless this exception is either SystemExit or KeyboardInterrupt
        except Exception as err:
            self.process.terminate()
            self._finalise_result()
            raise err

    def castep_full_phonon(self, calc_doc, seed):
        """ Perform a "full" phonon calculation on a system, i.e.
        first perform a relaxation in a standardised unit cell,
        then compute the dynamical matrix, then finally interpolate
        that dynamical matrix into dispersion curves and DOS.

        Parameters:
            calc_doc (dict): dictionary of structure and calculation
                parameters.
            seed (str): root seed for the calculation.

        Raises:
            RuntimeError: if any part of the calculation fails.

        """
        from matador.utils.cell_utils import cart2abc
        todo = {'relax': True, 'dynmat': True, 'dispersion': False, 'dos': False, 'thermodynamics': False}
        if calc_doc.get('task').lower() in ['phonon', 'thermodynamics']:
            if (('phonon_fine_kpoint_path' not in calc_doc and 'phonon_fine_kpoint_list' not in calc_doc) and
                    'phonon_fine_kpoint_path_spacing' in calc_doc):
                todo['dispersion'] = True
            if 'phonon_kpoint_mp_spacing' in calc_doc:
                todo['dos'] = True
            if calc_doc['task'].lower() == 'thermodynamics':
                todo['thermodynamics'] = True

        # always standardise the cell so that any phonon calculation can have
        # post-processing performed after the fact
        prim_doc, kpt_path = self._get_seekpath_compliant_input(
            calc_doc, calc_doc.get('phonon_fine_kpoint_path_spacing', 0.02))
        calc_doc.update(prim_doc)
        calc_doc['lattice_abc'] = cart2abc(calc_doc['lattice_cart'])

        if todo['dispersion']:
            calc_doc['phonon_fine_kpoint_list'] = kpt_path

        # always shift phonon grid
        if 'phonon_kpoint_mp_spacing' in calc_doc:
            from matador.utils.cell_utils import calc_mp_grid, shift_to_include_gamma
            grid = calc_mp_grid(calc_doc['lattice_cart'], calc_doc['phonon_kpoint_mp_spacing'])
            offset = shift_to_include_gamma(grid)
            if offset != [0, 0, 0]:
                calc_doc['phonon_kpoint_mp_offset'] = offset
                if self.verbosity > 2:
                    print('Set phonon MP grid offset to {}'.format(offset))

        # prepare to do pre-relax if there's no check file
        if os.path.isfile(seed + '.check'):
            todo['relax'] = False
            if self.verbosity > 0:
                print_notify('Restarting from {}.check, so not performing re-relaxation'.format(seed))

        if self.debug:
            print(todo)

        if todo['relax']:
            if self.verbosity > 0:
                print('Pre-relaxing structure...')
            success = self._castep_phonon_prerelax_only(calc_doc, seed,
                                                        intermediate=bool(sum([bool(todo[key]) for key in todo])))
            if success:
                todo['relax'] = False
                if self.verbosity > 0:
                    print_notify('Pre-relaxation complete.')
            else:
                raise RuntimeError('Pre-requisite geometry optimisation failed.')

        if todo['dynmat']:
            if self.verbosity > 0:
                print('Now computing phonon dynmat...')
            success = self._castep_phonon_dynmat_only(calc_doc, seed,
                                                      intermediate=bool(sum([bool(todo[key]) for key in todo])))
            if success:
                todo['dynmat'] = False
                if self.verbosity > 0:
                    print_notify('Dynmat complete.')
            else:
                raise RuntimeError('Phonon dynamical matrix calculation failed.')

        if todo['thermodynamics']:
            if self.verbosity > 0:
                print('Now performing phonon thermodynamics...')
            success = self._castep_phonon_thermodynamics_only(calc_doc, seed,
                                                              intermediate=bool(sum([bool(todo[key]) for key in todo])))
            if success:
                todo['thermodynamics'] = False
                if self.verbosity > 0:
                    print_notify('Thermodynamics complete.')
            else:
                raise RuntimeError('Phonon thermodynamics calculation failed.')

        if todo['dos']:
            if self.verbosity > 0:
                print('Now performing phonon DOS...')
            success = self._castep_phonon_dos_only(calc_doc, seed,
                                                   intermediate=bool(sum([bool(todo[key]) for key in todo])))
            if success:
                todo['dos'] = False
                if self.verbosity > 0:
                    print_notify('DOS complete.')
            else:
                raise RuntimeError('Phonon DOS calculation failed.')

        if todo['dispersion']:
            if self.verbosity > 0:
                print('Now performing phonon dispersion...')
            success = self._castep_phonon_dispersion_only(calc_doc, seed,
                                                          intermediate=bool(sum([bool(todo[key]) for key in todo])))
            if success:
                todo['dispersion'] = False
                if self.verbosity > 0:
                    print_notify('Dispersion complete.')
            else:
                raise RuntimeError('Phonon dispersion calculation failed.')

    def _get_seekpath_compliant_input(self, calc_doc, spacing):
        """ Return seekpath cell/kpoint path for the given cell and spacing.

        Parameters:
            calc_doc (dict): structural and calculation parameters.
            spacing (float): desired kpoint path spacing.

        Returns:
            (dict, list): dictionary containing the standardised unit cell
                and list containing the kpoints.

        """
        from matador.utils.cell_utils import get_seekpath_kpoint_path
        if self.verbosity >= 2:
            from matador.crystal import Crystal
            print('Old lattice:')
            print(Crystal(calc_doc))
        prim_doc, kpt_path, _ = get_seekpath_kpoint_path(calc_doc,
                                                         spacing=spacing,
                                                         debug=self.debug)

        if self.verbosity >= 2:
            print('New lattice:')
            print(Crystal(calc_doc))

        return prim_doc, kpt_path

    def scf(self, calc_doc, seed, keep=True, intermediate=False):
        """ Perform a single-shot (e.g. scf, bandstructure, NMR)
        calculation with CASTEP.

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
        from matador.utils.cell_utils import cart2abc
        try:
            self.cp_to_input(self.seed)

            # try to add a k/q-point path to cell, for spectral/phonon tasks
            elec_dispersion = False
            if 'spectral_task' in calc_doc and calc_doc['spectral_task'] == 'bandstructure':
                if 'spectral_kpoints_path' not in calc_doc and 'spectral_kpoints_list' not in calc_doc:
                    elec_dispersion = True

            if elec_dispersion:
                prim_doc, kpt_path = self._get_seekpath_compliant_input(
                    calc_doc, spacing=calc_doc.get('spectral_kpoints_path_spacing', 0.02))
                calc_doc.update(prim_doc)
                calc_doc['lattice_abc'] = cart2abc(calc_doc['lattice_cart'])
                calc_doc['spectral_kpoints_list'] = kpt_path

            if not self.custom_params:
                doc2param(calc_doc, seed, hash_dupe=False, overwrite=True)
            doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False, overwrite=True)

            # run CASTEP
            self.process = self.run_command(seed)
            self.process.communicate()

            # scrape dict but ignore the results
            _, success = castep2dict(seed + '.castep', db=False)
            # check for errors
            errors_present, errors = self._catch_castep_errors()
            if errors_present:
                if self.verbosity > 1:
                    print_warning('Failed to optimise {} as CASTEP crashed with error:'.format(seed))
                    print(errors)
                raise RuntimeError('SCF failed.')

            if not success:
                raise RuntimeError('Failed to scrape CASTEP file.')

            success = True

            if not intermediate:
                self.mv_to_completed(seed, keep=keep, completed_dir=self.paths['completed_dir'])
                if not keep:
                    self.tidy_up(seed)

            return success

        except Exception as err:
            if self.verbosity > 1:
                print_exc()
            self.mv_to_bad(seed)
            if not keep:
                self.tidy_up(seed)
            raise err

    @staticmethod
    def _validate_calc_doc(calc_doc, required, forbidden):
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

    def _castep_phonon_prerelax_only(self, calc_doc, seed, intermediate=False):
        """ Run a singleshot geometry optimisation before an SCF-style calculation.
        This is typically used to ensure phonon calculations start successfully.
        The phonon calculation will then be restarted from the .check file produced here.

        Parameters:
            calc_doc (dict): the structure to converge.
            seed (str): root filename of structure.

        Keyword arguments:
            final (bool): whether this is the final step in a calculation.

        """
        relax_doc = deepcopy(calc_doc)
        relax_doc['write_checkpoint'] = 'ALL'
        if 'geom_max_iter' not in relax_doc:
            relax_doc['geom_max_iter'] = 20
        relax_doc['task'] = 'geometryoptimisation'

        required = []
        forbidden = ['phonon_fine_kpoint_list',
                     'phonon_fine_kpoint_path',
                     'phonon_fine_kpoint_mp_spacing',
                     'phonon_fine_kpoint_path_spacing']

        self._validate_calc_doc(relax_doc, required, forbidden)

        return self.scf(relax_doc, seed, keep=True, intermediate=intermediate)

    def _castep_phonon_dynmat_only(self, calc_doc, seed, intermediate=False):
        """ Runs a singleshot phonon dynmat calculation, with no "fine_method" interpolation.

        Parameters:
            calc_doc (dict): the structure to converge.
            seed (str): root filename of structure.

        Keyword arguments:
            final (bool): whether this is the final step in a calculation.

        """
        relax_doc = deepcopy(calc_doc)
        relax_doc['write_checkpoint'] = 'ALL'
        relax_doc['continuation'] = 'default'
        relax_doc['task'] = 'phonon'

        required = []
        forbidden = ['phonon_fine_kpoint_list',
                     'phonon_fine_kpoint_path',
                     'phonon_fine_kpoint_mp_spacing',
                     'phonon_fine_kpoint_path_spacing']

        self._validate_calc_doc(relax_doc, required, forbidden)
        return self.scf(relax_doc, seed, keep=True, intermediate=intermediate)

    def _castep_phonon_dos_only(self, calc_doc, seed, intermediate=False):
        """ Runs a DOS interpolation on top of a completed
        phonon calculation.

        Parameters:
            calc_doc (dict): the structure to converge.
            seed (str): root filename of structure.

        Keyword arguments:
            final (bool): whether this is the final step in a calculation.

        """
        dos_doc = deepcopy(calc_doc)
        dos_doc['task'] = 'phonon'
        dos_doc['phonon_calculate_dos'] = True

        required = ['phonon_fine_kpoint_mp_spacing']
        forbidden = ['phonon_fine_kpoint_list',
                     'phonon_fine_kpoint_path',
                     'phonon_fine_kpoint_path_spacing']

        self._validate_calc_doc(dos_doc, required, forbidden)

        return self.scf(dos_doc, seed, keep=True, intermediate=intermediate)

    def _castep_phonon_dispersion_only(self, calc_doc, seed, intermediate=False):
        """ Runs a dispersion interpolation on top of a completed
        phonon calculation.

        Parameters:
            calc_doc (dict): the structure to converge.
            seed (str): root filename of structure.

        Keyword arguments:
            final (bool): whether this is the final step in a calculation.

        """
        disp_doc = deepcopy(calc_doc)
        disp_doc['task'] = 'phonon'
        disp_doc['phonon_calculate_dos'] = False

        required = ['phonon_fine_kpoint_list']
        forbidden = ['phonon_fine_kpoint_mp_spacing',
                     'phonon_fine_kpoint_path',
                     'phonon_fine_kpoint_path_spacing']

        self._validate_calc_doc(disp_doc, required, forbidden)

        return self.scf(disp_doc, seed, keep=True, intermediate=intermediate)

    def _castep_phonon_thermodynamics_only(self, calc_doc, seed, intermediate=False):
        """ Runs a "thermodynamics" interpolation on top of a completed
        phonon calculation, using the phonon_fine_kpoint_mp_grid.

        Parameters:
            calc_doc (dict): the structure to converge.
            seed (str): root filename of structure.

        Keyword arguments:
            final (bool): whether this is the final step in a calculation.

        """
        thermo_doc = deepcopy(calc_doc)
        thermo_doc['continuation'] = 'default'
        thermo_doc['task'] = 'thermodynamics'
        thermo_doc['phonon_calculate_dos'] = False

        required = ['phonon_fine_kpoint_mp_spacing']
        forbidden = ['phonon_fine_kpoint_list',
                     'phonon_fine_kpoint_path',
                     'phonon_fine_kpoint_path_spacing']

        self._validate_calc_doc(thermo_doc, required, forbidden)

        return self.scf(thermo_doc, seed, keep=True, intermediate=intermediate)

    def run_convergence_tests(self, calc_doc):
        """ Run kpoint and cutoff_energy convergence tests based on
        options passed to FullRelaxer.

        Parameters:
            calc_doc (dict): the structure to converge.

        Returns:
            bool: True unless every single calculation failed.

        """
        successes = []
        cached_cutoff = calc_doc['cut_off_energy']
        if self.conv_cutoff_bool:
            # run series of singlepoints for various cutoffs
            if self.verbosity > 1:
                print('Running cutoff convergence tests...')
            for cutoff in self.conv_cutoff:
                if self.verbosity > 1:
                    print('{} eV... '.format(cutoff), end='')
                calc_doc.update({'cut_off_energy': cutoff})
                self.paths['completed_dir'] = 'completed_cutoff'
                seed = self.seed + '_' + str(cutoff) + 'eV'
                success = self.scf(calc_doc, seed, keep=False)
                successes.append(success)
        if self.conv_kpt_bool:
            # run series of singlepoints for various cutoffs
            if self.verbosity > 1:
                print('Running cutoff convergence tests...')
            calc_doc['cut_off_energy'] = cached_cutoff
            for kpt in self.conv_kpt:
                if self.verbosity > 1:
                    print('{} 1/A... '.format(kpt), end='')
                calc_doc.update({'kpoints_mp_spacing': kpt})
                self.paths['completed_dir'] = 'completed_kpts'
                seed = self.seed + '_' + str(kpt) + 'A'
                success = self.scf(calc_doc, seed, keep=False)
                successes.append(success)
        return any(successes)

    def parse_executable(self, seed):
        """ Turn executable string into list with arguments to be executed.

        e.g.1:

            | self.executable = 'castep17'
            | seed = 'test'

            | returns
            | ['castep17', 'test']

        e.g.2:

            | self.executable = 'pw6.x -i $seed.in > $seed.out'
            | seed = 'test'

            | returns
            | ['pw6.x', '-i', 'test.in', '>' 'test.out']

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

        return command

    def test_exec(self):
        """ Test if <executable> --version returns a valid string.

        Raises:
            SystemExit: if executable not found.

        """
        try:
            proc = self.run_command('--version', exec_test=True)
        except FileNotFoundError:
            print('Unable to call mpirun/aprun/srun, currently selected: {}'.format(self.mpi_library))
            raise SystemExit('Please check initialistion of FullRelaxer object/CLI args.')

        out, errs = proc.communicate()
        if 'version' not in out.decode('utf-8') and errs is not None:
            err_string = 'Executable {} failed testing. Is it on your PATH?\n'.format(self.executable)
            print_failure(err_string)
            print('stdout: {}'.format(out.decode('utf-8')))
            print('stderr: {}'.format(errs.decode('utf-8')))
            raise SystemExit(err_string)

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
        if sum([self.archer, self.intel, self.slurm]) > 1:
            raise SystemExit('Conflicting command-line arguments for MPI library have been supplied, exiting.')
        elif self.archer:
            return 'archer'
        elif self.intel:
            return 'intel'
        elif self.slurm:
            return 'slurm'
        else:
            return self.detect_mpi(verbosity=self.verbosity)

    @staticmethod
    def detect_mpi(verbosity=0):
        """ Test which mpi library is being used when `mpirun`.

        Returns:
            mpi_library (str): 'intel', 'archer', or 'default'.

        """
        # check first for existence of mpirun command, then aprun if that fails
        try:
            try:
                mpi_version_string = str(sp.check_output('mpirun --version', shell=True))
            except sp.CalledProcessError:
                mpi_version_string = str(sp.check_output('aprun --version', shell=True))
        except Exception as exc:
            print_exc()
            raise exc
        if 'Intel' in mpi_version_string:
            mpi_version = 'intel'
        elif 'aprun' in mpi_version_string:
            mpi_version = 'archer'
        elif 'Open MPI' in mpi_version_string:
            mpi_version = 'default'
        else:
            if verbosity > 3:
                print('Could not detect MPI library, version string was:')
                print(mpi_version_string)
                print('Using default (OpenMPI-style) mpi calls.')
            mpi_version = 'default'

        if verbosity > 2:
            print('Using {} MPI library.'.format(mpi_version))

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
        doc2param(calc_doc, seed, hash_dupe=False)
        doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False)
        if self.verbosity > 2:
            print('Performing memcheck...')
        free_memory = float(virtual_memory().available) / 1024**2
        if self.maxmem is None:
            maxmem = 0.9 * free_memory
        else:
            maxmem = self.maxmem

        # check if cell is totally pathological, as CASTEP dryrun will massively underestimate mem
        if all([angle < 30 for angle in calc_doc['lattice_abc'][1]]):
            if self.debug:
                print('Cell is pathological...')
            return False

        if self.verbosity > 2:
            print('{:10}: {:8.0f} MB'.format('Available', maxmem))
        process = sp.Popen(['nice', '-n', '15', self.executable, '-d', seed])
        process.communicate()

        skip = False
        results, success = castep2dict(seed + '.castep', db=False)
        if not success:
            skip = True

        if 'estimated_mem_MB' not in results:
            skip = True
            if self.debug:
                print('CASTEP dryrun failed, this is probably a bad sign... skipping calculation')
                print(results)

        for _file in glob.glob(seed + '*'):
            if _file.endswith('.res'):
                continue
            else:
                os.remove(_file)

        if self.verbosity > 2:
            if 'estimated_mem_MB' in results:
                print('{:10}: {:8.0f} MB'.format('Estimate', results['estimated_mem_MB']))

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

        return process

    def _catch_castep_errors(self):
        """ Look for CASTEP error files and fallover appropriately.

        TO-DO: better parsing of harmless errors, e.g. LAPACK bugs.

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
                    error_str += ''.join(f.readlines())
                error_str += '\n'
                errors_present = True

        return errors_present, error_str

    def mv_to_bad(self, seed):
        """ Move all files associated with "seed" to bad_castep.

        Parameters:
            seed (str): filename of structure.

        """
        try:
            bad_dir = self.root_folder + '/bad_castep'
            if not os.path.exists(bad_dir):
                os.makedirs(bad_dir, exist_ok=True)
            if self.verbosity > 0:
                print('Something went wrong, moving files to bad_castep')
            seed_files = glob.glob(seed + '.*')
            for _file in seed_files:
                try:
                    shutil.copy(_file, bad_dir)
                    os.remove(_file)
                except Exception:
                    if self.verbosity > 1:
                        print_exc()
            # check root folder for any matching files and remove them
            fname = '{}/{}'.format(self.root_folder, seed)
            for ext in ['.res', '.res.lock', '.castep']:
                if os.path.isfile('{}{}'.format(fname, ext)):
                    os.remove('{}{}'.format(fname, ext))
        except Exception:
            if self.verbosity > 1:
                print_exc()

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
        if not os.path.exists(completed_dir):
            os.makedirs(completed_dir, exist_ok=True)
        if keep:
            seed_files = glob.glob(seed + '.*') + glob.glob(seed + '-out.cell')
            for _file in seed_files:
                shutil.copy(_file, completed_dir)
                os.remove(_file)
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
                    shutil.copy('{}{}'.format(seed, ext), completed_dir)
                    os.remove('{}{}'.format(seed, ext))
                except Exception:
                    if self.verbosity > 1:
                        print_exc()
            # check root folder for any matching files and remove them
            fname = '{}/{}'.format(self.root_folder, seed)
            for ext in file_exts + ['.res.lock']:
                if os.path.isfile('{}{}'.format(fname, ext)):
                    os.remove('{}{}'.format(fname, ext))

    def cp_to_input(self, seed, ext='res', glob_files=False):
        """ Copy initial cell and res to input folder.

        Parameters:
            seed (str): filename of structure.

        Keyword arguments:
            ext (str): file extension for structure.
            glob_files (bool): whether to glob all related seed files.

        """
        input_dir = self.root_folder + '/input'
        if not os.path.exists(input_dir):
            os.makedirs(input_dir, exist_ok=True)
        if glob_files:
            files = glob.glob('{}*'.format(seed))
            for f in files:
                if f.endswith('.lock'):
                    continue
                if not os.path.isfile(f):
                    shutil.copy('{}'.format(f), input_dir)
        else:
            if os.path.isfile('{}.{}'.format(seed, ext)):
                shutil.copy('{}.{}'.format(seed, ext), input_dir)

    def _setup_relaxation_dirs(self):
        """ Set up directories and files for relaxation. """
        if self.verbosity > 0:
            print_notify('Relaxing ' + self.seed)

        if self.compute_dir is not None:
            if not os.path.isdir(self.compute_dir):
                os.makedirs(self.compute_dir)

            # copy pspots and any intermediate calcs to compute_dir
            pspots = glob.glob('*.usp')
            for pspot in pspots:
                shutil.copy(pspot, self.compute_dir)

            # update res file with intermediate calculation if castep file is newer than res
            if os.path.isfile(self.seed + '.castep') and os.path.isfile(self.seed + '.res'):
                shutil.copy(self.seed + '.castep', self.compute_dir)
                if os.path.getmtime(self.seed + '.res') < os.path.getmtime(self.seed + '.castep'):
                    castep_dict, success = castep2dict(self.seed + '.castep', db=False)
                    if success:
                        self.res_dict.update(castep_dict)

            os.chdir(self.compute_dir)

        # copy initial res file to seed
        doc2res(self.res_dict, self.seed, info=False, hash_dupe=False, overwrite=True)
        self.cp_to_input(self.seed)

    def _update_input_files(self):
        """ Update the cell and param files for the next relaxation. """
        calc_doc = self.calc_doc
        # update cell
        if os.path.isfile(self.seed + '.cell'):
            os.remove(self.seed + '.cell')
        if self.kpts_1D:
            if self.verbosity > 0:
                print('Calculating 1D kpt grid...')
            n_kz = ceil(1 / (calc_doc['lattice_abc'][0][2] * self._target_spacing))
            if n_kz % 2 == 1:
                n_kz += 1
            calc_doc['kpoints_mp_grid'] = [1, 1, n_kz]
            if 'kpoints_mp_spacing' in calc_doc:
                del calc_doc['kpoints_mp_spacing']
        doc2cell(calc_doc, self.seed, hash_dupe=False, copy_pspots=False, spin=self.spin)

        # update param
        if self.custom_params:
            if self.verbosity > 0:
                print('Using custom param files...')
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
        for f in glob.glob(seed + '.*'):
            if not (f.endswith('.res') or f.endswith('.castep')):
                os.remove(f)

    def _update_output_files(self, opti_dict):
        """ Copy new data to output files and update
         the results dict.

        Parameter:
            opti_dict (dict): intermediate calculation results.

        """
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
        success = self.res_dict.get('optimised', False)
        if self.output_queue is not None:
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
        if self.verbosity > 2:
            print_notify('Ending process early for seed {}'.format(self.seed))
        process.terminate()
        if self.verbosity > 2:
            print_notify('Ended process early for seed {}'.format(self.seed))
        if self.compute_dir is not None:
            for f in glob.glob('{}.*'.format(self.seed)):
                shutil.copy(f, self.root_folder)
                os.remove(f)

        if os.path.isfile('{}/{}{}'.format(self.root_folder, self.seed, '.res.lock')):
            os.remove('{}/{}{}'.format(self.root_folder, self.seed, '.res.lock'))

    @staticmethod
    def remove_compute_dir_if_finished(compute_dir, debug=True):
        """ Delete the compute directory, provided it contains no
        calculation data.

        Parameters:
            compute_dir (str): path to compute directory.

        Returns:
            bool: True if folder was deleted as no res/castep files
                were found, otherwise False.
        """

        if not os.path.isdir(compute_dir):
            return False

        files = glob.glob(compute_dir + '/*')

        if debug:
            print('Checking {}'.format(compute_dir))
            print('Found {} files:'.format(len(files)))
            print(files)

        for fname in files:
            if fname.endswith('.res') or fname.endswith('.castep'):
                return False

        # remove files in directory, then delete directory
        for fname in files:
            if os.path.isfile(fname):
                try:
                    os.remove(fname)
                    if debug:
                        print('Removing {}'.format(fname))
                except FileNotFoundError:
                    if debug:
                        print('Failed to remove {}'.format(fname))

        if os.path.isdir(compute_dir):
            os.rmdir(compute_dir)
        if debug:
            print('Removed {}'.format(compute_dir))
        return True


class WalltimeError(Exception):
    """ Raise this when you don't want any more jobs to run
    because they're about to exceed the max walltime.

    """
    pass
