# -*- coding: utf-8 -*-

""" This module contains two classes:

* the FullRelaxer class for performing continuously restarted
geometry optimisation and SCF calculations in CASTEP, as well
as the execution of arbitrary programs with mpirun.

* the BatchRun class for running several independent FullRelaxer instances
on a folder of structures, without clashes.

"""

import os
import shutil
import subprocess as sp
import glob
import multiprocessing as mp
import sys
from copy import deepcopy
from traceback import print_exc, format_exception_only
from math import ceil
from collections import defaultdict
from psutil import virtual_memory

from matador.scrapers.castep_scrapers import cell2dict, param2dict
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

        """
        # set defaults and update class with desired values
        prop_defaults = {'paths': None, 'param_dict': None, 'cell_dict': None, 'mode': 'castep', 'executable': 'castep', 'memcheck': False,
                         'rough': 4, 'rough_iter': 2, 'fine_iter': 20, 'spin': False, 'redirect': None, 'reopt': False, 'compute_dir': None,
                         'custom_params': False, 'archer': False, 'maxmem': None, 'killcheck': True, 'kpts_1D': False,
                         'conv_cutoff': False, 'conv_kpt': False, 'debug': False, 'profile': False,
                         'slurm': False, 'intel': False, 'exec_test': True, 'start': True, 'verbosity': 0}
        self.__dict__.update(prop_defaults)
        self.__dict__.update(kwargs)

        if self.profile:
            import cProfile
            import pstats
            from sys import version_info
            from matador.version import __version__
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

        if self.paths is None:
            self.paths = {}
            self.paths['completed_dir'] = 'completed'
        else:
            assert 'completed_dir' in self.paths

        # run through CASTEP specific features
        if self.mode == 'castep':
            self.run_castep(res)
        # otherwise run generic script
        else:
            self.run_generic(res)

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

        Returns:
            bool: True if calculations were successful, False otherwise.
                In the case of convergence tests, this is always True
                unless every calculation fails.

        """
        if self.exec_test:
            self.test_exec()

        if self.kpts_1D:
            if 'kpoints_mp_spacing' not in self.cell_dict:
                raise RuntimeError('kpoints_mp_spacing not found, but kpts_1D requested...')
            self._target_spacing = deepcopy(self.cell_dict['kpoints_mp_spacing'])

        # read in initial structure and skip if failed
        if isinstance(res, str):
            self.res_dict, success = res2dict(res, db=False)
            if not success:
                if self.verbosity >= 1:
                    print(self.res_dict)
                    print_warning('Failed to parse res file {}'.format(res))
                return success
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
                    raise RuntimeError('You forgot your pseudos, you silly goose!')

        # do memcheck, if desired, and only continue if enough memory is free
        if self.memcheck:
            self.enough_memory = self.do_memcheck(calc_doc, self.seed)
        else:
            self.enough_memory = True
        if not self.enough_memory:
            return False

        # run convergence tests
        if any([self.conv_cutoff_bool, self.conv_kpt_bool]):
            return self.run_convergence_tests(calc_doc)

        # perform relaxation
        elif 'GEOMETRY' in calc_doc['task'].upper():
            # set up geom opt paramete
            self.max_iter = calc_doc['geom_max_iter']
            self.num_rough_iter = self.rough
            fine_iter = self.fine_iter
            rough_iter = self.rough_iter
            if 'geom_method' in calc_doc:
                if calc_doc['geom_method'].lower() == 'tpsd' and rough_iter < 3:
                    rough_iter = 3
            num_fine_iter = int(int(self.max_iter)/fine_iter)
            self.geom_max_iter_list = (self.num_rough_iter * [rough_iter])
            self.geom_max_iter_list.extend(num_fine_iter * [fine_iter])
            if not self.geom_max_iter_list:
                raise RuntimeError('Could not divide up relaxation; \
                                    consider increasing geom_max_iter')
            self.calc_doc = calc_doc

            # begin relaxation
            if self.start:
                success = self.relax()
                if self.compute_dir is not None:
                    # always cd back to root folder
                    os.chdir(self.root_folder)
                    self.remove_compute_dir_if_finished(self.compute_dir, debug=self.debug)
                return success

        # run in SCF mode, i.e. just call CASTEP on the seeds
        else:
            return self.scf(calc_doc, self.seed, keep=True)

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
            process = self.run_command(seed)
            process.communicate()
            if process.returncode != 0:
                self.mv_to_bad(seed)
                return False

            self.mv_to_completed(seed, keep=True, completed_dir=self.paths['completed_dir'])
            return True

        except(SystemExit, KeyboardInterrupt):
            print_exc()
            self.mv_to_bad(seed)
            raise SystemExit
        except Exception:
            print_exc()
            self.mv_to_bad(seed)
            return False

    def relax(self, output_queue=None):
        """ Set up a structural relaxation that is restarted intermittently
        in order to re-mesh the kpoint grid. Completed calculations are moved
        to the "completed" folder, and failures to "bad_castep".
            calc_doc (dict): structure to optimise.

        Keyword arguments:
            output_queue (multiprocessing.Queue): optional queue to push
                node and output dict data to upon completion.

        Returns:
            bool: True iff structure was optimised, False otherwise.

        """
        seed = self.seed
        calc_doc = self.calc_doc
        if self.verbosity > 1:
            print_notify('Relaxing ' + self.seed)
        geom_max_iter_list = self.geom_max_iter_list
        if self.debug:
            print(geom_max_iter_list)
            print(self.compute_dir)

        if self.compute_dir is not None:
            if not os.path.isdir(self.compute_dir):
                os.makedirs(self.compute_dir)

            # copy pspots and any intermediate calcs to compute_dir
            pspots = glob.glob('*.usp')
            for pspot in pspots:
                shutil.copy(pspot, self.compute_dir)

            # update res file with intermediate calculation if castep file is newer than res
            if os.path.isfile(self.seed + '.castep') and os.path.isfile(self.seed + '.res'):
                if os.path.getmtime(self.seed + '.res') < os.path.getmtime(self.seed + '.castep'):
                    shutil.copy(self.seed + '.castep', self.compute_dir)
                    castep_dict, success = castep2dict(self.seed + '.castep', db=False)
                    if success:
                        self.res_dict.update(castep_dict)

            os.chdir(self.compute_dir)

        # copy initial res file to seed
        doc2res(self.res_dict, self.seed, info=False, hash_dupe=False, overwrite=True)
        self.cp_to_input(self.seed)

        self.rerun = False
        for ind, num_iter in enumerate(geom_max_iter_list):
            if self.reopt and self.rerun:
                num_iter = 20
                if self.verbosity > 1:
                    print_notify('Performing one last iteration...')
            if self.verbosity > 1:
                if ind == 0:
                    if self.custom_params is not False:
                        print_notify('custom params: {}'.format(self.custom_params))
                    print_notify('Beginning geometry optimisation...')
                elif ind == self.num_rough_iter:
                    print_notify('Beginning fine geometry optimisation...')
            if self.killcheck:
                if os.path.isfile(self.seed + '.kill'):
                    os.remove(self.seed + '.kill')
                    if self.verbosity > 1:
                        print('Found {}.kill, ending job...'.format(self.seed))
                    if output_queue is not None:
                        output_queue.put(self.res_dict)
                        if self.debug:
                            print('wrote failed dict out to output_queue')
                    self.mv_to_bad(seed)
                    return False
            calc_doc['geom_max_iter'] = num_iter
            try:
                # delete any existing files and write new ones
                if os.path.isfile(seed + '.cell'):
                    os.remove(seed + '.cell')
                if self.kpts_1D:
                    if self.verbosity > 1:
                        print('Calculating 1D kpt grid...')
                    n_kz = ceil(1 / (calc_doc['lattice_abc'][0][2] * self._target_spacing))
                    if n_kz % 2 == 1:
                        n_kz += 1
                    calc_doc['kpoints_mp_grid'] = [1, 1, n_kz]
                    if 'kpoints_mp_spacing' in calc_doc:
                        del calc_doc['kpoints_mp_spacing']
                doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False, spin=self.spin)
                if self.custom_params:
                    if self.verbosity > 1:
                        print('Using custom param files...')
                if not self.custom_params:
                    if os.path.isfile(seed + '.param'):
                        os.remove(seed+'.param')
                    doc2param(calc_doc, seed, hash_dupe=False, spin=self.spin)

                # run CASTEP
                process = self.run_command(seed)
                process.communicate()

                # scrape new structure from castep file
                if not os.path.isfile(seed + '.castep'):
                    sys.exit('CASTEP file was not created, please check your executable: {}.'.format(self.executable))
                opti_dict, success = castep2dict(seed + '.castep', db=False, verbosity=self.verbosity)
                if self.debug:
                    print_notify('Intermediate calculation finished')
                    print(opti_dict)
                if not success and isinstance(opti_dict, str):
                    sys.exit('Failed to scrape CASTEP file...')

                # scrub keys that need to be rescraped
                keys_to_remove = ['kpoints_mp_spacing', 'kpoints_mp_grid',
                                  'species_pot', 'sedc_apply', 'sedc_scheme']
                for key in keys_to_remove:
                    if key in opti_dict:
                        del opti_dict[key]

                if self.reopt and self.rerun and not opti_dict['optimised']:
                    self.rerun = False
                if self.reopt and not self.rerun and opti_dict['optimised']:
                    # run once more to get correct symmetry
                    self.rerun = True
                    if os.path.isfile(seed+'.res'):
                        os.remove(seed+'.res')
                    doc2res(opti_dict, seed, hash_dupe=False)
                    if self.compute_dir is not None:
                        if os.path.isfile(seed+'.castep'):
                            shutil.copy(seed+'.castep', self.root_folder)
                        shutil.copy(seed+'.res', self.root_folder)

                elif (not self.reopt or self.rerun) and opti_dict['optimised']:
                    if self.verbosity > 1:
                        print_success('Successfully relaxed ' + seed)
                    # write res and castep file out to completed folder
                    if os.path.isfile(seed+'.res'):
                        os.remove(seed+'.res')
                    doc2res(opti_dict, seed, hash_dupe=False)
                    if self.compute_dir is not None:
                        if os.path.isfile(seed+'.castep'):
                            shutil.copy(seed+'.castep', self.root_folder)
                        shutil.copy(seed+'.res', self.root_folder)

                    # overwrite old data in res_dict with opti structure
                    # so that custom keys in initial res are still accessible
                    self.res_dict.update(opti_dict)
                    if output_queue is not None:
                        output_queue.put(self.res_dict)
                        if self.debug:
                            print('wrote relaxed dict out to output_queue')
                    self.mv_to_completed(seed, completed_dir=self.paths['completed_dir'])
                    # clean up rest of files
                    self.tidy_up(seed)
                    return True
                elif ind == len(geom_max_iter_list) - 1:
                    if self.verbosity > 1:
                        print_warning('Failed to optimise ' + seed)
                    # write final res file to bad_castep
                    if os.path.isfile(seed+'.res'):
                        os.remove(seed+'.res')
                    doc2res(opti_dict, seed, hash_dupe=False)
                    if self.compute_dir is not None:
                        shutil.copy(seed+'.res', self.root_folder)
                        if os.path.isfile(seed+'.castep'):
                            shutil.copy(seed+'.castep', self.root_folder)
                    self.res_dict.update(opti_dict)
                    if output_queue is not None:
                        output_queue.put(self.res_dict)
                        if self.debug:
                            print('wrote failed dict out to output_queue')
                    self.mv_to_bad(seed)
                    return False

                if self._catch_castep_errors(seed, opti_dict, mode='relax'):
                    if output_queue is not None:
                        output_queue.put(self.res_dict)
                        if self.debug:
                            print('wrote failed dict out to output_queue')
                    return False

                # update res file to latest step for restarts
                if os.path.isfile(seed+'.res'):
                    os.remove(seed+'.res')
                doc2res(opti_dict, seed, hash_dupe=False)

                if self.compute_dir is not None:
                    if os.path.isfile(seed+'.castep'):
                        shutil.copy(seed+'.castep', self.root_folder)
                    shutil.copy(seed+'.res', self.root_folder)

                # set atomic_init_spins with value from CASTEP file, if it exists
                if 'mulliken_spins' in calc_doc:
                    calc_doc['atomic_init_spins'] = calc_doc['mulliken_spins']

                # if writing out cell, use it for higher precision lattice_cart
                if calc_doc.get('write_cell_structure'):
                    try:
                        cell_dict, success = cell2dict(seed + '-out.cell', verbosity=self.verbosity, db=False, outcell=True)
                        if success:
                            opti_dict['lattice_cart'] = list(cell_dict['lattice_cart'])
                    except Exception:
                        if self.verbosity > 1:
                            print_exc()

                if self.debug:
                    print_notify('Restarting calculation with current state:')
                    print(calc_doc)
                if self.verbosity >= 2:
                    print(('num_iter: {:3d} | max F: {:5f} eV/A | stress: {: 5f} GPa | ' +
                           'cell volume: {:5f} A^3 | enthalpy per atom {:5f} eV')
                          .format(sum(self.geom_max_iter_list[:ind+1]),
                                  opti_dict['max_force_on_atom'],
                                  opti_dict['pressure'],
                                  opti_dict['cell_volume'],
                                  opti_dict['enthalpy_per_atom']))
                calc_doc.update(opti_dict)

            except(KeyboardInterrupt, FileNotFoundError, SystemExit):
                if self.verbosity > 1:
                    print_exc()
                    print_warning('Received exception, attempting to fail gracefully...')
                etype, evalue, etb = sys.exc_info()
                if self.verbosity > 1:
                    print(format_exception_only(etype, evalue))
                if self.debug:
                    print_exc()
                if self.verbosity > 1:
                    print('Killing CASTEP...')
                process.terminate()
                if self.verbosity > 1:
                    print_warning('Done!')
                    print('Tidying up...')
                self.mv_to_bad(seed)
                self.tidy_up(seed)
                if self.verbosity > 1:
                    print_warning('Done!')
                if output_queue is not None:
                    output_queue.put(self.res_dict)
                    if self.debug:
                        print('wrote failed dict out to output_queue')
                if self.compute_dir is not None:
                    os.chdir(self.root_folder)
                    self.remove_compute_dir_if_finished(self.compute_dir, debug=self.debug)
                return False

            except Exception:
                if self.verbosity > 1:
                    print_exc()
                process.terminate()
                self.mv_to_bad(seed)
                self.tidy_up(seed)
                if output_queue is not None:
                    output_queue.put(self.res_dict)
                    if self.debug:
                        print('wrote ll dict out to output_queue')
                if self.compute_dir is not None:
                    os.chdir(self.root_folder)
                    self.remove_compute_dir_if_finished(self.compute_dir, debug=self.debug)
                return False

    def scf(self, calc_doc, seed, keep=True):
        """ Perform an scf/bandstructure calculation with CASTEP. Files
        from completed runs are moved to "completed" and failed runs to
        "bad_castep".

        Parameters:
            calc_doc (dict): dictionary containing parameters and structure
            seed (str): structure filename

        Keyword arguments:
            keep (bool): whether to keep intermediate files e.g. .bands

        Returns:
            bool: True iff SCF completed successfully, False otherwise.

        """

        try:
            if self.verbosity > 1:
                print_notify('Calculating SCF ' + seed)
            if not self.custom_params:
                doc2param(calc_doc, seed, hash_dupe=False)
            self.cp_to_input(self.seed)

            if 'spectral_task' in calc_doc and calc_doc['spectral_task'] == 'bandstructure':
                if 'spectral_kpoints_path' not in calc_doc and 'spectral_kpoints_list' not in calc_doc:
                    from matador.utils.cell_utils import get_seekpath_kpoint_path, cart2abc
                    if self.verbosity >= 2:
                        from matador.crystal import Crystal
                        print('Old lattice:')
                        print(Crystal(calc_doc))
                    if calc_doc.get('spectral_kpoints_path_spacing') is None:
                        calc_doc['spectral_kpoints_path_spacing'] = 0.02

                    spacing = calc_doc['spectral_kpoints_path_spacing']
                    prim_doc, kpt_path, seekpath_results = get_seekpath_kpoint_path(calc_doc, spacing=spacing, debug=self.debug)
                    if self.verbosity >= 2:
                        print('New lattice:')
                        print(Crystal(calc_doc))
                    calc_doc.update(prim_doc)
                    calc_doc['lattice_abc'] = cart2abc(calc_doc['lattice_cart'])
                    calc_doc['spectral_kpoints_list'] = kpt_path

            doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False, overwrite=True)
            # run CASTEP
            process = self.run_command(seed)
            process.communicate()
            # scrape dict
            opti_dict, success = castep2dict(seed + '.castep', db=False)
            # check for errors
            if self._catch_castep_errors(seed, opti_dict, mode='scf'):
                return False
            if not success:
                return False

            self.mv_to_completed(seed, keep=keep, completed_dir=self.paths['completed_dir'])
            if not keep:
                self.tidy_up(seed)
            return True

        except(SystemExit, KeyboardInterrupt):
            if self.verbosity > 1:
                print_exc()
            self.mv_to_bad(seed)
            if not keep:
                self.tidy_up(seed)
            raise SystemExit
        except Exception:
            if self.verbosity > 1:
                print_exc()
            self.mv_to_bad(seed)
            if not keep:
                self.tidy_up(seed)
            return False

    def run_convergence_tests(self, calc_doc):
        """ Run kpoint and cutoff_energy convergence tests based on
        options passed to FullRelaxer.

        Parameters:
            calc_doc (dict): the structure to converge.

        Returns:
            bool: True unless every single calculation failed.

        """
        successes = []
        self.run_convergence_tests(calc_doc)
        self._cached_cutoff = calc_doc['cut_off_energy']
        if self.conv_cutoff_bool:
            # run series of singlepoints for various cutoffs
            for cutoff in self.conv_cutoff:
                calc_doc.update({'cut_off_energy': cutoff})
                self.paths['completed_dir'] = 'completed_cutoff'
                seed = self.seed + '_' + str(cutoff) + 'eV'
                success = self.scf(calc_doc, seed, keep=False)
                successes.append(success)
        if self.conv_kpt_bool:
            # run series of singlepoints for various cutoffs
            calc_doc['cut_off_energy'] = self._cached_cutoff
            for kpt in self.conv_kpt:
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
            self.redirect_filename = self.redirect.replace('$seed', seed)
        else:
            self.redirect_filename = None

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
            raise RuntimeError('Please check initialistion of FullRelaxer object/CLI args.')

        out, errs = proc.communicate()
        if 'version' not in out.decode('utf-8') and errs is not None:
            err_string = 'Executable {} failed testing. Is it on your PATH?\n'
            err_string += 'Error output: {}'.format(self.executable, errs.decode('utf-8'))
            print_failure(err_string)
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
            raise RuntimeError('Conflicting command-line arguments for MPI library have been supplied, exiting.')
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
            except sp.CalledProcessError as e:
                mpi_version_string = str(sp.check_output('aprun --version', shell=True))
        except Exception as e:
            print_exc()
            raise e
        if 'Intel' in mpi_version_string:
            if verbosity > 1:
                print('Detected Intel MPI library.')
            return 'intel'
        elif 'aprun' in mpi_version_string:
            if verbosity > 2:
                print('Detected ARCHER MPI library.')
            return 'archer'
        else:
            if verbosity > 1:
                print('Could not detect MPI library, version string was:')
                print(mpi_version_string)
                print('Using default (OpenMPI-style) mpi calls.')
            return 'default'

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
        if self.verbosity > 1:
            print('Performing memcheck...')
        free_memory = float(virtual_memory().available) / 1024**2
        if self.maxmem is None:
            maxmem = 0.9*free_memory
        else:
            maxmem = self.maxmem

        # check if cell is totally pathological, as CASTEP dryrun will massively underestimate mem
        if all([angle < 30 for angle in calc_doc['lattice_abc'][1]]):
            if self.debug:
                print('Cell is pathological...')
            return False

        if self.verbosity > 1:
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

        for _file in glob.glob(seed+'*'):
            if _file.endswith('.res'):
                continue
            else:
                os.remove(_file)

        if self.verbosity > 1:
            if 'estimated_mem_MB' in results:
                print('{:10}: {:8.0f} MB'.format('Estimate', results['estimated_mem_MB']))

        if skip or results['estimated_mem_MB'] > maxmem:
            if self.verbosity > 1:
                print('Not enough!')
            return False
        else:
            if self.verbosity > 1:
                print('Enough memory, proceeding...')
            return True

    def run_command(self, seed, exec_test=False):
        """ Calls executable on seed with desired number of cores.

        Parameters:
            seed (str): seedname to pass append to CASTEP command,
                e.g. <seed> or --version.

        Keyword arguments:
            exec_test (bool): run executable in test mode, with output
                piped to stdout.

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
                command = ['ssh', '{}'.format(self.node), 'cd', '{};'.format(cwd),
                           'mpirun', '-n', str(self.ncores)] + command
            else:
                command = ['nice', '-n', '15', 'mpirun', '-n', str(self.ncores)] + command
        else:
            if self.mpi_library == 'archer':
                command = ['aprun', '-n', str(self.ncores*self.nnodes),
                           '-N', str(self.ncores),
                           '-S', '12',
                           '-d', '1'] + command
            elif self.mpi_library == 'slurm':
                command = ['srun', '--exclusive', '-N', str(self.nnodes), '-n', str(self.ncores*self.nnodes)] + command
            elif self.mpi_library == 'intel':
                command = ['mpirun', '-n', str(self.ncores*self.nnodes),
                           '-ppn', str(self.ncores)] + command
            else:
                command = ['mpirun', '-n', str(self.ncores*self.nnodes),
                           '-npernode', str(self.ncores)] + command

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

        if self.redirect_filename is not None:
            redirect_file = open(self.redirect_filename, 'w')
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

    def _catch_castep_errors(self, seed, opti_dict, mode='relax'):
        """ Look for CASTEP error files and fallover appropriately.

        Paramters:
            seed (str): seedname of structure.

        Keyword arguments:
            mode (str): either `relax` or `scf`.

        Returns:
            bool: True if error files were found, otherwise False.

        """
        err_file = '{}*err'.format(seed)
        for globbed in glob.glob(err_file):
            if os.path.isfile(globbed):
                if self.verbosity > 1:
                    print_warning('Failed to optimise {} as CASTEP crashed with error:')
                    with open(globbed, 'r') as f:
                        print(f.readline())

                if mode == 'relax':
                    # write final res file for bad_castep
                    if os.path.isfile(seed+'.res'):
                        os.remove(seed+'.res')
                    self.res_dict.update(opti_dict)
                    doc2res(opti_dict, seed, info=False, hash_dupe=False)

                    if self.compute_dir is not None:
                        if os.path.isfile(seed+'.castep'):
                            shutil.copy(seed+'.castep', self.root_folder)
                        shutil.copy(seed+'.res', self.root_folder)

                self.mv_to_bad(seed)
                return True

        return False

    def mv_to_bad(self, seed):
        """ Move all files associated with "seed" to bad_castep.

        Parameters:
            seed (str): filename of structure.

        """
        try:
            bad_dir = self.root_folder + '/bad_castep'
            if not os.path.exists(bad_dir):
                os.makedirs(bad_dir, exist_ok=True)
            if self.verbosity > 1:
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
            if self.verbosity > 3:
                print(seed_files)
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
            for ext in file_exts:
                try:
                    shutil.copy('{}{}'.format(seed, ext), completed_dir)
                    os.remove('{}{}'.format(seed, ext))
                except Exception:
                    if self.verbosity > 1:
                        print_exc()
            # check root folder for any matching files and remove them
            fname = '{}/{}'.format(self.root_folder, seed)
            for ext in ['.res', '.res.lock', '.castep', '.param', '.cell', '-out.cell']:
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
        try:
            input_dir = self.root_folder + '/input'
            if not os.path.exists(input_dir):
                os.makedirs(input_dir, exist_ok=True)
            if glob_files:
                files = glob.glob('{}*'.format(seed))
                for f in files:
                    if f.endswith('.lock'):
                        continue
                    shutil.copy('{}'.format(f), input_dir)
            else:
                shutil.copy('{}.{}'.format(seed, ext), input_dir)
        except Exception:
            print_exc()

    @staticmethod
    def tidy_up(seed):
        """ Delete all created files before quitting.

        Parameters:

            seed (str): filename for structure.

        """
        for f in glob.glob(seed + '.*'):
            if not (f.endswith('.res') or f.endswith('.castep')):
                os.remove(f)

    @staticmethod
    def remove_compute_dir_if_finished(compute_dir, debug=False):
        """ Delete the compute directory, provided it contains no
        calculation data.

        Parameters:
            compute_dir (str): path to compute directory.

        Returns:
            bool: True if folder was deleted as no res/castep files
                were found, otherwise False.
        """

        if not os.path.isdir(compute_dir):
            if debug:
                print('No directory found called {}, so nothing to do...'.format(compute_dir))
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
            os.remove(fname)
            if debug:
                print('Removing {}'.format(fname))
        os.rmdir(compute_dir)
        if debug:
            print('Removed {}'.format(compute_dir))
        return True


class BatchRun:
    """ A class that implements the running of multiple generic jobs on
    a series of files without collisions with other nodesm using the
    FullRelaxer class. Jobs that have started are listed in jobs.txt,
    failed jobs are moved to bad_castep, completed jobs are moved to
    completed and listed in finished_cleanly.txt.

    Based on run.pl, run2.pl and PyAIRSS class CastepRunner.

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

        """
        # parse args, then co-opt them for passing directly into FullRelaxer
        prop_defaults = {'ncores': None, 'nprocesses': 1, 'nnodes': 1,
                         'executable': 'castep', 'no_reopt': False,
                         'redirect': None, 'debug': False, 'custom_params': False,
                         'verbosity': 0, 'archer': False, 'slurm': False,
                         'intel': False, 'conv_cutoff': False, 'conv_kpt': False,
                         'memcheck': False, 'maxmem': None, 'killcheck': True,
                         'kpts_1D': False, 'spin': False,
                         'rough': 4, 'rough_iter': 2, 'fine_iter': 20,
                         'limit': None, 'profile': False}
        self.args = {}
        self.args.update(prop_defaults)
        self.args.update(kwargs)
        self.debug = self.args.get('debug')
        self.seed = seed
        # if only one seed, check if it is a file, and if so treat
        # this run as a generic run, not a CASTEP cell/param run
        if len(self.seed) == 1:
            if '*' in self.seed[0]:
                self.seed = glob.glob(self.seed[0])
            elif not os.path.isfile(self.seed[0]):
                self.seed = self.seed[0]

        if isinstance(self.seed, str):
            self.mode = 'castep'
        else:
            self.mode = 'generic'

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

        # assign number of cores
        self.all_cores = mp.cpu_count()
        self.slurm_avail_tasks = os.environ.get('SLURM_NTASKS')
        if self.slurm_avail_tasks is not None:
            self.slurm_avail_tasks = int(self.slurm_avail_tasks)

        if self.args.get('ncores') is None:
            if self.slurm_avail_tasks is None:
                self.args['ncores'] = int(self.all_cores / self.nprocesses)
            else:
                self.args['ncores'] = int(self.slurm_avail_tasks / self.nprocesses)
        if self.args['nnodes'] < 1 or self.args['ncores'] < 1 or self.nprocesses < 1:
            sys.exit('Invalid number of cores, nodes or processes.')

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
            join (bool): whether or not to attach to FullRelaxer
                process. Useful for testing.

        """
        from random import sample
        procs = []
        for _ in range(self.nprocesses):
            procs.append(mp.Process(target=self.perform_new_calculations,
                                    args=(sample(self.file_lists['res'],
                                                 len(self.file_lists['res']))
                                          if self.mode == 'castep' else self.seed, )))
        try:
            for proc in procs:
                proc.start()
                if join:
                    proc.join()
        except(KeyboardInterrupt, SystemExit, RuntimeError):
            print_exc()
            for proc in procs:
                proc.terminate()
            sys.exit('Killing running jobs and exiting...')

    def perform_new_calculations(self, res_list):
        """ Perform all calculations that have not already
        failed or finished to completion.

        Parameters:
            res_list (:obj:`list` of :obj:`str`): list of structure filenames.

        """
        job_count = 0
        if isinstance(res_list, str):
            res_list = [res_list]
        for res in res_list:
            locked = os.path.isfile('{}.lock'.format(res))
            listed = self._check_jobs_file(res)
            running = any([listed, locked])
            if not running:

                # check we haven't reached job limit
                if job_count == self.limit:
                    raise SystemExit

                # write lock file
                if not os.path.isfile('{}.lock'.format(res)):
                    with open(res + '.lock', 'a') as job_file:
                        pass
                else:
                    print('Another node wrote this file when I wanted to, skipping...')
                    continue

                # write to jobs file
                with open(self.paths['jobs_fname'], 'a') as job_file:
                    job_file.write(res+'\n')

                # create full relaxer object for creation and running of job
                try:
                    job_count += 1
                    hostname = os.uname()[1]
                    relaxer = FullRelaxer(node=None, res=res,
                                          param_dict=self.param_dict,
                                          cell_dict=self.cell_dict,
                                          mode=self.mode, paths=self.paths, compute_dir=hostname,
                                          **self.args)
                    # if memory check failed, let other nodes have a go
                    if not relaxer.enough_memory:
                        with open(self.paths['memory_fname'], 'a') as job_file:
                            job_file.write(res+'\n')
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
                            job_file.write(res+'\n')
                    else:
                        with open(self.paths['failures_fname'], 'a') as job_file:
                            job_file.write(res+'\n')
                except(KeyboardInterrupt, SystemExit, RuntimeError):
                    print_exc()
                    raise SystemExit

    def generic_setup(self):
        """ Undo things that are set ready for CASTEP jobs... """
        self.cell_dict = None
        self.param_dict = None

    def castep_setup(self):
        """ Set up CASTEP jobs from res files, and $seed.cell/param. """
        # read cell/param files
        exts = ['cell', 'param']
        for ext in exts:
            if not os.path.isfile('{}.{}'.format(self.seed, ext)):
                sys.exit('Failed to find {} file, {}.{}'.format(ext, self.seed, ext))
        self.cell_dict, cell_success = cell2dict(self.seed + '.cell', db=False)
        if not cell_success:
            print(self.cell_dict)
            sys.exit('Failed to parse cell file')
        self.param_dict, param_success = param2dict(self.seed + '.param', db=False)
        if not param_success:
            print(self.param_dict)
            sys.exit('Failed to parse param file')

        # scan directory for files to run
        self.file_lists = defaultdict(list)
        self.file_lists['res'] = [file.name for file in os.scandir() if file.name.endswith('.res')]
        if len(self.file_lists['res']) < 1:
            error = ('run3 in CASTEP mode requires at least 1 res file in folder, found {}'
                     .format(len(self.file_lists['res'])))
            sys.exit(error)

        # do some prelim checks of parameters
        if 'GEOMETRY' in self.param_dict['task'].upper():
            if 'geom_max_iter' not in self.param_dict:
                raise RuntimeError('geom_max_iter is unset, please fix this.')
            elif int(self.param_dict['geom_max_iter']) <= 0:
                raise RuntimeError('geom_max_iter is only {}!'
                                   .format(self.param_dict['geom_max_iter']))

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
                raise RuntimeError('Missing cutoff.conv file')
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
                raise RuntimeError('Missing with conv.kpt file')
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


def reset_job_folder_and_count_remaining(debug=False):
    """ Remove all lock files and clean up jobs.txt
    ready for job restart.

    Note:
        This should be not called by a FullRelaxer instance, in case
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
                else:
                    f.write(line)
            f.truncate()
            flines = f.readlines()
            if debug:
                print('{} jobs remain in jobs.txt'.format(len(flines)))

    return len(res_list)
