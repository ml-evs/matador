# coding: utf-8
""" Contains the FullRelaxer class for continuously
restarted geometry optimisations. Previously part of run3.
"""
# matador modules
from matador.scrapers.castep_scrapers import cell2dict
from matador.scrapers.castep_scrapers import res2dict, castep2dict
from matador.utils.print_utils import print_success, print_warning, print_notify
from matador.export import doc2cell, doc2param, doc2res
# standard library
from os import makedirs, remove, devnull, getcwd
from os.path import isfile, exists
from shutil import copy
from copy import deepcopy
from traceback import print_exc, format_exception_only
from sys import exit, exc_info
from math import ceil
from psutil import virtual_memory
import subprocess as sp
import glob


class FullRelaxer:
    """ Perform full relxation of res input by first doing
    4 rough optimisations with only a few iterations, followed by
    4 larger optimisations with many iterations,
    e.g. 4 lots of 2 then 4 lots of geom_max_iter/4.

    Input:

        res           : str/dict, either filename or input structure dict
        param_dict    : dict, castep parameters
        cell_dict     : dict, castep cell input
        ncores        : int, number of cores for mpirun call
        nnodes        : int, number of nodes for mpirun call
        node          : str, node name to run on

    Keyword arguments:

        executable    : str, name of binary to execute (DEFAULT: castep)
        mode          : str, either 'castep' or 'generic' (DEFAULT: castep)
        custom_params : bool, use custom param file for each structure
        rough         : int, number of small "rough" calculations (DEFAULT: 4)
        spin          : bool, set spins in first calculation (DEFAULT: False)
        conv_cutoffs  : list(float) of cutoffs to use for SCF convergence test
        conv_kpts     : list(float) of kpt spacings to use for SCF convergence test
        kpts_1D       : bool, treat z-direction as special and create kpt_grid [1 1 n_kz].
        paths         : dict, folder names for output sorting
        archer        : bool, use aprun over mpirun
        redirect      : str, file to redirect stdout to (DEFAULT: /dev/null unless debug).
        bnl           : bool, use srun over mpirun
        start         : bool, begin calculation immediately or manually call it
        reopt         : bool, whether to optimise one more time after success (DEFAULT: false)
        memcheck      : bool, perform castep dryrun to estimate memory usage, do not proceed if fails
        maxmem        : int, maximum memory allowed in MB for memcheck
        killcheck     : bool, check for file called $seed.kill during operation, and kill if present

    """
    def __init__(self, res,
                 ncores, nnodes, node, paths=None,
                 param_dict=None, cell_dict=None,
                 mode='castep', executable='castep',
                 rough=None, spin=False, redirect=None,
                 reopt=False, custom_params=False, memcheck=False, maxmem=None, killcheck=True,
                 kpts_1D=False, conv_cutoff=None, conv_kpt=None, archer=False, bnl=False, intel_mpi=False,
                 start=True, verbosity=1, debug=False):
        """ Make the files to run the calculation and handle
        the calling of CASTEP itself.
        """
        self.ncores = ncores
        self.res = res
        self.archer = archer
        self.bnl = bnl
        self.intel_mpi = intel_mpi
        self.nnodes = nnodes
        self.node = node
        self.verbosity = verbosity
        self.memcheck = memcheck
        self.maxmem = maxmem
        self.killcheck = killcheck
        self.executable = executable
        self.mode = mode
        self.redirect = redirect
        self.custom_params = custom_params
        self.reopt = reopt
        self.debug = debug
        self.spin = spin
        self.start = start
        self.kpts_1D = kpts_1D
        self.conv_cutoff_bool = isinstance(conv_cutoff, list)
        self.conv_kpt_bool = isinstance(conv_kpt, list)
        self.conv_cutoff = conv_cutoff
        self.conv_kpt = conv_kpt
        self.enough_memory = True
        self.paths = paths
        if self.paths is None:
            self.paths = {}
            self.paths['completed_dir'] = 'completed'
        else:
            assert 'completed_dir' in self.paths

        self.success = None

        # run through CASTEP specific features
        if self.mode is 'castep':
            if self.kpts_1D:
                assert('kpoints_mp_spacing' in cell_dict)
                self.target_spacing = deepcopy(cell_dict['kpoints_mp_spacing'])

            # read in initial structure and skip if failed
            if isinstance(res, str):
                self.res_dict, success = res2dict(res, db=False)
                if not success:
                    if self.verbosity >= 1:
                        print(self.res_dict)
                        print_warning('Failed to parse res file ' + str(res))
                    self.success = False
            elif isinstance(res, dict):
                self.res_dict = res

            if self.success is None:
                calc_doc = deepcopy(self.res_dict)

                # set seed name
                assert isinstance(calc_doc['source'], list)
                self.seed = calc_doc['source'][0].replace('.res', '')

                # update global doc with cell and param dicts for folder
                calc_doc.update(cell_dict)
                calc_doc.update(param_dict)

                # check for pseudos
                for elem in self.res_dict['stoichiometry']:
                    if '|' not in calc_doc['species_pot'][elem[0]] and\
                            not isfile(calc_doc['species_pot'][elem[0]]):
                        exit('You forgot your pseudos, you silly goose!')

                # run convergence tests
                if any([self.conv_cutoff_bool, self.conv_kpt_bool]):
                    if self.conv_cutoff_bool:
                        # run series of singlepoints for various cutoffs
                        for cutoff in self.conv_cutoff:
                            calc_doc.update({'cut_off_energy': cutoff})
                            self.paths['completed_dir'] = 'completed_cutoff'
                            seed = self.seed + '_' + str(cutoff) + 'eV'
                            self.calc_doc = calc_doc
                            self.success = self.scf(calc_doc, seed, keep=False)
                    if self.conv_kpt_bool:
                        # run series of singlepoints for various cutoffs
                        for kpt in self.conv_kpt:
                            calc_doc.update({'kpoints_mp_spacing': kpt})
                            self.paths['completed_dir'] = 'completed_kpts'
                            seed = self.seed + '_' + str(kpt) + 'A'
                            self.calc_doc = calc_doc
                            self.success = self.scf(calc_doc, seed, keep=False)

                # run simple scf
                elif calc_doc['task'].upper() in ['SPECTRAL', 'SINGLEPOINT']:
                    # batch run density of states
                    self.calc_doc = calc_doc
                    self.success = self.scf(calc_doc, self.seed, keep=True)

                # perform relaxation
                else:
                    # set up geom opt parameters
                    self.max_iter = calc_doc['geom_max_iter']
                    self.num_rough_iter = rough if rough is not None else 4
                    fine_iter = 20
                    rough_iter = 2
                    if 'geom_method' in calc_doc:
                        if calc_doc['geom_method'].lower() == 'tpsd':
                            rough_iter = 3
                    num_fine_iter = int(int(self.max_iter)/fine_iter)
                    self.geom_max_iter_list = (self.num_rough_iter * [rough_iter])
                    self.geom_max_iter_list.extend(num_fine_iter * [fine_iter])
                    self.calc_doc = calc_doc

                    # do memcheck, if desired, and only continue if enough memory is free
                    if self.memcheck:
                        if self.verbosity > 1:
                            print('Trying to perform memcheck...')
                        self.enough_memory = self.do_memcheck(calc_doc, self.seed)
                    else:
                        self.enough_memory = True
                    if self.enough_memory:
                        # begin relaxation
                        if self.start:
                            self.success = self.relax()
        # otherwise run generic script
        else:
            self.seed = res
            if '.' in self.seed:
                self.seed = self.seed.split('.')[-2]
                self.input_ext = self.seed.split('.')[-1]
            else:
                self.input_ext = ''
            assert isinstance(self.seed, str)
            self.run_generic(self.seed)

    def relax(self, output_queue=None):
        """ Set up the calculation to perform 4 sets of two steps,
        then continue with the remainder of steps.

        Optional input:

            output_queue : push node and output dict to a multiprocessing queue (optional).

        Returns:

            True iff structure was optimised, False otherwise.
        """
        seed = self.seed
        calc_doc = self.calc_doc
        if self.verbosity > 1:
            print_notify('Relaxing ' + self.seed)
        geom_max_iter_list = self.geom_max_iter_list
        self.parse_executable(seed)
        # copy initial res file to seed
        if not isinstance(self.res, str):
            doc2res(self.res, self.seed, info=False, hash_dupe=False, overwrite=True)
            self.cp_to_input(self.seed)
        else:
            self.cp_to_input(self.seed)
            doc2res(self.res_dict, self.seed, info=False, hash_dupe=False, overwrite=True)
        self.rerun = False
        for ind, num_iter in enumerate(geom_max_iter_list):
            if self.reopt and self.rerun:
                num_iter = 20
                if self.verbosity > 1:
                    print_notify('Performing one last iteration...')
            if self.verbosity > 1:
                if ind == 0:
                    print_notify('custom params: {}'.format(self.custom_params))
                    print_notify('Beginning rough geometry optimisation...')
                elif ind == self.num_rough_iter:
                    print_notify('Beginning fine geometry optimisation...')
            if self.killcheck:
                if isfile(self.seed + '.kill'):
                    remove(self.seed + '.kill')
                    if self.verbosity > 1:
                        print('Found {}.kill, ending job...'.format(self.seed))
                    if output_queue is not None:
                        output_queue.put(self.res_dict)
                        if self.debug:
                            print('wrote failed dict out to output_queue')
                    self.mv_to_bad(seed)
                    return False
            if ind != 0:
                self.spin = False
            calc_doc['geom_max_iter'] = num_iter
            try:
                # delete any existing files and write new ones
                if isfile(seed + '.cell'):
                    remove(seed + '.cell')
                if self.kpts_1D:
                    if self.verbosity > 1:
                        print('Calculating 1D kpt grid...')
                    n_kz = ceil(1 / (calc_doc['lattice_abc'][0][2] * self.target_spacing))
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
                    if isfile(seed + '.param'):
                        remove(seed+'.param')
                    doc2param(calc_doc, seed, hash_dupe=False)
                # run CASTEP
                process = self.castep(seed)
                process.communicate()
                # scrape new structure from castep file
                opti_dict, success = castep2dict(seed + '.castep', db=False, verbosity=0)
                if self.debug:
                    print_notify('Intermediate calculation finished')
                    print(opti_dict)
                if not success and isinstance(opti_dict, str):
                    if self.verbosity > 1:
                        print_warning('Failed to scrape castep file...')
                    exit()
                try:
                    # delete any k-point and pspot information
                    del opti_dict['kpoints_mp_spacing']
                    del opti_dict['kpoints_mp_grid']
                    del opti_dict['species_pot']
                    del opti_dict['sedc_apply']
                    del opti_dict['sedc_scheme']
                except:
                    pass
                if self.reopt and self.rerun and not opti_dict['optimised']:
                    self.rerun = False
                if self.reopt and not self.rerun and opti_dict['optimised']:
                    # run once more to get correct symmetry
                    self.rerun = True
                    if isfile(seed+'.res'):
                        remove(seed+'.res')
                    doc2res(opti_dict, seed, hash_dupe=False)
                elif (not self.reopt or self.rerun) and opti_dict['optimised']:
                    if self.verbosity > 1:
                        print_success('Successfully relaxed ' + seed)
                    # write res and castep file out to completed folder
                    if isfile(seed+'.res'):
                        remove(seed+'.res')
                    doc2res(opti_dict, seed, hash_dupe=False)
                    self.opti_dict = deepcopy(opti_dict)
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
                    if isfile(seed+'.res'):
                        remove(seed+'.res')
                    doc2res(opti_dict, seed, hash_dupe=False)
                    self.res_dict.update(opti_dict)
                    if output_queue is not None:
                        output_queue.put(self.res_dict)
                        if self.debug:
                            print('wrote failed dict out to output_queue')
                    self.mv_to_bad(seed)
                    return False
                err_file = seed + '*.err'
                for globbed in glob.glob(err_file):
                    if isfile(globbed):
                        if self.verbosity > 1:
                            print_warning('Failed to optimise ' + seed + ' CASTEP crashed.')
                        # write final res file to bad_castep
                        if isfile(seed+'.res'):
                            remove(seed+'.res')
                        self.res_dict.update(opti_dict)
                        if output_queue is not None:
                            output_queue.put(self.res_dict)
                            if self.debug:
                                print('wrote failed dict out to output_queue')
                        doc2res(opti_dict, seed, info=False, hash_dupe=False)
                        self.mv_to_bad(seed)
                        return False

                # update res file to latest step for restarts
                if isfile(seed+'.res'):
                    remove(seed+'.res')
                doc2res(opti_dict, seed, hash_dupe=False)
                # remove atomic_init_spins from calc_doc if there
                if 'atomic_init_spins' in calc_doc:
                    del calc_doc['atomic_init_spins']
                # if writing out cell, use it for higher precision lattice_cart
                if calc_doc.get('write_cell_structure'):
                    cell_dict, success = cell2dict(seed + '-out.cell', db=False, outcell=True)
                    opti_dict['lattice_cart'] = list(cell_dict['lattice_cart'])
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
                etype, evalue, etb = exc_info()
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
                return False
            except:
                if self.verbosity > 1:
                    print_exc()
                process.terminate()
                self.mv_to_bad(seed)
                self.tidy_up(seed)
                if output_queue is not None:
                    output_queue.put(self.res_dict)
                    if self.debug:
                        print('wrote ll dict out to output_queue')
                return False

    def scf(self, calc_doc, seed, keep=True):
        """ Perform only the scf calculation without relaxation.  """
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
                        print('Old lattice:')
                        for i in range(3):
                            print(calc_doc['lattice_cart'][i])
                    if calc_doc.get('spectral_kpoints_path_spacing') is None:
                        spacing = 0.02
                    else:
                        spacing = calc_doc['spectral_kpoints_path_spacing']
                    prim_doc, kpt_path, seekpath_results = get_seekpath_kpoint_path(calc_doc, spacing=spacing, debug=False)
                    if self.verbosity >= 2:
                        print('New lattice:')
                        for i in range(3):
                            print(prim_doc['lattice_cart'][i])
                            print(seekpath_results)
                    calc_doc.update(prim_doc)
                    calc_doc['lattice_abc'] = cart2abc(calc_doc['lattice_cart'])
                    calc_doc['spectral_kpoints_list'] = kpt_path

            doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False, overwrite=True)
            self.parse_executable(seed)
            # run CASTEP
            process = self.castep(seed)
            process.communicate()
            # scrape dict
            opti_dict, success = castep2dict(seed + '.castep', db=False)
            err_file = seed + '.*err'
            for globbed in glob.glob(err_file):
                if isfile(globbed):
                    if self.verbosity > 1:
                        print_warning('Failed to optimise ' + seed + ' CASTEP crashed.')
                    # write final res file to bad_castep
                    self.mv_to_bad(seed)
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
        except:
            if self.verbosity > 1:
                print_exc()
            self.mv_to_bad(seed)
            if not keep:
                self.tidy_up(seed)
            return False

    def run_generic(self, seed):
        """ Run a generic command on the given seed. """
        try:
            self.cp_to_input(seed, ext=self.input_ext, glob_files=True)
            self.parse_executable(seed)
            process = self.castep(seed)
            process.communicate()
            if process.returncode != 0:
                self.mv_to_bad(seed)
                return False
            else:
                self.mv_to_completed(seed, keep=True, completed_dir=self.paths['completed_dir'])
                return True
        except(SystemExit, KeyboardInterrupt):
            print_exc()
            self.mv_to_bad(seed)
            raise SystemExit
        except:
            print_exc()
            self.mv_to_bad(seed)
            return False

    def parse_executable(self, seed):
        """ Turn executable list with arguments into
        command to execute, setting the self.command
        and self.redirect_file variables.

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

        Input:

            seed: str, filename (including extension) to replace $seed with in command.

        """
        if isinstance(self.executable, str):
            executable = self.executable.split()
        command = []
        found_seed = False
        for ind, item in enumerate(executable):
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

        self.command = command

    def do_memcheck(self, calc_doc, seed):
        """ Perform a CASTEP dryrun to estimate memory usage.

        Returns:

            True if the memory estimate is <90% of node RAM,
            otherwise False.

        """
        doc2param(calc_doc, seed, hash_dupe=False)
        doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False)
        if self.debug:
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

        if self.debug:
            print('{:10}: {:8.0f} MB'.format('Available', maxmem))
        process = sp.Popen(['nice', '-n', '15', self.executable, '-d', seed])
        process.communicate()
        results, success = castep2dict(seed + '.castep', db=False)

        skip = False
        if 'estimated_mem_MB' not in results:
            skip = True
            if self.debug:
                print('CASTEP dryrun failed, this is probably a bad sign... skipping calculation')
                print(results)

        for _file in glob.glob(seed+'*'):
            if _file.endswith('.res'):
                continue
            else:
                remove(_file)

        if self.debug:
            if 'estimated_mem_MB' in results:
                print('{:10}: {:8.0f} MB'.format('Estimate', results['estimated_mem_MB']))

        if skip or results['estimated_mem_MB'] > maxmem:
            if self.debug:
                print('Not enough!')
            return False
        else:
            if self.debug:
                print('Enough memory, proceeding...')
            return True

    def castep(self, seed):
        """ Calls executable on desired seed with desired number of cores. """
        if self.nnodes is None or self.nnodes == 1:
            if self.ncores == 1 and self.node is None:
                command = ['nice', '-n', '15'] + self.command
            elif self.archer:
                command = ['aprun', '-n', str(self.ncores)] + self.command
            elif self.bnl:
                command = ['srun', '--exclusive', '-N', '1', '-n', str(self.ncores)] + self.command
            elif self.intel_mpi:
                command = ['mpirun', '-n', str(self.ncores), '-ppn', str(self.ncores)] + self.command
            elif self.node is not None:
                cwd = getcwd()
                command = ['ssh', '{}'.format(self.node), 'cd', '{};'.format(cwd),
                           'mpirun', '-n', str(self.ncores)] + self.command
            else:
                command = ['nice', '-n', '15', 'mpirun', '-n', str(self.ncores)] + self.command
        else:
            if self.archer:
                command = ['aprun', '-n', str(self.ncores*self.nnodes),
                           '-N', str(self.ncores),
                           '-S', '12',
                           '-d', '1'] + self.command
            elif self.bnl:
                command = ['srun', '--exclusive', '-N', str(self.nnodes), '-n', str(self.ncores*self.nnodes)] + self.command
            elif self.intel_mpi:
                command = ['mpirun', '-n', str(self.ncores*self.nnodes),
                           '-ppn', str(self.ncores)] + self.command
            else:
                command = ['mpirun', '-n', str(self.ncores*self.nnodes),
                           '-npernode', str(self.ncores)] + self.command

        if self.debug:
            if self.redirect_filename is not None:
                redirect_file = open(self.redirect_filename, 'w')
                process = sp.Popen(command, shell=False, stdout=redirect_file)
                redirect_file.close()
            else:
                process = sp.Popen(command, shell=False)
        else:
            dev_null = open(devnull, 'w')
            if self.redirect_filename is not None:
                redirect_file = open(self.redirect_filename, 'w')
                process = sp.Popen(command, shell=False, stdout=redirect_file, stderr=dev_null)
                redirect_file.close()
            else:
                process = sp.Popen(command, shell=False, stdout=dev_null, stderr=dev_null)
            dev_null.close()

        return process

    def mv_to_bad(self, seed):
        """ Move all associated files to bad_castep. """
        try:
            if not exists('bad_castep'):
                makedirs('bad_castep', exist_ok=True)
            if self.verbosity > 1:
                print('Something went wrong, moving files to bad_castep')
            seed_files = glob.glob(seed + '.*')
            for _file in seed_files:
                try:
                    copy(_file, 'bad_castep')
                    remove(_file)
                except:
                    if self.verbosity > 1:
                        print_exc()
                    pass
        except:
            if self.verbosity > 1:
                print_exc()
            pass
        return

    def mv_to_completed(self, seed, completed_dir='completed', keep=False):
        """ Move all associated files to completed. """
        if not exists(completed_dir):
            makedirs(completed_dir, exist_ok=True)
        if keep:
            seed_files = glob.glob(seed + '.*') + glob.glob(seed + '-out.cell')
            if self.debug > 3:
                print(seed_files)
            for _file in seed_files:
                copy(_file, completed_dir)
                remove(_file)
        else:
            file_exts = ['.castep']
            if self.kpts_1D:
                file_exts.append('.param')
            if not self.conv_kpt_bool and not self.conv_cutoff_bool:
                file_exts.append('.res')
            if self.calc_doc.get('write_cell_structure'):
                file_exts.append('-out.cell')
            for ext in file_exts:
                try:
                    copy('{}{}'.format(seed, ext), completed_dir)
                    remove('{}{}'.format(seed, ext))
                except:
                    if self.verbosity > 1:
                        print_exc()
                    pass
        return

    def cp_to_input(self, seed, ext='res', glob_files=False):
        """ Copy initial cell and res to input folder. """
        try:
            if not exists('input'):
                makedirs('input', exist_ok=True)
            if glob_files:
                files = glob.glob('{}*'.format(seed))
                for f in files:
                    if f.endswith('.lock'):
                        continue
                    copy('{}'.format(f), 'input')
            else:
                copy('{}.{}'.format(seed, ext), 'input')
        except:
            print_exc()
            if self.verbosity > 1:
                print_exc()
            pass
        return

    def tidy_up(self, seed):
        """ Delete all run3 created files before quitting. """
        for f in glob.glob(seed + '.*'):
            if not (f.endswith('.res') or f.endswith('.castep')):
                remove(f)
        return
