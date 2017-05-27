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
from os import makedirs, remove, system, devnull, getcwd, getpid
from os.path import isfile, exists
from shutil import copy2
from copy import deepcopy
from traceback import print_exc, format_exception_only
from sys import exit, exc_info
import sys
import subprocess as sp
import glob


class FullRelaxer:
    """ Perform full relxation of res input by first doing
    4 rough optimisations with only a few iterations, followed by
    4 larger optimisations with many iterations,
    e.g. 4 lots of 2 then 4 lots of geom_max_iter/4.

    Input:

        ncores      : number of cores for mpirun call
        nnodes      : number of nodes for mpirun call (DEPCRECATED)
        node        : node name to run on
        res         : either filename or input structure dict
        param_dict  : dict of castep parameters
        cell_dict   : dict of castep cell input
        executable  : name of binary to execute (DEFAULT: castep)
        rough       : number of small "rough" calculations (DEFAULT: 4)
        spin        : set spins in first calculation (DEFAULT: False)
        conv_cutoff : read cutoffs from cutoff.conv and run them all
        conv_kpt    : read kpt spacings kpt.conv and run them all
        archer      : use aprun over mpirun
        bnl         : use srun over mpirun
        start       : begin calculation immediately or manually call it
        redirect    : redirect all output to pid.file

    """
    def __init__(self, ncores, nnodes, node, res, param_dict, cell_dict,
                 executable='castep', rough=None, spin=False,
                 reopt=True,
                 conv_cutoff=None, conv_kpt=None,
                 archer=False, bnl=False,
                 start=True, redirect=False,
                 verbosity=0, debug=False):
        """ Make the files to run the calculation and handle
        the calling of CASTEP itself.
        """
        self.ncores = ncores
        self.res = res
        self.archer = archer
        self.bnl = bnl
        self.nnodes = nnodes
        self.node = node
        self.verbosity = verbosity
        self.executable = executable
        self.reopt = reopt
        self.debug = debug
        self.spin = spin
        self.start = start
        self.conv_cutoff_bool = True if conv_cutoff is not None else False
        self.conv_kpt_bool = True if conv_kpt is not None else False
        if self.conv_cutoff_bool:
            self.conv_cutoff = conv_cutoff
        if self.conv_kpt_bool:
            self.conv_kpt = conv_kpt
        self.success = None
        if redirect:
            self.redirect = True
            sys.stdout = open(str(getpid()) + '.out', 'w')

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

            if self.conv_cutoff_bool:
                # run series of singlepoints for various cutoffs
                for cutoff in self.conv_cutoff:
                    calc_doc.update({'cut_off_energy': cutoff})
                    seed = self.seed + '_' + str(cutoff) + 'eV'
                    self.success = self.scf(calc_doc, seed, keep=False)

            elif self.conv_kpt_bool:
                # run series of singlepoints for various cutoffs
                for kpt in self.conv_kpt:
                    calc_doc.update({'kpoints_mp_spacing': kpt})
                    seed = self.seed + '_' + str(kpt) + 'A'
                    self.success = self.scf(calc_doc, seed, keep=False)

            elif calc_doc['task'].upper() in ['SPECTRAL', 'SINGLEPOINT']:
                # batch run density of states
                self.success = self.scf(calc_doc, self.seed, keep=True)
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

                # begin relaxation
                if self.start:
                    self.success = self.relax()

    def relax(self, output_queue=None):
        """ Set up the calculation to perform 4 sets of two steps,
        then continue with the remainder of steps.

        Optional input:

        output_queue : push node and output dict to a multiprocessing queue (optional).

        Returns:

        True iff structure was optimised, false otherwise.
        """
        seed = self.seed
        calc_doc = self.calc_doc
        if self.verbosity >= 1:
            print_notify('Relaxing ' + self.seed)
        geom_max_iter_list = self.geom_max_iter_list
        # copy initial res file to seed
        self.cp_to_input(self.seed)
        if not isinstance(self.res, str):
            doc2res(self.res, self.seed, info=False, hash_dupe=False)

        self.rerun = False
        for ind, num_iter in enumerate(geom_max_iter_list):
            if self.reopt and self.rerun:
                num_iter = 20
                if self.verbosity >= 1:
                    print_notify('Performing one last iteration...')
            if self.verbosity >= 1:
                if ind == 0:
                        print_notify('Beginning rough geometry optimisation...')
                elif ind == self.num_rough_iter:
                    print_notify('Beginning fine geometry optimisation...')
            if ind != 0:
                self.spin = False
            calc_doc['geom_max_iter'] = num_iter
            try:
                # delete any existing files
                if isfile(seed + '.param'):
                    remove(seed+'.param')
                if isfile(seed + '.cell'):
                    remove(seed+'.cell')
                # write new param and cell
                doc2param(calc_doc, seed, hash_dupe=False)
                doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False, spin=self.spin)
                # run CASTEP
                process = self.castep(seed)
                process.communicate()
                # scrape new structure from castep file
                opti_dict, success = castep2dict(seed + '.castep', db=False)
                if self.debug:
                    print_notify('Intermediate calculation finished')
                    print(opti_dict)
                if not success and isinstance(opti_dict, str):
                    if self.verbosity >= 1:
                        print_warning('Failed to scrape castep file...')
                    exit()
                try:
                    # delete any k-point and pspot information
                    del opti_dict['kpoints_mp_spacing']
                    del opti_dict['kpoints_mp_grid']
                    del opti_dict['species_pot']
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
                    if self.verbosity >= 1:
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
                    self.mv_to_completed(seed)
                    if calc_doc.get('write_cell_structure'):
                        system('mv ' + seed + '-out.cell' + ' completed/' + seed + '-out.cell')
                    # clean up rest of files
                    self.tidy_up(seed)
                    return True
                elif ind == len(geom_max_iter_list) - 1:
                    if self.verbosity >= 1:
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
                        if self.verbosity >= 1:
                            print_warning('Failed to optimise ' + seed + ' CASTEP crashed.')
                        # write final res file to bad_castep
                        if isfile(seed+'.res'):
                            remove(seed+'.res')
                        self.res_dict.update(opti_dict)
                        if output_queue is not None:
                            output_queue.put(self.res_dict)
                            if self.debug:
                                print('wrote failed dict out to output_queue')
                        doc2res(opti_dict, seed, hash_dupe=False)
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
                    print(('num_iter: {:3d} | max F: {:5f} eV/A | stress: {: 5f} GPa | '
                          + 'cell volume: {:5f} A^3 | enthalpy per atom {:5f} eV')
                          .format(sum(self.geom_max_iter_list[:ind+1]),
                                  opti_dict['max_force_on_atom'],
                                  opti_dict['pressure'],
                                  opti_dict['cell_volume'],
                                  opti_dict['enthalpy_per_atom']))
                calc_doc.update(opti_dict)

            except(KeyboardInterrupt, SystemExit):
                if self.verbosity >= 1:
                    print_exc()
                    print_warning('Received exception, attempting to fail gracefully...')
                etype, evalue, etb = exc_info()
                if self.verbosity >= 1:
                    print(format_exception_only(etype, evalue))
                if self.debug:
                    print_exc()
                if self.verbosity >= 1:
                    print('Killing CASTEP...')
                process.terminate()
                if self.verbosity >= 1:
                    print_warning('Done!')
                    print('Tidying up...', end=' ')
                self.mv_to_bad(seed)
                self.tidy_up(seed)
                if self.verbosity >= 1:
                    print_warning('Done!')
                return False
            except:
                if self.verbosity >= 1:
                    print_exc()
                process.terminate()
                self.mv_to_bad(seed)
                self.tidy_up(seed)
                return False

    def scf(self, calc_doc, seed, keep=True):
        """ Perform only the scf calculation without relaxation.  """
        try:
            if self.verbosity >= 1:
                print_notify('Calculating SCF ' + self.seed)
            doc2param(calc_doc, seed, hash_dupe=False)
            doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False)
            # run CASTEP
            process = self.castep(seed)
            process.communicate()
            # scrape dict
            opti_dict, success = castep2dict(seed + '.castep', db=False)
            err_file = seed + '.*err'
            for globbed in glob.glob(err_file):
                if isfile(globbed):
                    if self.verbosity >= 1:
                        print_warning('Failed to optimise ' + seed + ' CASTEP crashed.')
                    # write final res file to bad_castep
                    self.mv_to_bad(seed)
                    return False
            self.mv_to_completed(seed, keep)
            if not keep:
                self.tidy_up(seed)
            return True
        except(SystemExit, KeyboardInterrupt):
            if self.verbosity >= 1:
                print_exc()
            self.mv_to_bad(seed)
            if not keep:
                self.tidy_up(seed)
            raise SystemExit
        except:
            if self.verbosity >= 1:
                print_exc()
            if not keep:
                self.mv_to_bad(seed)
            self.tidy_up(seed)
            return False

    def castep(self, seed):
        """ Calls CASTEP on desired seed with desired number of cores.
        """
        if self.nnodes is None or self.nnodes == 1:
            if self.ncores == 1 and self.node is None:
                process = sp.Popen(['nice', '-n', '15', self.executable, seed])
            elif self.archer:
                process = sp.Popen(['aprun', '-n', str(self.ncores),
                                    self.executable, seed])
            elif self.bnl:
                command = ['srun', '-n', str(self.ncores), self.executable, seed]
                if self.debug:
                    print(command)
                process = sp.Popen(command)
            elif self.node is not None:
                cwd = getcwd()
                if self.debug:
                    process = sp.Popen(['ssh', '{}'.format(self.node),
                                        'cd', '{};'.format(cwd),
                                        'mpirun', '-n', str(self.ncores),
                                        self.executable, seed],
                                       shell=False)
                else:
                    dev_null = open(devnull, 'w')
                    process = sp.Popen(['ssh', '{}'.format(self.node),
                                        'cd', '{};'.format(cwd),
                                        'mpirun', '-n', str(self.ncores),
                                        self.executable, seed],
                                       shell=False, stdout=dev_null, stderr=dev_null)
                    dev_null.close()
            else:
                if self.debug:
                    process = sp.Popen(['nice', '-n', '15', 'mpirun', '-n', str(self.ncores),
                                        self.executable, seed])
                else:
                    dev_null = open(devnull, 'w')
                    process = sp.Popen(['nice', '-n', '15', 'mpirun', '-n', str(self.ncores),
                                        self.executable, seed], stdout=dev_null, stderr=dev_null)
                    dev_null.close()
        else:
            raise NotImplementedError
            if self.archer:
                command = ['aprun', '-n', str(self.ncores*self.nnodes),
                           '-N', str(self.ncores),
                           '-S', '12',
                           '-d', '1',
                           self.executable, seed]
                if self.verbosity >= 1:
                    print(command)
                process = sp.Popen(command)
            else:
                if self.verbosity >= 1:
                    print(['mpirun', '-n', str(self.ncores*self.nnodes),
                           '-ppn', str(self.ncores),
                           self.executable, seed])
                process = sp.Popen(['mpirun', '-n', str(self.ncores*self.nnodes),
                                    '-ppn', str(self.ncores),
                                    self.executable, seed])
        return process

    def mv_to_bad(self, seed):
        """ Move all associated files to bad_castep. """
        try:
            if not exists('bad_castep'):
                makedirs('bad_castep', exist_ok=True)
            if self.verbosity >= 1:
                print('Something went wrong, moving files to bad_castep')
            system('mv ' + seed + '* bad_castep')
        except:
            if self.verbosity > 0:
                print_exc()
            pass
        return

    def mv_to_completed(self, seed, keep=False):
        """ Move all associated files to completed. """
        if not exists('completed'):
            makedirs('completed', exist_ok=True)
        if keep:
            system('mv ' + seed + '*' + ' completed')
        else:
            system('mv ' + seed + '.castep' + ' completed/' + seed + '.castep')
            system('mv ' + seed + '.res' + ' completed/' + seed + '.res')
        return

    def cp_to_input(self, seed):
        """ Copy initial cell and res to input folder. """
        try:
            if not exists('input'):
                makedirs('input', exist_ok=True)
            _ = copy2('{}.res'.format(seed), 'input')
        except:
            if self.verbosity > 0:
                print_exc()
            pass
        return

    def tidy_up(self, seed):
        """ Delete all run3 created files before quitting. """
        for f in glob.glob(seed + '.*'):
            if not (f.endswith('.res') or f.endswith('.castep')):
                remove(f)
        return
