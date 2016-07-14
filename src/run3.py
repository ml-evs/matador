#!/usr/bin/python
# coding :utf-8
""" Run many CASTEP calculations from multiple res,
a single cell and a single param file.

Jobs in progress are listed in jobs.txt, failed
jobs are moved to bad_castep, completed jobs are
moved to completed and listed in finished_cleanly.txt.

Based on run.pl, run2.pl and PyAIRSS class CastepRunner.

Matthew Evans 2016
"""

from __future__ import print_function
from os import walk, makedirs, remove, system
from os.path import isfile, exists
from collections import defaultdict
from scrapers.castep_scrapers import cell2dict, param2dict
from scrapers.castep_scrapers import res2dict, castep2dict
from print_utils import print_failure, print_success, print_warning, print_notify
from export import doc2cell, doc2param, doc2res
from traceback import print_exc
import argparse
import multiprocessing as mp
import subprocess as sp
import glob


class BatchRun:
    """ A class that implements the running of multiple CASTEP jobs from a series
    of .res files and single cell and param files.
    """
    def __init__(self, *args, **kwargs):
        """ Check directory has valid contents and prepare log files and directories if
        not already prepared, then begin running calculations.
        """
        # analyse parallel allocation
        self.args = kwargs
        self.debug = self.args.get('debug')
        self.all_cores = mp.cpu_count()
        self.seed = self.args.get('seed')
        self.limit = self.args.get('limit')
        print(self.limit)
        valid = True
        if self.args.get('nprocesses') is not None:
            self.nprocesses = self.args['nprocesses']
        else:
            self.nprocesses = 1
        if self.args.get('ncores') is not None:
            self.ncores = self.args['ncores']
        else:
            self.ncores = self.all_cores / self.nprocesses
        if self.args.get('nnodes') is not None:
            self.nnodes = self.args['nnodes']
            print_warning('Attempting to run over multiple nodes, please ensure \
                           that you are using Intel MPI or this will produce \
                           unexpected behaviour!')
        else:
            self.nnodes = 1 
        try:
            assert self.nnodes >= 1
            assert self.ncores >= 1
            assert self.nprocess >= 1
        except(AssertionError):
            print_failure('Invalid number of cores, nodes or processes.')
            exit()
        if self.ncores*self.nprocesses > self.all_cores:
            print_warning('Requested more cores (' + str(self.ncores*self.nprocesses) +
                          ') than available (' + str(self.all_cores) + ').' +
                          '(Maybe you are using a queueing system).')
        # scan directory for files to run
        self.file_lists = defaultdict(list)
        for root, dirs, files in walk('.'):
            for file in files:
                if root == '.':
                    if file.endswith('.res'):
                        self.file_lists['res'].append(file)
        if len(self.file_lists['res']) < 1:
            valid = False
            print_failure('run3 requires at least 1 res file in folder, found' +
                          str(len(self.file_lists['res'])))
        self.cell_dict, cell_success = cell2dict(self.seed + '.cell', db=False)
        if not cell_success:
            valid = False
            print_failure('Failed to parse cell file')
        self.param_dict, param_success = param2dict(self.seed + '.param', db=False)
        if not param_success:
            valid = False
            print_failure('Failed to parse cell file')
        if int(self.param_dict['geom_max_iter']) < 20:
            valid = False
            print_failure('geom_max_iter is only ' +
                          str(self.param_dict['geom_max_iter']) + '... quitting.')
        if self.args.get('conv_cutoff'):
            try:
                if isfile('cutoff.conv'):
                    with open('cutoff.conv', 'r') as f:
                        flines = f.readlines()
                        self.cutoffs = []
                        for line in flines:
                            self.cutoffs.append(int(line))
                else:
                    raise RuntimeError
            except:
                valid = False
                print_failure('Error with cutoff.conv file.')
        else:
            self.cutoffs = None
        if not valid:
            exit('Exiting...')
        # delete source from cell and param
        del self.cell_dict['source']
        del self.param_dict['source']
        # prepare folders and text files
        self.paths = dict()
        self.paths['completed_dir'] = 'completed'
        self.paths['failed_dir'] = 'bad_castep'
        self.paths['jobs_fname'] = 'jobs.txt'
        self.paths['completed_fname'] = 'finished_cleanly.txt'
        if not isfile(self.paths['jobs_fname']):
            with open(self.paths['jobs_fname'], 'a'):
                pass
        if not isfile(self.paths['completed_fname']):
            with open(self.paths['completed_fname'], 'a'):
                pass

    def spawn(self):
        # begin running calculations that have not already started
        # initialise pool of processes
        procs = []
        for ind in range(self.nprocesses):
            procs.append(mp.Process(target=self.perform_new_calculations,
                         args=(self.file_lists['res'], self.paths)))
        try:
            for proc in procs:
                proc.start()
        except(KeyboardInterrupt, SystemExit, RuntimeError):
            for proc in procs:
                proc.terminate()
            exit('Killing running jobs and exiting...')

    def perform_new_calculations(self, res_list, paths):
        """ Perform all calculations that have not already
        failed or finished to completion. """
        job_count = 0
        for res in res_list:
            running = False
            with open(paths['jobs_fname'], 'rw') as job_file:
                flines = job_file.readlines()
                for line in flines:
                    if res in line:
                        running = True
                        break
            if not running:
                if self.limit is not None:
                    if job_count == self.limit:
                        print(self.limit)
                        raise SystemExit
                with open(paths['jobs_fname'], 'a') as job_file:
                    job_file.write(res+'\n')
                # create full relaxer object for creation and running of job
                try:
                    FullRelaxer(paths=self.paths,
                                ncores=self.ncores,
                                nnodes=self.nnodes,
                                res=res,
                                param_dict=self.param_dict,
                                cell_dict=self.cell_dict,
                                debug=self.debug,
                                conv_cutoff=self.cutoffs)
                    with open(paths['completed_fname'], 'a') as job_file:
                        job_file.write(res+'\n')
                    job_count += 1
                except(KeyboardInterrupt, SystemExit, RuntimeError):
                    raise SystemExit
        return


class FullRelaxer:
    """ Perform full relxation of res input by first doing
    4 rough optimisations with only a few iterations, followed by
    4 larger optimisations with many iterations,
    e.g. 4 lots of 2 then 4 lots of geom_max_iter/4.
    """
    def __init__(self, paths, ncores, nnodes, res, param_dict, cell_dict,
                 debug=False, conv_cutoff=None):
        """ Make the files to run the calculation and handle
        the calling of CASTEP itself.
        """
        self.paths = paths
        self.ncores = ncores
        self.nnodes = nnodes
        self.executable = 'castep'
        self.debug = debug
        self.conv_cutoff_bool = True if conv_cutoff is not None else False
        if self.conv_cutoff_bool:
            self.conv_cutoff = conv_cutoff
        # read in initial structure
        res_dict, success = res2dict(res, db=False)
        calc_doc = res_dict
        # set seed name
        self.seed = calc_doc['source'][0].replace('.res', '')
        print_notify('Relaxing ' + self.seed)
        # update global doc with cell and param dicts for folder
        calc_doc.update(cell_dict)
        calc_doc.update(param_dict)
        calc_doc['task'] = 'geometryoptimization'
        # set up geom opt parameters
        self.max_iter = calc_doc['geom_max_iter']
        self.num_rough_iter = 4
        fine_iter = 20
        rough_iter = 2
        num_fine_iter = int(self.max_iter)/fine_iter
        self.geom_max_iter_list = (self.num_rough_iter * [rough_iter])
        self.geom_max_iter_list.extend(num_fine_iter * [fine_iter])
        if self.conv_cutoff_bool:
            for cutoff in self.conv_cutoff:
                calc_doc.update({'cut_off_energy': cutoff})
                seed = self.seed + '_' + str(cutoff) + 'eV'
                self.success = self.relax(calc_doc, seed)
        else:
            self.success = self.relax(calc_doc, self.seed)
        if not success:
            self.mv_to_bad(self.seed)

    def relax(self, calc_doc, seed):
        """ Set up the calculation to perform 4 sets of two steps,
        then continue with the remainder of steps.
        """
        geom_max_iter_list = self.geom_max_iter_list
        # relax structure
        print(geom_max_iter_list)
        # copy initial res file to seed
        self.cp_to_input(self.seed)
        for ind, num_iter in enumerate(geom_max_iter_list):
            if ind == 0:
                print_notify('Beginning rough geometry optimisation...')
            elif ind == self.num_rough_iter:
                print_notify('Beginning fine geometry optimisation...')
            calc_doc['geom_max_iter'] = num_iter
            try:
                # delete any existing files
                if isfile(seed + '.param'):
                    remove(seed+'.param')
                if isfile(seed + '.cell'):
                    remove(seed+'.cell')
                # write new param and cell
                doc2param(calc_doc, seed, hash_dupe=False)
                doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False)
                # run CASTEP
                process = self.castep(seed)
                process.communicate()
                # scrape new structure from castep file
                opti_dict, success = castep2dict(seed + '.castep', db=False)
                try:
                    # delete any k-point information
                    del opti_dict['kpoint_mp_spacing']
                    del opti_dict['kpoint_mp_grid']
                except:
                    pass
                if not success and opti_dict == '':
                    print_warning('Failed to scrape castep file...')
                    return False
                if opti_dict['optimised'] == True:
                    if not exists('completed'):
                        makedirs('completed')
                    print_success('Successfully relaxed ' + seed)

                    # write res and castep file out to completed folder
                    doc2res(opti_dict, 'completed/' + seed, hash_dupe=False)
                    system('mv ' + seed + '.castep' + ' completed/' + seed + '.castep')
                    system('mv ' + seed + '.param' + ' completed/' + seed + '.param')
                    system('mv ' + seed + '.cell' + ' completed/' + seed + '.cell')
                    if calc_doc.get('write_cell_structure') == 'true':
                        system('mv ' + seed + '-out.cell' + ' completed/' + seed + '-out.cell')
                    # clean up rest of files
                    self.tidy_up(seed)
                    return True
                elif ind == len(geom_max_iter_list) - 1:
                    print_warning('Failed to optimise ' + seed)
                    self.mv_to_bad(seed)
                    # write final res file to bad_castep
                    doc2res(opti_dict, 'bad_castep/' + seed, hash_dupe=False)
                else:
                    err_file = seed + '*.err'
                    for globbed in glob.glob(err_file):
                        if isfile(globbed):
                            print_warning('Failed to optimise ' + seed + ' CASTEP crashed.')
                            # write final res file to bad_castep
                            doc2res(opti_dict, 'bad_castep/' + seed, hash_dupe=False)
                            self.mv_to_bad(seed)
                            return False
                    # update res file to latest step for restarts
                    doc2res(opti_dict, seed, hash_dupe=False)
                    calc_doc.update(opti_dict)
            except(SystemExit, KeyboardInterrupt):
                self.mv_to_bad(seed)
                self.tidy_up(seed)
                raise SystemExit
            except:
                print_exc()
                self.mv_to_bad(seed)
                self.tidy_up(seed)
                return False

    def castep(self, seed):
        """ Calls CASTEP on desired seed with desired number of cores.
        """
        if self.nnodes == 1:
            if self.ncores == 1:
                process = sp.Popen(['nice', '-n', '15', self.executable, seed])
            else:
                process = sp.Popen(['nice', '-n', '15', 'mpirun', '-n', str(self.ncores),
                                    self.executable, seed])
        elif self.nnodes > 1:
            process = sp.Popen(['mpirun', '-n', str(self.ncores*self.nnodes), '-ppn', str(self.ncores)])
        return process

    def mv_to_bad(self, seed):
        """ Move all associated files to bad_castep. """
        if not exists('bad_castep'):
            makedirs('bad_castep')
        print('Something went wrong, moving files to bad_castep')
        system('mv ' + seed + '* bad_castep')
        return

    def cp_to_input(self, seed):
        """ Copy initial cell and res to input folder. """
        if not exists('input'):
            makedirs('input')
        system('cp ' + seed + '.res input')
        return

    def tidy_up(self, seed):
        """ Delete all run3 created files before quitting. """
        system('rm ' + seed + '*')
        return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='run3',
        description='Run multiple CASTEP geometry optmizations from a series of .res \
                     files and single cell and param files. The geometry optimization will \
                     be split into four chunks of 2 iteratiosn, followed by chunks of 20 iterations, \
                     until geom_max_iter is reached in the param file. \
                     Successful runs will be moved to completed, crashes/failures will go to bad_castep \
                     and initial res files will go into input. Running jobs will be listed in jobs.txt \
                     and those that completed cleanly will be listed in finished_cleanly.txt.',
        epilog='Written by Matthew Evans (2016), based primarily on run.pl and run2.pl \
                by Chris Pickard and Andrew Morris and PyAIRSS CastepRunner by Jamie Wynn.')
    parser.add_argument('seed', type=str,
                        help='cell and param seed to use as template for calculations')
    parser.add_argument('-nc', '--ncores', type=int,
                        help='number of cores CASTEP per job [DEFAULT=cpu_count/nprocesses]')
    parser.add_argument('-np', '--nprocesses', type=int,
                        help='number of concurrent calculations, i.e. number \
                              of concurrent mpiruns [DEFAULT=1]')
    parser.add_argument('-nn', '--nnodes', type=int,
                        help='number of nodes per job, i.e. number of nodes \
                              using -nc cores [DEFAULT=1]. REQUIRES Intel MPI as \
                              found on e.g. Darwin HPC.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='debug output')
    parser.add_argument('--conv_cutoff', action='store_true',
                        help='run all res files at cutoff defined in cutoff.conv file')
    parser.add_argument('-l', '--limit', type=int,
                        help='limit to n structures per run')
    args = parser.parse_args()
    runner = BatchRun(ncores=args.ncores, nprocesses=args.nprocesses, nnodes=args.nnodes, debug=args.debug,
                      seed=args.seed, conv_cutoff=args.conv_cutoff, limit=args.limit)
    try:
        runner.spawn()
    except(KeyboardInterrupt, SystemExit):
        exit('Exiting top-level...')
