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
from export import doc2cell, doc2param, doc2res
from traceback import print_exc
import bson.json_util as json
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
        valid = True
        if self.args.get('nprocesses') is not None:
            self.nprocesses = self.args['nprocesses']
        else:
            self.nprocesses = 1
        if self.args.get('ncores') is not None:
            self.ncores = self.args['ncores']
        else:
            self.ncores = self.all_cores / self.nprocesses
        if self.ncores*self.nprocesses > self.all_cores:
            print('Requested more cores (' + str(self.ncores*self.nprocesses) +
                  ') than available (' + str(self.all_cores) + ').')
            exit('Exiting...')
        # scan directory for files to run
        self.file_lists = defaultdict(list)
        for root, dirs, files in walk('.'):
            for file in files:
                if root == '.':
                    if file.endswith('.res'):
                        self.file_lists['res'].append(file)
        if len(self.file_lists['res']) < 1:
            valid = False
            print('run3 requires at least 1 res file in folder, found',
                  len(self.file_lists['res']))
        self.cell_dict, cell_success = cell2dict(self.seed + '.cell', db=False)
        if not cell_success:
            valid = False
            print('Failed to parse cell file')
        self.param_dict, param_success = param2dict(self.seed + '.param', db=False)
        if not param_success:
            valid = False
            print('Failed to parse cell file')
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
        for res in res_list:
            running = False
            with open(paths['jobs_fname'], 'rw') as job_file:
                flines = job_file.readlines()
                for line in flines:
                    if res in line:
                        running = True
                        break
            if not running:
                with open(paths['jobs_fname'], 'a') as job_file:
                    job_file.write(res+'\n')
                # create full relaxer object for creation and running of job
                try:
                    FullRelaxer(paths=self.paths,
                                ncores=self.ncores,
                                res=res,
                                param_dict=self.param_dict,
                                cell_dict=self.cell_dict,
                                debug=self.debug)
                    print('Completed', res)
                except(KeyboardInterrupt, SystemExit, RuntimeError):
                    raise SystemExit
        return

class FullRelaxer:
    """ Perform full relxation of res input by first doing 
    4 rough optimisations with only a few iterations, followed by
    4 larger optimisations with many iterations,
    e.g. 4 lots of 2 then 4 lots of geom_max_iter/4.
    """
    def __init__(self, paths, ncores, res, param_dict, cell_dict, debug):
        """ Make the files to run the calculation and handle
        the calling of CASTEP itself.
        """
        self.paths = paths
        self.ncores = ncores
        self.executable = 'castep'
        self.debug = debug
        # read in initial structure
        res_dict, success = res2dict(res, db=False)
        calc_doc = res_dict
        # set seed name
        self.seed = calc_doc['source'][0].replace('.res', '')
        print('Relaxing', self.seed)
        # update global doc with cell and param dicts for folder
        calc_doc.update(cell_dict)
        calc_doc.update(param_dict)
        calc_doc['task'] = 'geometryoptimization'
        if self.debug:
            print(json.dumps(calc_doc, indent=2))
        self.success = self.relax(calc_doc)
        if not success:
            self.mv_to_bad()

    def relax(self, calc_doc):
        """ Set up the calculation to perform 4 sets of two steps,
        then continue with the remainder of steps.
        """
        seed = self.seed
        geom_max_iter_list = [2, 2, 2, 2, calc_doc['geom_max_iter']]
        # relax structure
        for ind, num_iter in enumerate(geom_max_iter_list):
            calc_doc['geom_max_iter'] = num_iter
            try:
                # delete any existing files
                if isfile(self.seed + '.param'):
                    remove(seed+'.param')
                if isfile(seed + '.cell'):
                    remove(seed+'.cell')
                # write new param and cell
                doc2param(calc_doc, seed, hash_dupe=False)
                doc2cell(calc_doc, seed, hash_dupe=False, copy_pspots=False)
                # run CASTEP
                process = self.castep()
                process.communicate()
                # scrape new structure from castep file
                opti_dict, success = castep2dict(seed + '.castep', db=False)
                if not success and opti_dict == '':
                    print('Failed to scrape castep file...')
                    return False
                if opti_dict['optimised'] == True:
                    if not exists('completed'):
                        makedirs('completed')
                    print('Successfully relaxed', seed)
                    # write res and castep file out to completed folder
                    doc2res(opti_dict, 'completed/' + seed)
                    system('mv ' + seed + '.castep' + ' completed/' + seed + '.castep')
                    system('mv ' + seed + '.param' + ' completed/' + seed + '.param')
                    system('mv ' + seed + '.cell' + ' completed/' + seed + '.cell')
                    system('mv ' + seed + '-out.cell' + ' completed/' + seed + '-out.cell')
                    # clean up rest of files
                    self.tidy_up()
                    return True
                elif ind == len(geom_max_iter_list) - 1:
                    print('Failed to optimise', self.seed)
                    self.mv_to_bad()
                else:
                    err_file = seed + '*.err'
                    for globbed in glob.glob(err_file):
                        if isfile(globbed):
                            self.mv_to_bad()
                            return False
                    calc_doc.update(opti_dict)
            except:
                print_exc()
                self.mv_to_bad()
                self.tidy_up()
                return False

    def mv_to_bad(self):
        """ Move all associated files to bad_castep. """
        if not exists('bad_castep'):
            makedirs('bad_castep')
        print('Something went wrong, moving files to bad_castep')
        system('mv ' + self.seed + '* bad_castep')
        return

    def tidy_up(self):
        """ Delete all run3 created files before quitting. """
        system('rm ' + self.seed + '*')
        return

    def castep(self):
        """ Calls CASTEP on desired seed with desired number of cores. 
        Errors piped to /dev/null for now...
        """
        if self.ncores == 1:
            process = sp.Popen(['nice', '-n', '15', self.executable, self.seed])
        else:
            process = sp.Popen(['nice', '-n', '15', 'mpirun', '-n', str(self.ncores),
                                self.executable, self.seed])
        return process


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='run3',
        description='Run multiple CASTEP geometry optmizations from a series of .res \
                     files and single cell and param files. The calculation will be \
                     restarted after 3 geom_opt steps by default, or otherwise by \
                     value of geom_max_iter in the param file.',
        epilog='Written by Matthew Evans (2016), based primarily on run.pl and run2.pl \
                by Chris Pickard and Andrew Morris and PyAIRSS CastepRunner by Jamie Wynn.')
    parser.add_argument('seed', type=str,
                        help='cell and param seed to use as template for calculations')
    parser.add_argument('-nc', '--ncores', type=int,
                        help='number of cores CASTEP per job [DEFAULT=cpu_count/nprocesses]')
    parser.add_argument('-np', '--nprocesses', type=int,
                        help='number of concurrent calculations [DEFAULT=1]')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='debug output')
    args = parser.parse_args()
    runner = BatchRun(ncores=args.ncores, nprocesses=args.nprocesses, debug=args.debug, seed=args.seed)
    try:
        runner.spawn()
    except(KeyboardInterrupt, SystemExit):
        exit('Exiting top-level...')
