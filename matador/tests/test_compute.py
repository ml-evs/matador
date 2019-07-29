#!/usr/bin/env python

""" Some tests for high-throughput calculations. """

import unittest
import subprocess as sp
import glob
import time

from shutil import copy
from os import getcwd, uname
from os.path import realpath, isfile, isdir


from matador.compute.errors import (
    CalculationError, MaxMemoryEstimateExceeded, CriticalError,
    WalltimeError, InputError
)
from matador.compute import ComputeTask, BatchRun, reset_job_folder
from matador.scrapers.castep_scrapers import cell2dict, param2dict, res2dict, castep2dict

HOSTNAME = uname()[1]
PATHS_TO_DEL = ['completed', 'bad_castep', 'input', 'logs', HOSTNAME]
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'
TMP_DIR = 'tmp_test'
ROOT_DIR = getcwd()
VERBOSITY = 4
EXECUTABLE = 'castep'


try:
    with open('/dev/null', 'w') as devnull:
        sp.Popen([EXECUTABLE, '--version'], stdout=devnull, stderr=devnull).communicate()
    if VERBOSITY > 0:
        print('Successfully detected CASTEP')
    CASTEP_PRESENT = True
except FileNotFoundError:
    if VERBOSITY > 0:
        print('Failed to detect CASTEP')
    CASTEP_PRESENT = False

try:
    with open('/dev/null', 'w') as devnull:
        sp.Popen(['mpirun', '--version'], stdout=devnull, stderr=devnull).communicate()
    if VERBOSITY > 0:
        print('Successfully detected mpirun')
    MPI_PRESENT = True
except FileNotFoundError:
    if VERBOSITY > 0:
        print('Failed to detect mpirun')
    MPI_PRESENT = False

if CASTEP_PRESENT and MPI_PRESENT:
    NCORES = 4
else:
    NCORES = 1


class ComputeTest(unittest.TestCase):
    """ Run tests equivalent to using the run3 script for
    various artificial setups.

    """
    def tearDown(self):
        from shutil import rmtree
        from os import chdir
        chdir(REAL_PATH)
        if isdir(TMP_DIR):
            rmtree(TMP_DIR)

        chdir(ROOT_DIR)

    def setUp(self):
        from os import chdir, makedirs
        chdir(REAL_PATH)
        makedirs(TMP_DIR, exist_ok=True)
        chdir(TMP_DIR)

    def test_missing_exec(self):
        """ Ensure failure if exec misses. """
        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s

        node = None
        nnodes = None
        seed = REAL_PATH + '/data/structures/LiAs_testcase.res'

        fall_over = False

        try:
            ComputeTask(ncores=NCORES, nnodes=nnodes, node=node,
                        res=seed, param_dict=param_dict, cell_dict=cell_dict,
                        verbosity=VERBOSITY, killcheck=True,
                        reopt=False, executable='THIS WAS MEANT TO FAIL, DON\'T WORRY',
                        start=True)
        except CriticalError:
            fall_over = True

        self.assertTrue(fall_over)

    def test_file_not_written(self):
        """ Run a calculation with an executable that only does "sleep" and
        check that run3 will stop the calculation early as no file is written.

        """
        seed = REAL_PATH + 'data/symmetry_failure/Sb.res'
        cell_dict, s = cell2dict(REAL_PATH + '/data/symmetry_failure/KSb.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/symmetry_failure/KSb.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = REAL_PATH + 'data/missing_file_test/monkey_patch_sleep.sh'
        node = None

        relaxer = ComputeTask(ncores=NCORES, nnodes=None, node=node,
                              res=seed, param_dict=param_dict, cell_dict=cell_dict,
                              verbosity=VERBOSITY, executable=executable,
                              exec_test=False, compute_dir=None, polltime=1,
                              start=False)
        errored = False

        try:
            relaxer.relax()
        except CalculationError:
            errored = True

        self.assertTrue(relaxer.final_result is None)
        self.assertTrue(errored)

    def test_old_file(self):
        """ Run a calculation with an executable that only does "sleep", in the
        presence of a file that was written previouisly, and check that run3
        will stop the calculation early as no file is written.

        """
        seed = REAL_PATH + 'data/symmetry_failure/Sb.res'
        with open('Sb.castep', 'w') as f:
            f.write('I am a CASTEP file, for sure.')

        cell_dict, s = cell2dict(REAL_PATH + '/data/symmetry_failure/KSb.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/symmetry_failure/KSb.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = REAL_PATH + 'data/missing_file_test/monkey_patch_sleep.sh'
        node = None

        relaxer = ComputeTask(ncores=NCORES, nnodes=None, node=node,
                              res=seed, param_dict=param_dict, cell_dict=cell_dict,
                              verbosity=VERBOSITY, executable=executable,
                              exec_test=False, compute_dir=None, polltime=1,
                              start=False)
        errored = False

        try:
            relaxer.relax()
        except CalculationError:
            errored = True

        self.assertTrue(relaxer.final_result is None)
        self.assertTrue(errored)

    def test_faked_error_recovery(self):
        """ Run a calculation that *should* throw a symmetry error, and try to
        recover from the error. If CASTEP is not present, monkey patch such that
        ComputeTask copies the output files it would have expected.

        """
        seed = REAL_PATH + 'data/symmetry_failure/Sb.res'
        cell_dict, s = cell2dict(REAL_PATH + '/data/symmetry_failure/KSb.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/symmetry_failure/KSb.param', verbosity=VERBOSITY, db=False)
        assert s
        ncores = 1
        executable = REAL_PATH + 'data/symmetry_failure/monkey_patch_move.sh'
        node = None
        errored = False

        relaxer = ComputeTask(ncores=ncores, nnodes=None, node=node,
                              res=seed, param_dict=param_dict, cell_dict=cell_dict,
                              verbosity=VERBOSITY, executable=executable,
                              exec_test=False, compute_dir=None,
                              start=False)
        errored = False
        try:
            relaxer.relax()
        except CalculationError:
            errored = True

        bad_castep_exists = isdir('bad_castep')
        completed_exists = isdir('completed')

        self.assertTrue(errored)
        self.assertTrue(bad_castep_exists)
        self.assertFalse(completed_exists)
        self.assertTrue(relaxer.final_result is None)
        self.assertTrue(errored)
        self.assertEqual(relaxer._num_retries, 3)
        self.assertTrue('symmetry_generate' not in relaxer.calc_doc)
        self.assertTrue('snap_to_symmetry' not in relaxer.calc_doc)
        self.assertTrue('symmetry_tol' not in relaxer.calc_doc)

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_relax_to_queue(self):
        """ Mimic GA and test Queue relaxations. """
        import multiprocessing as mp

        newborn, s = res2dict(REAL_PATH + '/data/structures/LiAs_testcase.res', verbosity=VERBOSITY, db=False)
        assert s
        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s

        node = None
        executable = 'castep'
        newborn['source'] = [REAL_PATH + '/data/GA_TESTCASE.res']

        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        queue = mp.Queue()
        relaxer = ComputeTask(ncores=NCORES, nnodes=None, node=node,
                              res=newborn, param_dict=param_dict, cell_dict=cell_dict,
                              verbosity=VERBOSITY, killcheck=True,
                              reopt=False, executable=executable, output_queue=queue,
                              start=False)
        # store proc object with structure ID, node name, output queue and number of cores
        proc = (1, node, mp.Process(target=relaxer.relax), NCORES)
        proc[2].start()
        while proc[2].is_alive():
            time.sleep(1)

        result, success = castep2dict('completed/GA_TESTCASE.castep')
        queue_result = queue.get()

        match_dict = dict()
        for key in queue_result:
            if key in ['source', 'site_occupancy', 'geom_iter']:
                continue
            match_dict[key] = (queue_result[key] == result[key])
            if not match_dict[key]:
                print(key, queue_result[key], result[key])

        completed_exists = isfile('completed/GA_TESTCASE.res')
        input_exists = isfile('input/GA_TESTCASE.res')

        self.assertTrue(completed_exists, "couldn't find output file!")
        self.assertTrue(input_exists, "couldn't find copy of input file!")
        self.assertTrue(success, "couldn't parse output file!")
        self.assertTrue(all([match_dict[key] for key in match_dict]))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_relax_to_file(self):
        """ Relax structure from file to file. """
        seed = '_Li.res'
        copy(REAL_PATH + 'data/structures/Li.res', '_Li.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = 'castep'
        node = None

        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        ComputeTask(ncores=NCORES, nnodes=None, node=node,
                    res=seed, param_dict=param_dict, cell_dict=cell_dict,
                    verbosity=VERBOSITY, killcheck=True,
                    reopt=True, executable=executable, exec_test=False,
                    start=True)

        print('Process completed!')

        completed_exists = isfile('completed/_Li.res')
        self.assertTrue(completed_exists, "couldn't find output file!")

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_failed_relaxation(self):
        """ Set a relaxation up to fail. """
        seed = '_LiAs_testcase.res'
        copy(REAL_PATH + 'data/structures/LiAs_testcase_bad.res', '_LiAs_testcase.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s
        param_dict['geom_max_iter'] = 3
        executable = 'castep'
        node = None

        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        relaxer = ComputeTask(ncores=NCORES, nnodes=None, node=node,
                              res=seed, param_dict=param_dict, cell_dict=cell_dict,
                              verbosity=VERBOSITY, killcheck=True, memcheck=False,
                              reopt=True, executable=executable, rough=0, fine_iter=3,
                              start=False)

        errored = False
        try:
            relaxer.relax()
        except CalculationError:
            errored = True
        self.assertTrue(errored, 'error not raised!')

        bad_exists = isfile('bad_castep/_LiAs_testcase.res')

        num = reset_job_folder()

        self.assertTrue(bad_exists, "couldn't find output file!")
        self.assertEqual(num, 0)

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_dont_restart_completed_calc(self):
        """ Set a relaxation up to fail. """

        copy(REAL_PATH + 'data/no_steps_left_todo/cache/NaP_intermediates_stopped_early.res', '.')
        copy(REAL_PATH + 'data/no_steps_left_todo/cache/NaP_intermediates_stopped_early.castep', '.')

        cell_dict, s = cell2dict(REAL_PATH + 'data/no_steps_left_todo/NaP.cell', verbosity=VERBOSITY, db=False)
        self.assertTrue(s)
        param_dict, s = param2dict(REAL_PATH + 'data/no_steps_left_todo/NaP.param', verbosity=VERBOSITY, db=False)
        self.assertTrue(s)
        executable = 'castep'
        node = None
        seed = 'NaP_intermediates_stopped_early'
        errored = False

        try:
            ComputeTask(ncores=NCORES, nnodes=None, node=node,
                        res=seed, param_dict=param_dict, cell_dict=cell_dict,
                        verbosity=VERBOSITY, killcheck=True, memcheck=False,
                        reopt=True, executable=executable,
                        start=True)
        except CalculationError:
            errored = True

        print('Process completed!')

        bad_exists = []
        bad_exists.append(isfile('bad_castep/NaP_intermediates_stopped_early.res'))
        bad_exists.append(isfile('bad_castep/NaP_intermediates_stopped_early.castep'))
        bad_exists = all(bad_exists)
        good_exists = all(isdir(path) for path in ['input', 'bad_castep', 'logs'])

        self.assertTrue(errored)
        self.assertTrue(bad_exists)
        self.assertTrue(errored)
        self.assertTrue(good_exists)

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_memcheck(self):
        """ Test the memory checker will not proceed with huge jobs. """

        copy(REAL_PATH + 'data/structures/LiAs_testcase.res', '_LiAs_testcase.res')
        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        cell_dict, s = cell2dict(REAL_PATH + 'data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        self.assertTrue(s)
        param_dict, s = param2dict(REAL_PATH + 'data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        self.assertTrue(s)

        raised_error = False
        try:
            ComputeTask(ncores=NCORES, nnodes=None, node=None,
                        res='_LiAs_testcase', param_dict=param_dict, cell_dict=cell_dict,
                        verbosity=VERBOSITY, killcheck=True, memcheck=True, maxmem=1,
                        start=False)
        except MaxMemoryEstimateExceeded:
            raised_error = True

        files_that_should_not_exist = ['_LiAs_testcase.res.lock', 'jobs.txt']
        folders_that_should_exist = ['logs']
        folders_that_should_not_exist = ['bad_castep', 'input', 'completed']

        correct_files = all([not isfile(_file) for _file in files_that_should_not_exist])
        correct_folders = all([isdir(folder) for folder in folders_that_should_exist])
        correct_folders *= all([not isdir(folder) for folder in folders_that_should_not_exist])

        self.assertTrue(raised_error)
        self.assertTrue(correct_folders)
        self.assertTrue(correct_files)

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_batch_relax(self):
        """ Batch relax structures from file to file. """

        copy(REAL_PATH + 'data/structures/LiC.res', '_LiC.res')
        copy(REAL_PATH + 'data/LiC_tests/LiC.cell', '.')
        copy(REAL_PATH + 'data/LiC_tests/LiC.param', '.')

        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/C_00PBE.usp', '.')

        runner = BatchRun(seed=['LiC'], debug=False, no_reopt=True, verbosity=VERBOSITY, ncores=NCORES, executable=EXECUTABLE)
        runner.spawn(join=False)

        completed_exists = isfile('completed/_LiC.res')
        cruft = glob.glob('_LiC*')
        cruft_doesnt_exist = bool(len(cruft))

        num = reset_job_folder()

        self.assertEqual(num, 0)
        self.assertTrue(completed_exists, "couldn't find output file!")
        self.assertFalse(cruft_doesnt_exist, "found some cruft {}".format(cruft))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_scf(self):
        """ Perform SCF on structure from file. """
        seed = '_LiC.res'
        copy(REAL_PATH + 'data/structures/LiC.res', '_LiC.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiC_tests/LiC_scf.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiC_tests/LiC_scf.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = 'castep'
        node = None

        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/C_00PBE.usp', '.')

        ComputeTask(ncores=NCORES, nnodes=None, node=node,
                    res=seed, param_dict=param_dict, cell_dict=cell_dict,
                    verbosity=VERBOSITY, killcheck=True,
                    reopt=True, executable=executable,
                    start=True)

        completed_exists = [isfile('completed/_LiC.res'),
                            isfile('completed/_LiC.castep'),
                            isfile('completed/_LiC-out.cell')]

        self.assertTrue(all(completed_exists), "couldn't find output files!")

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_convergence_runner(self):
        """ Check that convergence tests run to completion. """
        copy(REAL_PATH + 'data/structures/Li.res', '_LiAs_testcase.res')
        copy(REAL_PATH + 'data/LiAs_tests/LiAs_scf.cell', '.')
        copy(REAL_PATH + 'data/LiAs_tests/LiAs_scf.param', '.')
        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        with open('kpt.conv', 'w') as f:
            f.write('0.08\n')
            f.write('0.07')

        with open('cutoff.conv', 'w') as f:
            f.write('300\n')
            f.write('400')

        runner = BatchRun(seed=['LiAs_scf'], debug=False,
                          conv_cutoff=True, conv_kpt=True,
                          verbosity=VERBOSITY, ncores=NCORES, nprocesses=2, executable=EXECUTABLE)
        runner.spawn(join=False)

        dirs_exist = [isdir(_dir) for _dir in ['completed_kpts', 'completed_cutoff']]
        bad_castep_exist = isdir('bad_castep')

        files_that_should_not_exist = ['_LiAs_testcase_300eV.res', '_LiAs_testcase_400eV.res',
                                       '_LiAs_testcase_0.08A.res', '_LiAs_testcase_0.07A.res']
        do_bad_files_exist = [isfile(_file) for _file in files_that_should_not_exist]

        results = ['completed_cutoff/_LiAs_testcase_{}eV.castep'.format(cutoff) for cutoff in [300, 400]]
        results += ['completed_kpts/_LiAs_testcase_{}A.castep'.format(kpts) for kpts in [0.08, 0.07]]
        files_exist = [isfile(_file) for _file in results]

        self.assertTrue(all(dirs_exist))
        self.assertFalse(bad_castep_exist)
        self.assertTrue(all(files_exist))
        self.assertFalse(any(do_bad_files_exist))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_batch_failed_scf(self):
        """ Check that SCF failures don't kill everything... """

        copy(REAL_PATH + 'data/structures/Li.res', '_Li.res')
        copy(REAL_PATH + 'data/structures/LiC.res', '_LiC.res')
        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/C_00PBE.usp', '.')
        copy(REAL_PATH + 'data/fail_scf/LiC_scf.cell', '.')
        copy(REAL_PATH + 'data/fail_scf/LiC_scf.param', '.')

        runner = BatchRun(seed=['LiC_scf'], debug=False, no_reopt=True,
                          verbosity=VERBOSITY, ncores=NCORES, executable=EXECUTABLE)
        runner.spawn()

        completed_folder_exists = isdir('completed')
        bad_castep_folder_exists = isdir('bad_castep')

        seeds = ['_Li.res', '_LiC.res']
        output_files_exist = all([isfile('bad_castep/{}'.format(seed)) for seed in seeds])

        cruft = glob.glob('_Li*')
        cruft_doesnt_exist = bool(len(cruft))

        num = reset_job_folder()

        self.assertEqual(num, 0)
        self.assertFalse(completed_folder_exists, "couldn't find output file!")
        self.assertTrue(bad_castep_folder_exists, "couldn't find bad_castep")
        self.assertTrue(output_files_exist, "couldn't find both outputs")
        self.assertFalse(cruft_doesnt_exist, "found some cruft {}".format(cruft))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_batch_max_walltime_threaded(self):
        """ Check that WallTimeErrors do kill everything... """

        copy(REAL_PATH + 'data/structures/LiAs_testcase.res', '.')
        copy(REAL_PATH + 'data/structures/LiAs_testcase_bad.res', '.')
        copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')
        copy(REAL_PATH + 'data/max_walltime/LiAs.cell', '.')
        copy(REAL_PATH + 'data/max_walltime/LiAs.param', '.')

        walltime_error = False
        runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True,
                          verbosity=VERBOSITY, ncores=2, nprocesses=2, executable=EXECUTABLE,
                          max_walltime=15, polltime=1)
        try:
            runner.spawn()
        except WalltimeError:
            walltime_error = True

        castep_exists = isfile('LiAs_testcase.castep')
        bad_castep_exists = isfile('LiAs_testcase_bad.castep')
        res_exists = isfile('LiAs_testcase.res') and isfile('LiAs_testcase_bad.res')
        lock_exists = isfile('LiAs_testcase.res.lock') and isfile('LiAs_testcase_bad.res.lock')
        compute_dir_exist = isdir(HOSTNAME)

        self.assertTrue(walltime_error, 'Walltime error was not raised')
        self.assertTrue(castep_exists, 'Could not find castep file!')
        self.assertTrue(bad_castep_exists, 'Could not find bad castep file!')
        self.assertTrue(res_exists, 'Could not find res file!')
        self.assertFalse(compute_dir_exist, 'Compute dir not cleaned up!')
        self.assertFalse(lock_exists, 'Lock file was not deleted!')

    def test_generic_batch(self):
        """ Run a calculation in generic mode, check that res files
        are cycled over.

        """

        executable = 'echo'
        files = glob.glob(REAL_PATH + '/data/structures/*.res')
        for _file in files:
            copy(_file, '.')

        runner = BatchRun(seed='*.res', debug=False, mode='generic',
                          verbosity=4, ncores=1, executable=executable)

        runner.spawn(join=True)

        completed_files_exist = all([isfile('completed/' + _file.split('/')[-1]) for _file in files])
        txt_files_exist = all([isfile(_file) for _file in ['jobs.txt', 'finished_cleanly.txt']])
        dirs = ['completed', 'input', 'logs']
        dirs_exist = all([isdir(_dir) for _dir in dirs])

        logs = glob.glob('logs/*.log')
        num_logs = len(logs)
        log_lines = []
        for log in logs:
            with open(log, 'r') as f:
                log_lines.append(len(f.readlines()))

        self.assertEqual(num_logs, len(files), msg='Not enough log files!')
        self.assertTrue(all(lines > 5 for lines in log_lines), msg='Log files were too short!')
        self.assertTrue(completed_files_exist)
        self.assertTrue(dirs_exist)
        self.assertTrue(txt_files_exist)

    def test_batch_nothing_todo(self):
        """ Check that nothing is done when there's nothing to do... """
        for file in glob.glob(REAL_PATH + 'data/nothing_to_do/*.*'):
            copy(file, '.')

        runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True,
                          verbosity=VERBOSITY, ncores=2, nprocesses=2, executable=EXECUTABLE)
        start = time.time()
        runner.spawn()
        elapsed = time.time() - start
        self.assertTrue(elapsed < 10, 'Sluggish to quit!')

    def test_res_name_collision(self):
        """ Check that run3 safely falls over if there is a file called <seed>.res. """
        copy(REAL_PATH + 'data/file_collision/LiAs.cell', '.')
        copy(REAL_PATH + 'data/file_collision/LiAs.res', '.')
        copy(REAL_PATH + 'data/file_collision/LiAs.param', '.')

        failed_safely = False
        try:
            runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True,
                              verbosity=VERBOSITY, ncores=2, nprocesses=2, executable=EXECUTABLE)
            runner.spawn()
        except InputError:
            failed_safely = True

        self.assertTrue(failed_safely)

    def test_missing_basics(self):
        """" Check that run3 falls over when e.g. xc_functional is missing. """
        for file in glob.glob(REAL_PATH + 'data/misisng_basics/*.*'):
            copy(file, '.')

        tests = ['missing_cutoff', 'missing_kpts', 'missing_xc', 'missing_pspot']
        errors = [False for test in tests]
        for ind, test in enumerate(tests):
            try:
                runner = BatchRun(seed=['LiC_' + test], debug=False, no_reopt=True,
                                  verbosity=VERBOSITY, executable=EXECUTABLE)
                runner.spawn()
            except InputError:
                errors[ind] = True

        self.assertTrue(all(errors))


if __name__ == '__main__':
    unittest.main()
