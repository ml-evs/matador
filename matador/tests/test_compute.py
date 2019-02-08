#!/usr/bin/env python

""" Some tests for high-throughput calculations. """

import unittest
import subprocess as sp
import shutil
import glob
import os
from os.path import realpath
from matador.compute.errors import CalculationError, MaxMemoryEstimateExceeded

HOSTNAME = os.uname()[1]
PATHS_TO_DEL = ['completed', 'bad_castep', 'input', 'logs', HOSTNAME]
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'
ROOT_DIR = os.getcwd()
VERBOSITY = 4
EXECUTABLE = 'castep'
print(80*'*')

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

    def test_missing_exec(self):
        """ Ensure failure if exec misses. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict
        from matador.compute.compute import CriticalError
        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s

        node = None
        nnodes = None
        seed = REAL_PATH + '/data/structures/LiAs_testcase.res'

        fall_over = False

        try:
            FullRelaxer(ncores=NCORES, nnodes=nnodes, node=node,
                        res=seed, param_dict=param_dict, cell_dict=cell_dict,
                        debug=False, verbosity=VERBOSITY, killcheck=True,
                        reopt=False, executable='THIS WAS MEANT TO FAIL, DON\'T WORRY',
                        start=True)
        except CriticalError:
            fall_over = True

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'As_00PBE.usp', 'LiAs_testcase.res']
        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        self.assertTrue(fall_over)

    def test_file_not_written(self):
        """ Run a calculation with an executable that only does "sleep" and
        check that run3 will stop the calculation early as no file is written.

        """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict
        os.chdir(REAL_PATH)
        if not os.path.isdir(REAL_PATH + '/missing_file_test'):
            os.makedirs(REAL_PATH + '/missing_file_test')
        os.chdir(REAL_PATH + '/missing_file_test')

        seed = REAL_PATH + 'data/symmetry_failure/Sb.res'
        cell_dict, s = cell2dict(REAL_PATH + '/data/symmetry_failure/KSb.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/symmetry_failure/KSb.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = REAL_PATH + 'data/missing_file_test/monkey_patch_sleep.sh'
        node = None

        relaxer = FullRelaxer(ncores=NCORES, nnodes=None, node=node,
                              res=seed, param_dict=param_dict, cell_dict=cell_dict,
                              debug=True, verbosity=VERBOSITY, executable=executable,
                              exec_test=False, compute_dir=None, polltime=1,
                              start=False)
        errored = False

        try:
            relaxer.relax()
        except CalculationError:
            errored = True

        os.chdir(REAL_PATH)
        from shutil import rmtree
        rmtree('missing_file_test')

        os.chdir(ROOT_DIR)
        self.assertTrue(relaxer.final_result is None)
        self.assertTrue(errored)

    def test_old_file(self):
        """ Run a calculation with an executable that only does "sleep", in the
        presence of a file that was written previouisly, and check that run3
        will stop the calculation early as no file is written.

        """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict
        os.chdir(REAL_PATH)
        if not os.path.isdir(REAL_PATH + '/missing_file_test'):
            os.makedirs(REAL_PATH + '/missing_file_test')
        os.chdir(REAL_PATH + '/missing_file_test')

        seed = REAL_PATH + 'data/symmetry_failure/Sb.res'
        with open('Sb.castep', 'w') as f:
            f.write('I am a CASTEP file, for sure.')

        cell_dict, s = cell2dict(REAL_PATH + '/data/symmetry_failure/KSb.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/symmetry_failure/KSb.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = REAL_PATH + 'data/missing_file_test/monkey_patch_sleep.sh'
        node = None

        relaxer = FullRelaxer(ncores=NCORES, nnodes=None, node=node,
                              res=seed, param_dict=param_dict, cell_dict=cell_dict,
                              debug=True, verbosity=VERBOSITY, executable=executable,
                              exec_test=False, compute_dir=None, polltime=1,
                              start=False)
        errored = False

        try:
            relaxer.relax()
        except CalculationError:
            errored = True

        os.chdir(REAL_PATH)
        from shutil import rmtree
        rmtree('missing_file_test')

        os.chdir(ROOT_DIR)
        self.assertTrue(relaxer.final_result is None)
        self.assertTrue(errored)

    def test_faked_error_recovery(self):
        """ Run a calculation that *should* throw a symmetry error, and try to
        recover from the error. If CASTEP is not present, monkey patch such that
        FullRelaxer copies the output files it would have expected.

        """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict
        os.chdir(REAL_PATH)
        if not os.path.isdir(REAL_PATH + '/fake_symmetry_test'):
            os.makedirs(REAL_PATH + '/fake_symmetry_test')
        os.chdir(REAL_PATH + '/fake_symmetry_test')

        seed = REAL_PATH + 'data/symmetry_failure/Sb.res'
        cell_dict, s = cell2dict(REAL_PATH + '/data/symmetry_failure/KSb.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/symmetry_failure/KSb.param', verbosity=VERBOSITY, db=False)
        assert s
        ncores = 1
        executable = REAL_PATH + 'data/symmetry_failure/monkey_patch_move.sh'
        node = None
        errored = False

        relaxer = FullRelaxer(ncores=ncores, nnodes=None, node=node,
                              res=seed, param_dict=param_dict, cell_dict=cell_dict,
                              debug=True, verbosity=VERBOSITY, executable=executable,
                              exec_test=False, compute_dir=None,
                              start=False)
        errored = False
        try:
            relaxer.relax()
        except CalculationError:
            errored = True

        bad_castep_exists = os.path.isdir('bad_castep')
        completed_exists = os.path.isdir('completed')

        os.chdir(REAL_PATH)
        from shutil import rmtree
        rmtree('fake_symmetry_test')

        os.chdir(ROOT_DIR)
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
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import res2dict, cell2dict, param2dict, castep2dict
        import multiprocessing as mp
        from time import sleep

        newborn, s = res2dict(REAL_PATH + '/data/structures/LiAs_testcase.res', verbosity=VERBOSITY, db=False)
        assert s
        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s

        node = None
        executable = 'castep'
        newborn['source'] = [REAL_PATH + '/data/GA_TESTCASE.res']

        os.chdir(REAL_PATH)
        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        queue = mp.Queue()
        relaxer = FullRelaxer(ncores=NCORES, nnodes=None, node=node,
                              res=newborn, param_dict=param_dict, cell_dict=cell_dict,
                              debug=False, verbosity=VERBOSITY, killcheck=True,
                              reopt=False, executable=executable, output_queue=queue,
                              start=False)
        # store proc object with structure ID, node name, output queue and number of cores
        proc = (1, node, mp.Process(target=relaxer.relax), NCORES)
        proc[2].start()
        while proc[2].is_alive():
            sleep(1)

        result, success = castep2dict('completed/GA_TESTCASE.castep')
        queue_result = queue.get()

        match_dict = dict()
        for key in queue_result:
            if key in ['source', 'site_occupancy', 'geom_iter']:
                continue
            match_dict[key] = (queue_result[key] == result[key])
            if not match_dict[key]:
                print(key, queue_result[key], result[key])

        print('Process completed!')

        completed_exists = os.path.isfile('completed/GA_TESTCASE.res')
        input_exists = os.path.isfile('input/GA_TESTCASE.res')

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'As_00PBE.usp']
        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        os.chdir(ROOT_DIR)
        self.assertTrue(completed_exists, "couldn't find output file!")
        self.assertTrue(input_exists, "couldn't find copy of input file!")
        self.assertTrue(success, "couldn't parse output file!")
        self.assertTrue(all([match_dict[key] for key in match_dict]))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_relax_to_file(self):
        """ Relax structure from file to file. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict

        os.chdir(REAL_PATH)
        seed = '_Li.res'
        shutil.copy(REAL_PATH + 'data/structures/Li.res', REAL_PATH + '_Li.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        FullRelaxer(ncores=NCORES, nnodes=None, node=node,
                    res=seed, param_dict=param_dict, cell_dict=cell_dict,
                    debug=False, verbosity=VERBOSITY, killcheck=True,
                    reopt=True, executable=executable,
                    start=True)

        print('Process completed!')

        completed_exists = os.path.isfile('completed/_Li.res')

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'As_00PBE.usp', REAL_PATH + 'data/_Li.res']
        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        os.chdir(ROOT_DIR)
        self.assertTrue(completed_exists, "couldn't find output file!")

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_failed_relaxation(self):
        """ Set a relaxation up to fail. """
        from matador.compute import FullRelaxer, reset_job_folder
        from matador.scrapers.castep_scrapers import cell2dict, param2dict

        seed = '_LiAs_testcase.res'
        os.chdir(REAL_PATH)
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase_bad.res', REAL_PATH + '_LiAs_testcase.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_tests/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_tests/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s
        param_dict['geom_max_iter'] = 3
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        relaxer = FullRelaxer(ncores=NCORES, nnodes=None, node=node,
                              res=seed, param_dict=param_dict, cell_dict=cell_dict,
                              debug=False, verbosity=VERBOSITY, killcheck=True, memcheck=False,
                              reopt=True, executable=executable, rough=0, fine_iter=3,
                              start=False)

        errored = False
        try:
            relaxer.relax()
        except CalculationError:
            errored = True
        self.assertTrue(errored, 'error not raised!')

        bad_exists = os.path.isfile('bad_castep/_LiAs_testcase.res')

        num = reset_job_folder()

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'As_00PBE.usp', REAL_PATH + 'data/_LiAs_testcase.res']
        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        os.chdir(ROOT_DIR)
        self.assertTrue(bad_exists, "couldn't find output file!")
        self.assertEqual(num, 0)

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_dont_restart_completed_calc(self):
        """ Set a relaxation up to fail. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict

        os.chdir(REAL_PATH + 'data/no_steps_left_todo')

        shutil.copy('cache/NaP_intermediates_stopped_early.res', '.')
        shutil.copy('cache/NaP_intermediates_stopped_early.castep', '.')

        cell_dict, s = cell2dict('NaP.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict('NaP.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = 'castep'
        node = None
        seed = 'NaP_intermediates_stopped_early'
        errored = False

        try:
            FullRelaxer(ncores=NCORES, nnodes=None, node=node,
                        res=seed, param_dict=param_dict, cell_dict=cell_dict,
                        debug=False, verbosity=VERBOSITY, killcheck=True, memcheck=False,
                        reopt=True, executable=executable,
                        start=True)
        except CalculationError:
            errored = True

        print('Process completed!')

        bad_exists = []
        bad_exists.append(os.path.isfile('bad_castep/NaP_intermediates_stopped_early.res'))
        bad_exists.append(os.path.isfile('bad_castep/NaP_intermediates_stopped_early.castep'))

        good_exists = []

        paths = ['input', 'bad_castep', 'logs']
        for path in paths:
            if os.path.isdir(path):
                good_exists.append(True)
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)
            else:
                good_exists.append(False)

        for _file in glob.glob('*.usp'):
            os.remove(_file)

        os.chdir(ROOT_DIR)
        self.assertTrue(errored)
        self.assertTrue(all(bad_exists))
        self.assertTrue(errored)
        self.assertTrue(all(good_exists))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_memcheck(self):
        """ Test the memory checker will not proceed with huge jobs. """
        from matador.scrapers.castep_scrapers import cell2dict, param2dict
        from matador.compute import FullRelaxer
        os.chdir(REAL_PATH)
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase.res', REAL_PATH + '/_LiAs_testcase.res')
        shutil.copy(REAL_PATH + 'data/LiAs_tests/LiAs.cell', REAL_PATH + '/LiAs.cell')
        shutil.copy(REAL_PATH + 'data/LiAs_tests/LiAs.param', REAL_PATH + '/LiAs.param')
        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', REAL_PATH + '/Li_00PBE.usp')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', REAL_PATH + '/As_00PBE.usp')

        cell_dict, s = cell2dict(REAL_PATH + 'LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + 'LiAs.param', verbosity=VERBOSITY, db=False)
        assert s

        raised_error = False
        try:
            FullRelaxer(ncores=NCORES, nnodes=None, node=None,
                        res='_LiAs_testcase', param_dict=param_dict, cell_dict=cell_dict,
                        debug=False, verbosity=VERBOSITY, killcheck=True, memcheck=True, maxmem=1,
                        start=False)
        except MaxMemoryEstimateExceeded:
            raised_error = True

        files_to_del = ['_LiAs_testcase.res', 'LiAs.cell', 'LiAs.param', 'Li_00PBE.usp', 'As_00PBE.usp']
        files_that_should_not_exist = ['_LiAs_testcase.res.lock', 'jobs.txt']
        folders_that_should_exist = ['logs']
        folders_that_should_not_exist = ['bad_castep', 'input', 'completed']

        correct_files = all([not os.path.isfile(_file) for _file in files_that_should_not_exist])
        correct_folders = all([os.path.isdir(folder) for folder in folders_that_should_exist])
        correct_folders *= all([not os.path.isdir(folder) for folder in folders_that_should_not_exist])

        for _f in files_to_del + files_that_should_not_exist:
            if os.path.isfile(_f):
                os.remove(_f)

        for _folder in folders_that_should_exist + folders_that_should_not_exist:
            if os.path.isdir(_folder):
                for _file in glob.glob(_folder + '/*'):
                    os.remove(_file)
                os.removedirs(_folder)

        os.chdir(ROOT_DIR)

        self.assertTrue(raised_error)
        self.assertTrue(correct_folders)
        self.assertTrue(correct_files)

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_batch_relax(self):
        """ Batch relax structures from file to file. """
        from matador.compute import BatchRun, reset_job_folder

        shutil.copy(REAL_PATH + 'data/structures/LiC.res', REAL_PATH + '_LiC.res')
        shutil.copy(REAL_PATH + 'data/LiC_tests/LiC.cell', REAL_PATH + 'LiC.cell')
        shutil.copy(REAL_PATH + 'data/LiC_tests/LiC.param', REAL_PATH + 'LiC.param')

        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', REAL_PATH + 'Li_00PBE.usp')
        shutil.copy(REAL_PATH + 'data/pspots/C_00PBE.usp', REAL_PATH + 'C_00PBE.usp')

        os.chdir(REAL_PATH)
        runner = BatchRun(seed=['LiC'], debug=False, no_reopt=True, verbosity=VERBOSITY, ncores=NCORES, executable=EXECUTABLE)
        runner.spawn(join=False)

        completed_exists = os.path.isfile('completed/_LiC.res')
        cruft = glob.glob('_LiC*')
        cruft_doesnt_exist = bool(len(cruft))

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'C_00PBE.usp',
                 REAL_PATH + '/LiC.cell',
                 REAL_PATH + '/LiC.param',
                 'jobs.txt', 'finished_cleanly.txt', 'failures.txt']

        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        num = reset_job_folder()
        os.chdir(ROOT_DIR)

        self.assertEqual(num, 0)
        self.assertTrue(completed_exists, "couldn't find output file!")
        self.assertFalse(cruft_doesnt_exist, "found some cruft {}".format(cruft))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_scf(self):
        """ Perform SCF on structure from file. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict

        os.chdir(REAL_PATH)
        seed = '_LiC.res'
        shutil.copy(REAL_PATH + 'data/structures/LiC.res', REAL_PATH + '_LiC.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiC_tests/LiC_scf.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiC_tests/LiC_scf.param', verbosity=VERBOSITY, db=False)
        assert s
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/C_00PBE.usp', '.')

        FullRelaxer(ncores=NCORES, nnodes=None, node=node,
                    res=seed, param_dict=param_dict, cell_dict=cell_dict,
                    debug=False, verbosity=VERBOSITY, killcheck=True,
                    reopt=True, executable=executable,
                    start=True)

        print('Process completed!')

        completed_exists = [os.path.isfile('completed/_LiC.res'),
                            os.path.isfile('completed/_LiC.castep'),
                            os.path.isfile('completed/_LiC-out.cell')]

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'C_00PBE.usp', REAL_PATH + 'data/_LiC.res']
        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        os.chdir(ROOT_DIR)
        self.assertTrue(all(completed_exists), "couldn't find output files!")

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_convergence_runner(self):
        """ Check that convergence tests run to completion. """
        from matador.compute import BatchRun
        os.chdir(REAL_PATH)
        shutil.copy(REAL_PATH + 'data/structures/Li.res', REAL_PATH + '_LiAs_testcase.res')
        shutil.copy(REAL_PATH + 'data/LiAs_tests/LiAs_scf.cell', REAL_PATH + 'LiAs_scf.cell')
        shutil.copy(REAL_PATH + 'data/LiAs_tests/LiAs_scf.param', REAL_PATH + 'LiAs_scf.param')
        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', REAL_PATH + 'Li_00PBE.usp')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', REAL_PATH + 'As_00PBE.usp')

        with open(REAL_PATH + 'kpt.conv', 'w') as f:
            f.write('0.08\n')
            f.write('0.07')

        with open(REAL_PATH + 'cutoff.conv', 'w') as f:
            f.write('300\n')
            f.write('400')

        runner = BatchRun(seed=['LiAs_scf'], debug=False,
                          conv_cutoff=True, conv_kpt=True,
                          verbosity=VERBOSITY, ncores=NCORES, nprocesses=2, executable=EXECUTABLE)
        runner.spawn(join=False)

        dirs_exist = [os.path.isdir(_dir) for _dir in ['completed_kpts', 'completed_cutoff']]
        bad_castep_exist = os.path.isdir('bad_castep')

        files_that_should_not_exist = ['_LiAs_testcase_300eV.res', '_LiAs_testcase_400eV.res',
                                       '_LiAs_testcase_0.08A.res', '_LiAs_testcase_0.07A.res']
        do_bad_files_exist = [os.path.isfile(_file) for _file in files_that_should_not_exist]

        results = ['completed_cutoff/_LiAs_testcase_{}eV.castep'.format(cutoff) for cutoff in [300, 400]]
        results += ['completed_kpts/_LiAs_testcase_{}A.castep'.format(kpts) for kpts in [0.08, 0.07]]
        files_exist = [os.path.isfile(_file) for _file in results]

        to_clean = ['jobs.txt', 'failures.txt', 'finished_cleanly.txt', 'kpt.conv', 'cutoff.conv',
                    'Li_00PBE.usp', 'As_00PBE.usp', '_LiAs_testcase.res', '_LiAs_testcase.res.lock',
                    'LiAs_scf.cell', 'LiAs_scf.param']

        for _file in to_clean:
            if os.path.isfile(REAL_PATH + _file):
                os.remove(REAL_PATH + _file)

        shutil.rmtree('completed_cutoff')
        shutil.rmtree('completed_kpts')
        shutil.rmtree('input')

        os.chdir(ROOT_DIR)

        self.assertTrue(all(dirs_exist))
        self.assertFalse(bad_castep_exist)
        self.assertTrue(all(files_exist))
        self.assertFalse(any(do_bad_files_exist))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_batch_failed_scf(self):
        """ Check that SCF failures don't kill everything... """
        from matador.compute import BatchRun, reset_job_folder

        os.chdir(REAL_PATH + 'data/fail_scf')
        shutil.copy(REAL_PATH + 'data/structures/Li.res', REAL_PATH + 'data/fail_scf/' + '_Li.res')
        shutil.copy(REAL_PATH + 'data/structures/LiC.res', REAL_PATH + 'data/fail_scf/' + '_LiC.res')
        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/C_00PBE.usp', '.')
        runner = BatchRun(seed=['LiC_scf'], debug=False, no_reopt=True,
                          verbosity=VERBOSITY, ncores=NCORES, executable=EXECUTABLE)
        runner.spawn()

        completed_folder_exists = os.path.isdir('completed')
        bad_castep_folder_exists = os.path.isdir('bad_castep')

        seeds = ['_Li.res', '_LiC.res']
        output_files_exist = all([os.path.isfile('bad_castep/{}'.format(seed)) for seed in seeds])

        cruft = glob.glob('_Li*')
        cruft_doesnt_exist = bool(len(cruft))

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'C_00PBE.usp',
                 'jobs.txt', 'finished_cleanly.txt', 'failures.txt']

        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        num = reset_job_folder()
        os.chdir(ROOT_DIR)

        self.assertEqual(num, 0)
        self.assertFalse(completed_folder_exists, "couldn't find output file!")
        self.assertTrue(bad_castep_folder_exists, "couldn't find bad_castep")
        self.assertTrue(output_files_exist, "couldn't find both outputs")
        self.assertFalse(cruft_doesnt_exist, "found some cruft {}".format(cruft))

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def test_batch_max_walltime_threaded(self):
        """ Check that WallTimeErrors do kill everything... """
        from matador.compute import BatchRun
        from matador.compute.compute import WalltimeError

        os.chdir(REAL_PATH + 'data/max_walltime')
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase.res', REAL_PATH + 'data/max_walltime/' + 'LiAs_testcase.res')
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase_bad.res', REAL_PATH + 'data/max_walltime/' + 'LiAs_testcase_bad.res')
        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')
        walltime_error = False
        runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True,
                          verbosity=VERBOSITY, ncores=2, nprocesses=2, executable=EXECUTABLE,
                          max_walltime=15, polltime=1)
        try:
            runner.spawn()
        except WalltimeError:
            walltime_error = True

        castep_exists = os.path.isfile('LiAs_testcase.castep')
        bad_castep_exists = os.path.isfile('LiAs_testcase_bad.castep')
        res_exists = os.path.isfile('LiAs_testcase.res') and os.path.isfile('LiAs_testcase_bad.res')
        lock_exists = os.path.isfile('LiAs_testcase.res.lock') and os.path.isfile('LiAs_testcase_bad.res.lock')
        compute_dir_exist = os.path.isdir(HOSTNAME)

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'As_00PBE.usp',
                 'jobs.txt', 'finished_cleanly.txt', 'failures.txt']

        paths += glob.glob('LiAs_testcase.*')
        paths += glob.glob('LiAs_testcase_bad.*')

        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        os.chdir(ROOT_DIR)

        self.assertTrue(walltime_error, 'Walltime error was not raised')
        self.assertTrue(castep_exists, 'Could not find castep file!')
        self.assertTrue(bad_castep_exists, 'Could not find bad castep file!')
        self.assertTrue(res_exists, 'Could not find res file!')
        self.assertFalse(compute_dir_exist, 'Compute dir not cleaned up!')
        self.assertFalse(lock_exists, 'Lock file was not deleted!')

    def test_batch_nothing_todo(self):
        """ Check that nothing is done when there's nothing to do... """
        from matador.compute import BatchRun
        import time

        os.chdir(REAL_PATH + 'data/nothing_to_do')
        runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True,
                          verbosity=VERBOSITY, ncores=2, nprocesses=2, executable=EXECUTABLE)
        start = time.time()
        runner.spawn()
        elapsed = time.time() - start
        os.chdir(ROOT_DIR)
        self.assertTrue(elapsed < 10, 'Sluggish to quit!')

    def test_res_name_collision(self):
        """ Check that run3 safely falls over if there is a file called <seed>.res. """
        from matador.compute import BatchRun
        from matador.compute.compute import InputError
        os.chdir(REAL_PATH + 'data/file_collision')
        try:
            runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True,
                              verbosity=VERBOSITY, ncores=2, nprocesses=2, executable=EXECUTABLE)
            runner.spawn()
        except InputError:
            failed_safely = True

        for path in PATHS_TO_DEL:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        os.chdir(ROOT_DIR)

        self.assertTrue(failed_safely)

    def test_missing_basics(self):
        """" Check that run3 falls over when e.g. xc_functional is missing. """
        from matador.compute import BatchRun
        from matador.compute.compute import InputError
        os.chdir(REAL_PATH + 'data/missing_basics')
        tests = ['missing_cutoff', 'missing_kpts', 'missing_xc', 'missing_pspot']
        errors = [False for test in tests]
        for ind, _ in enumerate(tests):
            try:
                runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True,
                                  verbosity=VERBOSITY, ncores=2, nprocesses=2, executable=EXECUTABLE)
                runner.spawn()
            except InputError:
                errors[ind] = True

        os.chdir(ROOT_DIR)
        self.assertTrue(all(errors))


if __name__ == '__main__':
    unittest.main()
