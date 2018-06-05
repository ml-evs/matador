#!/usr/bin/env python
import unittest
from os.path import realpath
import subprocess as sp
import os
import shutil
import glob

HOSTNAME = os.uname()[1]
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'
ROOT_DIR = os.getcwd()
VERBOSITY = 0
EXECUTABLE = 'castep'

try:
    with open('/dev/null', 'w') as f:
        proc = sp.Popen([EXECUTABLE, '--version'], stdout=f, stderr=f)
        proc.communicate()
    if VERBOSITY > 0:
        print('Successfully detected CASTEP')
    CASTEP_PRESENT = True
except FileNotFoundError:
    if VERBOSITY > 0:
        print('Failed to detect CASTEP')
    CASTEP_PRESENT = False

try:
    with open('/dev/null', 'w') as f:
        proc = sp.Popen(['mpirun', '--version'], stdout=f, stderr=f)
        proc.communicate()
    if VERBOSITY > 0:
        print('Successfully detected mpirun')
    MPI_PRESENT = True
except FileNotFoundError:
    if VERBOSITY > 0:
        print('Failed to detect mpirun')
    MPI_PRESENT = False


class ComputeTest(unittest.TestCase):
    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def testRelaxToQueue(self):
        """ Mimic GA and test Queue relaxations. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import res2dict, cell2dict, param2dict, castep2dict
        import multiprocessing as mp
        from time import sleep

        newborn, s = res2dict(REAL_PATH + '/data/structures/LiAs_testcase.res', verbosity=VERBOSITY, db=False)
        assert s
        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s

        ncores = 4
        node = None
        executable = 'castep'
        newborn['source'] = [REAL_PATH + '/data/_LiAs_testcase.res']

        os.chdir(REAL_PATH)
        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        queue = mp.Queue()
        relaxer = FullRelaxer(ncores=ncores, nnodes=None, node=node,
                              res=newborn, param_dict=param_dict, cell_dict=cell_dict,
                              debug=False, verbosity=VERBOSITY, killcheck=True,
                              reopt=False, executable=executable, output_queue=queue,
                              start=False)
        # store proc object with structure ID, node name, output queue and number of cores
        proc = (1, node, mp.Process(target=relaxer.relax), ncores)
        proc[2].start()
        while proc[2].is_alive():
            sleep(1)

        result, success = castep2dict('completed/_LiAs_testcase.castep')
        queue_result = queue.get()

        match_dict = dict()
        for key in queue_result:
            if key == 'source':
                continue
            match_dict[key] = queue_result[key] == result[key]

        print('Process completed!')

        completed_exists = os.path.isfile('completed/_LiAs_testcase.res')
        input_exists = os.path.isfile('input/_LiAs_testcase.res')

        paths = ['completed', 'input', 'bad_castep']
        for path in paths:
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
    def testMissingExec(self):
        """ Ensure failure if exec misses. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict
        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s

        node = None
        nnodes = None
        seed = 'data/_LiAs_testcase.res'

        fall_over = False

        try:
            FullRelaxer(ncores=4, nnodes=nnodes, node=node,
                        res=seed, param_dict=param_dict, cell_dict=cell_dict,
                        debug=False, verbosity=VERBOSITY, killcheck=True,
                        reopt=False, executable='THIS WAS MEANT TO FAIL, DON\'T WORRY',
                        start=True)
        except SystemExit:
            fall_over = True

        paths = ['completed', 'input', 'bad_castep']
        for path in paths:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'As_00PBE.usp']
        for path in paths:
            if os.path.isfile(path):
                os.remove(path)

        self.assertTrue(fall_over)

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def testRelaxToFile(self):
        """ Relax structure from file to file. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict

        os.chdir(REAL_PATH)
        seed = '_LiAs_testcase.res'
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase.res', REAL_PATH + '_LiAs_testcase.res')
        assert os.path.isfile(REAL_PATH + '_LiAs_testcase.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s
        ncores = 4
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        FullRelaxer(ncores=ncores, nnodes=None, node=node,
                    res=seed, param_dict=param_dict, cell_dict=cell_dict,
                    debug=False, verbosity=VERBOSITY, killcheck=True,
                    reopt=True, executable=executable,
                    start=True)

        print('Process completed!')

        completed_exists = os.path.isfile('completed/_LiAs_testcase.res')

        paths = ['completed', 'input', 'bad_castep']
        for path in paths:
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
        self.assertTrue(completed_exists, "couldn't find output file!")

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def testFailedRelaxation(self):
        """ Set a relaxation up to fail. """
        from matador.compute import FullRelaxer, reset_job_folder
        from matador.scrapers.castep_scrapers import cell2dict, param2dict

        seed = '_LiAs_testcase.res'
        os.chdir(REAL_PATH)
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase_bad.res', REAL_PATH + '_LiAs_testcase.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s
        param_dict['geom_max_iter'] = 3
        ncores = 4
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        FullRelaxer(ncores=ncores, nnodes=None, node=node,
                    res=seed, param_dict=param_dict, cell_dict=cell_dict,
                    debug=False, verbosity=VERBOSITY, killcheck=True, memcheck=True,
                    reopt=True, executable=executable, rough=0, fine_iter=3,
                    start=True)

        print('Process completed!')

        bad_exists = os.path.isfile('bad_castep/_LiAs_testcase.res')

        num = reset_job_folder()

        paths = ['completed', 'input', 'bad_castep']
        for path in paths:
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
    def testBatchRelax(self):
        """ Batch relax structures from file to file. """
        from matador.compute import BatchRun, reset_job_folder

        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase.res', REAL_PATH + '_LiAs_testcase.res')
        shutil.copy(REAL_PATH + 'data/LiAs.cell', REAL_PATH + 'LiAs.cell')
        shutil.copy(REAL_PATH + 'data/LiAs.param', REAL_PATH + 'LiAs.param')

        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', REAL_PATH + 'Li_00PBE.usp')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', REAL_PATH + 'As_00PBE.usp')

        os.chdir(REAL_PATH)
        runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True, verbosity=VERBOSITY, ncores=4, executable=EXECUTABLE)
        runner.spawn(join=False)

        completed_exists = os.path.isfile('completed/_LiAs_testcase.res')
        cruft = glob.glob('_LiAs_testcase*')
        cruft_doesnt_exist = bool(len(cruft))

        paths = ['completed', 'input', 'bad_castep']
        for path in paths:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'As_00PBE.usp',
                 REAL_PATH + '/LiAs.cell',
                 REAL_PATH + '/LiAs.param',
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
    def testSCF(self):
        """ Perform SCF on structure from file. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import cell2dict, param2dict

        os.chdir(REAL_PATH)
        seed = '_LiAs_testcase.res'
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase_bad.res', REAL_PATH + '_LiAs_testcase.res')
        assert os.path.isfile(REAL_PATH + '_LiAs_testcase.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_scf.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_scf.param', verbosity=VERBOSITY, db=False)
        assert s
        ncores = 4
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')

        FullRelaxer(ncores=ncores, nnodes=None, node=node,
                    res=seed, param_dict=param_dict, cell_dict=cell_dict,
                    debug=False, verbosity=VERBOSITY, killcheck=True,
                    reopt=True, executable=executable,
                    start=True)

        print('Process completed!')

        completed_exists = [os.path.isfile('completed/_LiAs_testcase.res'),
                            os.path.isfile('completed/_LiAs_testcase.castep'),
                            os.path.isfile('completed/_LiAs_testcase-out.cell')]

        paths = ['completed', 'input', 'bad_castep']
        for path in paths:
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
        self.assertTrue(all(completed_exists), "couldn't find output files!")

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def testConvergenceRunner(self):
        """ Check that convergence tests run to completion. """
        from matador.compute import BatchRun
        os.chdir(REAL_PATH)
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase.res', REAL_PATH + '_LiAs_testcase.res')
        shutil.copy(REAL_PATH + 'data/LiAs_scf.cell', REAL_PATH + 'LiAs_scf.cell')
        shutil.copy(REAL_PATH + 'data/LiAs_scf.param', REAL_PATH + 'LiAs_scf.param')
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
                          verbosity=VERBOSITY, ncores=4, nprocesses=2, executable=EXECUTABLE)
        runner.spawn(join=False)

        dirs_exist = [os.path.isdir(_dir) for _dir in ['completed_kpts', 'completed_cutoff']]
        bad_castep_exist = os.path.isdir('bad_castep')

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

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def testBatchFailedSCF(self):
        """ Check that SCF failures don't kill everything... """
        from matador.compute import BatchRun, reset_job_folder

        os.chdir(REAL_PATH + 'data/fail_scf')
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase_bad.res', REAL_PATH + 'data/fail_scf/' + 'LiAs_testcase_bad.res')
        shutil.copy(REAL_PATH + 'data/structures/LiAs_testcase.res', REAL_PATH + 'data/fail_scf/' + 'LiAs_testcase.res')
        shutil.copy(REAL_PATH + 'data/pspots/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/pspots/As_00PBE.usp', '.')
        runner = BatchRun(seed=['LiAs_scf'], debug=False, no_reopt=True, verbosity=VERBOSITY, ncores=4, executable=EXECUTABLE)
        runner.spawn(join=False)

        completed_folder_exists = os.path.isdir('completed')
        bad_castep_folder_exists = os.path.isdir('bad_castep')

        seeds = ['LiAs_testcase.res', 'LiAs_testcase_bad.res']
        output_files_exist = all([os.path.isfile('bad_castep/{}'.format(seed)) for seed in seeds])

        cruft = glob.glob('LiAs_testcase*')
        cruft_doesnt_exist = bool(len(cruft))

        paths = ['completed', 'input', 'bad_castep']
        for path in paths:
            if os.path.isdir(path):
                files = glob.glob(path + '/*')
                for file in files:
                    os.remove(file)
                os.removedirs(path)

        paths = ['Li_00PBE.usp', 'As_00PBE.usp',
                 REAL_PATH + '/LiAs.cell',
                 REAL_PATH + '/LiAs.param',
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
    def testBatchMaxWallTimeThreaded(self):
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
                          max_walltime=5, polltime=1)
        try:
            runner.spawn()
        except WalltimeError as err:
            walltime_error = True

        castep_exists = os.path.isfile('LiAs_testcase.castep')
        bad_castep_exists = os.path.isfile('LiAs_testcase_bad.castep')
        res_exists = os.path.isfile('LiAs_testcase.res') and os.path.isfile('LiAs_testcase_bad.res')
        lock_exists = os.path.isfile('LiAs_testcase.res.lock') and os.path.isfile('LiAs_testcase_bad.res.lock')
        compute_dir_doesnt_exist = os.path.isdir(HOSTNAME)

        paths = ['completed', 'input', 'bad_castep', HOSTNAME]
        for path in paths:
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
        self.assertFalse(lock_exists, 'Lock file was not deleted!')

    @unittest.skipIf((not CASTEP_PRESENT or not MPI_PRESENT), 'castep or mpirun executable not found in PATH')
    def testBatchNothingToDo(self):
        """ Check that nothing is done when there's nothing to do... """
        from matador.compute import BatchRun
        import time

        os.chdir(REAL_PATH + 'data/nothing_to_do')
        runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True,
                          verbosity=VERBOSITY, ncores=2, nprocesses=2, executable=EXECUTABLE)
        start = time.time()
        runner.spawn()
        elapsed = time.time() - start
        self.assertTrue(elapsed < 10, 'Sluggish to quit!')


if __name__ == '__main__':
    unittest.main()