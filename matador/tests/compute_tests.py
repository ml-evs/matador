#!/usr/bin/env python
import unittest
from os.path import realpath
import subprocess as sp
import os
import shutil
import glob

hostname = os.uname()[1]
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'
ROOT_DIR = os.getcwd()
VERBOSITY = 0

try:
    with open('/dev/null', 'w') as f:
        proc = sp.Popen(['castep', '--version'], stdout=f, stderr=f)
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

        newborn, s = res2dict(REAL_PATH + '/data/LiAs_testcase.res', verbosity=VERBOSITY, db=False)
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
        shutil.copy(REAL_PATH + 'data/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/As_00PBE.usp', '.')

        os.chdir(REAL_PATH)

        relaxer = FullRelaxer(ncores=ncores, nnodes=None, node=node,
                              res=newborn, param_dict=param_dict, cell_dict=cell_dict,
                              debug=False, verbosity=VERBOSITY, killcheck=True,
                              reopt=False, executable=executable,
                              start=False)
        queue = mp.Queue()
        # store proc object with structure ID, node name, output queue and number of cores
        proc = (1,
                node,
                mp.Process(target=relaxer.relax,
                           args=(queue,)),
                ncores)
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
                        reopt=False, executable='definitely_fake_executable',
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
        shutil.copy(REAL_PATH + 'data/LiAs_testcase.res', REAL_PATH + '_LiAs_testcase.res')
        assert os.path.isfile(REAL_PATH + '_LiAs_testcase.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s
        ncores = 4
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/As_00PBE.usp', '.')

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
        from matador.compute import FullRelaxer, reset_job_folder_and_count_remaining
        from matador.scrapers.castep_scrapers import cell2dict, param2dict

        seed = '_LiAs_testcase.res'
        os.chdir(REAL_PATH)
        shutil.copy(REAL_PATH + 'data/LiAs_testcase_bad.res', REAL_PATH + '_LiAs_testcase.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs.param', verbosity=VERBOSITY, db=False)
        assert s
        param_dict['geom_max_iter'] = 3
        ncores = 4
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/As_00PBE.usp', '.')

        FullRelaxer(ncores=ncores, nnodes=None, node=node,
                    res=seed, param_dict=param_dict, cell_dict=cell_dict,
                    debug=False, verbosity=VERBOSITY, killcheck=True, memcheck=True,
                    reopt=True, executable=executable, rough=0, fine_iter=3,
                    start=True)

        print('Process completed!')

        bad_exists = os.path.isfile('bad_castep/_LiAs_testcase.res')

        num = reset_job_folder_and_count_remaining()

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
        """ Relax structure from file to file. """
        from matador.compute import BatchRun, reset_job_folder_and_count_remaining

        shutil.copy(REAL_PATH + 'data/LiAs_testcase.res', REAL_PATH + '_LiAs_testcase.res')
        shutil.copy(REAL_PATH + 'data/LiAs.cell', REAL_PATH + 'LiAs.cell')
        shutil.copy(REAL_PATH + 'data/LiAs.param', REAL_PATH + 'LiAs.param')

        shutil.copy(REAL_PATH + 'data/Li_00PBE.usp', REAL_PATH + 'Li_00PBE.usp')
        shutil.copy(REAL_PATH + 'data/As_00PBE.usp', REAL_PATH + 'As_00PBE.usp')

        os.chdir(REAL_PATH)
        runner = BatchRun(seed=['LiAs'], debug=False, no_reopt=True, verbosity=VERBOSITY, ncores=4)
        runner.spawn(join=True)

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

        num = reset_job_folder_and_count_remaining()
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
        shutil.copy(REAL_PATH + 'data/LiAs_testcase_bad.res', REAL_PATH + '_LiAs_testcase.res')
        assert os.path.isfile(REAL_PATH + '_LiAs_testcase.res')

        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs_scf.cell', verbosity=VERBOSITY, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs_scf.param', verbosity=VERBOSITY, db=False)
        assert s
        ncores = 4
        executable = 'castep'
        node = None

        shutil.copy(REAL_PATH + 'data/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/As_00PBE.usp', '.')

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


if __name__ == '__main__':
    unittest.main()
