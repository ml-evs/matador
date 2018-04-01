#!/usr/bin/env python
import unittest
from os.path import realpath
import os
import shutil
import glob

hostname = os.uname()[1]
REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


@unittest.skipIf(hostname != 'cluster2', 'You probably don\'t have CASTEP!')
class ComputeTest(unittest.TestCase):
    def testRelaxToQueue(self):
        """ Mimic GA and test Queue relaxations. """
        from matador.compute import FullRelaxer
        from matador.scrapers.castep_scrapers import res2dict, cell2dict, param2dict
        import multiprocessing as mp
        from time import sleep

        newborn, s = res2dict(REAL_PATH + '/data/LiAs_testcase.res', verbosity=5, db=False)
        assert s
        cell_dict, s = cell2dict(REAL_PATH + '/data/LiAs.cell', verbosity=5, db=False)
        assert s
        param_dict, s = param2dict(REAL_PATH + '/data/LiAs.param', verbosity=5, db=False)
        assert s
        print('GEOM_MAX_ITER: {}'.format(param_dict['geom_max_iter']))
        ncores = 16
        verbosity = 10
        executable = 'castep'
        node = 'node1'

        newborn['source'] = [REAL_PATH + '/data/_LiAs_testcase.res']

        shutil.copy(REAL_PATH + 'data/Li_00PBE.usp', '.')
        shutil.copy(REAL_PATH + 'data/As_00PBE.usp', '.')

        relaxer = FullRelaxer(ncores=ncores, nnodes=None, node=node,
                              res=newborn, param_dict=param_dict, cell_dict=cell_dict,
                              debug=False, verbosity=verbosity, killcheck=True,
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
            sleep(10)

        print('Process completed!')

        assert os.path.isfile('completed/LiAs_testcase.res')

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
