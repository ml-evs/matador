#!/usr/bin/env python

""" Some tests for high-throughput calculations. """

import unittest
import os
from os.path import isfile, isdir
import glob
import copy
import shutil
import time
import warnings
import multiprocessing as mp

import psutil

from matador.utils.errors import (
    CalculationError,
    MaxMemoryEstimateExceeded,
    CriticalError,
    WalltimeError,
    InputError,
)
from matador.compute import ComputeTask, BatchRun, reset_job_folder
from matador.compute.slurm import SlurmQueueManager
from matador.compute.pbs import PBSQueueManager
from matador.scrapers.castep_scrapers import (
    cell2dict,
    param2dict,
    res2dict,
    castep2dict,
)
from .utils import REAL_PATH, MatadorUnitTest, detect_program

HOSTNAME = os.uname()[1]
PATHS_TO_DEL = ["completed", "bad_castep", "input", "logs", HOSTNAME]
VERBOSITY = 10
EXECUTABLE = "castep"
RUN_SLOW_TESTS = HOSTNAME == "cluster2"

CASTEP_PRESENT = detect_program(EXECUTABLE)
MPI_PRESENT = detect_program("mpirun")

ACTUAL_NCORES = psutil.cpu_count(logical=False)
if CASTEP_PRESENT and MPI_PRESENT:
    NCORES = max(psutil.cpu_count(logical=False) - 2, 1)
else:
    NCORES = 1


class ComputeTest(MatadorUnitTest):
    """Run tests equivalent to using the run3 script for
    various artificial setups.

    """

    def test_missing_exec(self):
        """Ensure failure if exec misses."""
        cell_dict, s = cell2dict(
            REAL_PATH + "/data/LiAs_tests/LiAs.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/LiAs_tests/LiAs.param", verbosity=VERBOSITY, db=False
        )
        assert s

        node = None
        nnodes = None
        seed = REAL_PATH + "/data/structures/LiAs_testcase.res"

        fall_over = False

        try:
            ComputeTask(
                ncores=NCORES,
                nnodes=nnodes,
                node=node,
                res=seed,
                param_dict=param_dict,
                cell_dict=cell_dict,
                verbosity=VERBOSITY,
                killcheck=True,
                reopt=False,
                executable="THIS WAS MEANT TO FAIL, DON'T WORRY",
                start=True,
            )
        except CriticalError:
            fall_over = True

        self.assertTrue(fall_over)

    def test_file_not_written(self):
        """Run a calculation with an executable that only does "sleep" and
        check that run3 will stop the calculation early as no file is written.

        """
        seed = REAL_PATH + "data/symmetry_failure/Sb.res"
        cell_dict, s = cell2dict(
            REAL_PATH + "/data/symmetry_failure/KSb.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/symmetry_failure/KSb.param",
            verbosity=VERBOSITY,
            db=False,
        )
        assert s
        executable = REAL_PATH + "data/missing_file_test/monkey_patch_sleep.sh"
        node = None

        relaxer = ComputeTask(
            ncores=NCORES,
            nnodes=None,
            node=node,
            res=seed,
            param_dict=param_dict,
            cell_dict=cell_dict,
            verbosity=VERBOSITY,
            executable=executable,
            exec_test=False,
            polltime=1,
            start=False,
        )
        with self.assertRaises(CalculationError):
            relaxer.begin()

        self.assertTrue(relaxer.final_result is None)

    def test_old_file(self):
        """Run a calculation with an executable that only does "sleep", in the
        presence of a file that was written previouisly, and check that run3
        will stop the calculation early as no file is written.

        """
        seed = REAL_PATH + "data/symmetry_failure/Sb.res"
        with open("Sb.castep", "w") as f:
            f.write("I am a CASTEP file, for sure.")

        cell_dict, s = cell2dict(
            REAL_PATH + "/data/symmetry_failure/KSb.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/symmetry_failure/KSb.param",
            verbosity=VERBOSITY,
            db=False,
        )
        assert s
        executable = REAL_PATH + "data/missing_file_test/monkey_patch_sleep.sh"
        node = None

        relaxer = ComputeTask(
            ncores=NCORES,
            nnodes=None,
            node=node,
            res=seed,
            param_dict=param_dict,
            cell_dict=cell_dict,
            verbosity=VERBOSITY,
            executable=executable,
            exec_test=False,
            polltime=1,
            start=False,
        )

        with self.assertRaises(CalculationError):
            relaxer.begin()
        self.assertTrue(relaxer.final_result is None)

    def test_faked_error_recovery(self):
        """Run a calculation that *should* throw a symmetry error, and try to
        recover from the error. If CASTEP is not present, monkey patch such that
        ComputeTask copies the output files it would have expected.

        """
        seed = REAL_PATH + "data/symmetry_failure/Sb.res"
        cell_dict, s = cell2dict(
            REAL_PATH + "/data/symmetry_failure/KSb.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/symmetry_failure/KSb.param",
            verbosity=VERBOSITY,
            db=False,
        )
        assert s
        ncores = 1
        executable = REAL_PATH + "data/symmetry_failure/monkey_patch_move.sh"
        node = None

        relaxer = ComputeTask(
            ncores=ncores,
            nnodes=None,
            node=node,
            res=seed,
            param_dict=param_dict,
            cell_dict=cell_dict,
            verbosity=VERBOSITY,
            executable=executable,
            exec_test=False,
            polltime=1,
            start=False,
        )

        with self.assertRaises(CalculationError):
            relaxer.begin()

        bad_castep_exists = isdir("bad_castep")
        completed_exists = isdir("completed")

        self.assertTrue(bad_castep_exists)
        self.assertFalse(completed_exists)
        self.assertTrue(relaxer.final_result is None)
        self.assertEqual(relaxer._num_retries, 3)
        self.assertTrue("symmetry_generate" not in relaxer.calc_doc)
        self.assertTrue("snap_to_symmetry" not in relaxer.calc_doc)
        self.assertTrue("symmetry_tol" not in relaxer.calc_doc)

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_relax_to_queue(self):
        """Mimic GA and test Queue relaxations."""

        newborn, s = res2dict(
            REAL_PATH + "/data/structures/LiAs_testcase.res",
            verbosity=VERBOSITY,
            db=False,
        )
        assert s
        cell_dict, s = cell2dict(
            REAL_PATH + "/data/LiAs_tests/LiAs.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/LiAs_tests/LiAs.param", verbosity=VERBOSITY, db=False
        )
        assert s

        node = None
        executable = "castep"
        newborn["source"] = [REAL_PATH + "/data/GA_TESTCASE.res"]

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/As_00PBE.usp", ".")

        queue = mp.Queue()
        relaxer = ComputeTask(
            ncores=NCORES,
            nnodes=None,
            node=node,
            res=newborn,
            param_dict=param_dict,
            cell_dict=cell_dict,
            verbosity=VERBOSITY,
            killcheck=True,
            reopt=False,
            executable=executable,
            output_queue=queue,
            start=True,
        )
        # store proc object with structure ID, node name, output queue and number of cores
        proc = (1, node, mp.Process(target=relaxer.relax), NCORES)
        proc[2].start()
        while proc[2].is_alive():
            time.sleep(1)

        result, success = castep2dict("completed/GA_TESTCASE.castep")
        queue_result = queue.get()

        match_dict = dict()
        for key in queue_result:
            if key in ["source", "site_occupancy", "geom_iter"]:
                continue
            match_dict[key] = queue_result[key] == result[key]
            if not match_dict[key]:
                print(key, queue_result[key], result[key])

        completed_exists = isfile("completed/GA_TESTCASE.res")
        input_exists = isfile("input/GA_TESTCASE.res")

        self.assertTrue(completed_exists, "couldn't find output file!")
        self.assertTrue(input_exists, "couldn't find shutil.copy of input file!")
        self.assertTrue(success, "couldn't parse output file!")
        self.assertTrue(all([match_dict[key] for key in match_dict]))

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_relax_to_file(self):
        """Relax structure from file to file."""
        seed = "_Li.res"
        shutil.copy(REAL_PATH + "data/structures/Li.res", "_Li.res")

        cell_dict, s = cell2dict(
            REAL_PATH + "/data/LiAs_tests/LiAs.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/LiAs_tests/LiAs.param", verbosity=VERBOSITY, db=False
        )
        assert s
        executable = "castep"
        node = None

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/As_00PBE.usp", ".")

        ComputeTask(
            ncores=NCORES,
            nnodes=None,
            node=node,
            res=seed,
            param_dict=param_dict,
            cell_dict=cell_dict,
            verbosity=VERBOSITY,
            killcheck=True,
            reopt=True,
            executable=executable,
            exec_test=False,
            start=True,
        )

        print("Process completed!")

        completed_exists = isfile("completed/_Li.res")
        base_files_exist = [
            isfile("_Li.res"),
            isfile("_Li.res.lock"),
            isfile("_Li.castep"),
        ]
        self.assertTrue(completed_exists, "couldn't find output file!")
        self.assertFalse(any(base_files_exist), "couldn't clean input files")

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_failed_relaxation(self):
        """Set a relaxation up to fail."""
        seed = "_LiAs_testcase.res"
        shutil.copy(
            REAL_PATH + "data/structures/LiAs_testcase_bad.res", "_LiAs_testcase.res"
        )

        cell_dict, s = cell2dict(
            REAL_PATH + "/data/LiAs_tests/LiAs.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/LiAs_tests/LiAs.param", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict["geom_max_iter"] = 3
        executable = "castep"
        node = None

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/As_00PBE.usp", ".")

        relaxer = ComputeTask(
            ncores=NCORES,
            nnodes=None,
            node=node,
            res=seed,
            param_dict=param_dict,
            cell_dict=cell_dict,
            verbosity=VERBOSITY,
            killcheck=True,
            memcheck=False,
            reopt=True,
            executable=executable,
            rough=0,
            fine_iter=3,
            start=False,
        )

        with self.assertRaises(CalculationError):
            relaxer.begin()

        bad_exists = isfile("bad_castep/_LiAs_testcase.res")

        num = reset_job_folder()

        self.assertTrue(bad_exists, "couldn't find output file!")
        self.assertEqual(num, 0)

    def test_dont_restart_completed_calc(self):
        """Set a relaxation up to fail."""

        shutil.copy(
            REAL_PATH
            + "data/no_steps_left_todo/cache/NaP_intermediates_stopped_early.res",
            ".",
        )
        shutil.copy(
            REAL_PATH
            + "data/no_steps_left_todo/cache/NaP_intermediates_stopped_early.castep",
            ".",
        )

        cell_dict, s = cell2dict(
            REAL_PATH + "data/no_steps_left_todo/NaP.cell",
            verbosity=VERBOSITY,
            db=False,
        )
        self.assertTrue(s)
        param_dict, s = param2dict(
            REAL_PATH + "data/no_steps_left_todo/NaP.param",
            verbosity=VERBOSITY,
            db=False,
        )
        self.assertTrue(s)
        executable = "castep"
        node = None
        seed = "NaP_intermediates_stopped_early.res"
        with self.assertRaises(CalculationError):
            ComputeTask(
                ncores=NCORES,
                nnodes=None,
                node=node,
                res=seed,
                param_dict=param_dict,
                cell_dict=cell_dict,
                verbosity=VERBOSITY,
                killcheck=True,
                memcheck=False,
                exec_test=False,
                reopt=True,
                executable=executable,
                max_walltime=5,
                start=True,
            )

        print("Process completed!")

        bad_exists = []
        bad_exists.append(isfile("bad_castep/NaP_intermediates_stopped_early.res"))
        bad_exists.append(isfile("bad_castep/NaP_intermediates_stopped_early.castep"))
        bad_exists = all(bad_exists)
        good_exists = all(isdir(path) for path in ["input", "bad_castep", "logs"])

        self.assertTrue(bad_exists)
        self.assertTrue(good_exists)

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_memcheck(self):
        """Test the memory checker will not proceed with huge jobs."""

        shutil.copy(
            REAL_PATH + "data/structures/LiAs_testcase.res", "_LiAs_testcase.res"
        )
        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/As_00PBE.usp", ".")

        cell_dict, s = cell2dict(
            REAL_PATH + "data/LiAs_tests/LiAs.cell", verbosity=VERBOSITY, db=False
        )
        self.assertTrue(s)
        param_dict, s = param2dict(
            REAL_PATH + "data/LiAs_tests/LiAs.param", verbosity=VERBOSITY, db=False
        )
        self.assertTrue(s)

        with self.assertRaises(MaxMemoryEstimateExceeded):
            ComputeTask(
                ncores=NCORES,
                nnodes=None,
                node=None,
                res="_LiAs_testcase.res",
                param_dict=param_dict,
                cell_dict=cell_dict,
                verbosity=VERBOSITY,
                killcheck=True,
                memcheck=True,
                maxmem=1,
                start=True,
            )

        files_that_should_not_exist = ["_LiAs_testcase.res.lock", "jobs.txt"]
        folders_that_should_exist = ["logs"]
        folders_that_should_not_exist = ["bad_castep", "input", "completed"]

        correct_files = all(
            [not isfile(_file) for _file in files_that_should_not_exist]
        )
        correct_folders = all([isdir(folder) for folder in folders_that_should_exist])
        correct_folders *= all(
            [not isdir(folder) for folder in folders_that_should_not_exist]
        )

        self.assertTrue(correct_folders)
        self.assertTrue(correct_files)

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    @unittest.skipIf(not RUN_SLOW_TESTS, "this is a slow test, skipping")
    def test_batch_relax(self):
        """Batch relax structures from file to file."""

        shutil.copy(REAL_PATH + "data/structures/LiC.res", "_LiC.res")
        shutil.copy(REAL_PATH + "data/LiC_tests/LiC.cell", ".")
        shutil.copy(REAL_PATH + "data/LiC_tests/LiC.param", ".")

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/C_00PBE.usp", ".")

        runner = BatchRun(
            seed=["LiC"],
            debug=False,
            no_reopt=True,
            verbosity=VERBOSITY,
            ncores=NCORES,
            executable=EXECUTABLE,
        )
        runner.spawn(join=False)

        completed_exists = isfile("completed/_LiC.res")
        cruft = glob.glob("_LiC*")
        cruft_doesnt_exist = bool(len(cruft))

        num = reset_job_folder()

        self.assertEqual(num, 0)
        self.assertTrue(completed_exists, "couldn't find output file!")
        self.assertFalse(cruft_doesnt_exist, "found some cruft {}".format(cruft))

    def test_batch_queues(self):
        """Test the scraping of queuing environments."""

        shutil.copy(REAL_PATH + "data/structures/LiC.res", "_LiC.res")
        shutil.copy(REAL_PATH + "data/LiC_tests/LiC.cell", ".")
        shutil.copy(REAL_PATH + "data/LiC_tests/LiC.param", ".")

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/C_00PBE.usp", ".")

        slurm_env = {
            "SLURM_NTASKS": "120",
            "SLURM_JOB_ID": "123456",
            "SLURM_ARRAY_TASK_ID": "123123",
            "SLURM_MEM_PER_CPU": "1024",
            "SLURM_RANDOM_STRING": "hello",
            "blah": "abc",
        }

        pbs_env = {"PBS_TASKNUM": "120", "PBS_JOB_ID": "999", "PBS_ARRAYID": "123123"}

        old_env = copy.deepcopy(os.environ)
        os.environ.update(slurm_env)

        runner = BatchRun(
            seed=["LiC"],
            debug=False,
            ncores=NCORES,
            verbosity=VERBOSITY,
            executable=EXECUTABLE,
        )

        self.assertEqual(runner.args["ncores"], NCORES)
        self.assertEqual(type(runner.queue_mgr), SlurmQueueManager)
        self.assertEqual(runner.maxmem, 1024 * 120)
        self.assertEqual(runner.queue_mgr.max_memory, 1024 * 120)
        self.assertEqual(runner.queue_mgr.array_id, 123123)
        self.assertEqual(runner.queue_mgr.env["SLURM_RANDOM_STRING"], "hello")
        self.assertTrue("blah" not in runner.queue_mgr.env)

        runner = BatchRun(
            seed=["LiC"],
            debug=False,
            ncores=None,
            verbosity=VERBOSITY,
            executable=EXECUTABLE,
        )

        self.assertEqual(runner.args["ncores"], 120)
        self.assertEqual(type(runner.queue_mgr), SlurmQueueManager)
        self.assertEqual(runner.maxmem, 1024 * 120)
        self.assertEqual(runner.queue_mgr.max_memory, 1024 * 120)
        self.assertEqual(runner.queue_mgr.array_id, 123123)

        os.environ = copy.deepcopy(old_env)
        os.environ.update(pbs_env)

        runner = BatchRun(
            seed=["LiC"],
            debug=False,
            ncores=None,
            verbosity=VERBOSITY,
            executable=EXECUTABLE,
        )

        print(runner.queue_mgr)

        self.assertEqual(runner.args["ncores"], 120)
        self.assertEqual(type(runner.queue_mgr), PBSQueueManager)
        self.assertEqual(runner.maxmem, None)
        self.assertEqual(runner.queue_mgr.max_memory, None)
        self.assertEqual(runner.queue_mgr.array_id, 123123)

        os.environ.update(slurm_env)
        with self.assertRaises(RuntimeError):
            runner = BatchRun(
                seed=["LiC"],
                debug=False,
                ncores=None,
                verbosity=VERBOSITY,
                executable=EXECUTABLE,
            )

        os.environ = copy.deepcopy(old_env)

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_scf(self):
        """Perform SCF on structure from file."""
        seed = "_LiC.res"
        shutil.copy(REAL_PATH + "data/structures/LiC.res", "_LiC.res")

        cell_dict, s = cell2dict(
            REAL_PATH + "/data/LiC_tests/LiC_scf.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/LiC_tests/LiC_scf.param", verbosity=VERBOSITY, db=False
        )
        assert s
        executable = "castep"
        node = None

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/C_00PBE.usp", ".")

        ComputeTask(
            ncores=NCORES,
            nnodes=None,
            node=node,
            res=seed,
            param_dict=param_dict,
            cell_dict=cell_dict,
            verbosity=VERBOSITY,
            killcheck=True,
            reopt=True,
            executable=executable,
            compute_dir="/tmp/compute_test",
            start=True,
        )

        completed_exists = [
            isfile("completed/_LiC.res"),
            isfile("completed/_LiC.castep"),
            isfile("completed/_LiC-out.cell"),
        ]
        base_file_exists = [
            isfile("_LiC.res"),
            isfile("_LiC.castep"),
            isfile("_LiC.res.lock"),
        ]

        self.assertFalse(any(base_file_exists), "failed to clean up files!")
        self.assertTrue(all(completed_exists), "couldn't find output files!")

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_scf_max_walltime(self):
        """Perform SCF on structure from file."""
        seed = "_LiC.res"
        shutil.copy(REAL_PATH + "data/structures/LiC.res", "_LiC.res")

        cell_dict, s = cell2dict(
            REAL_PATH + "/data/LiC_tests/LiC_scf.cell", verbosity=VERBOSITY, db=False
        )
        assert s
        param_dict, s = param2dict(
            REAL_PATH + "/data/LiC_tests/LiC_scf.param", verbosity=VERBOSITY, db=False
        )
        assert s
        executable = "castep"
        node = None

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/C_00PBE.usp", ".")

        with self.assertRaises(WalltimeError):
            ComputeTask(
                ncores=NCORES,
                nnodes=None,
                node=node,
                res=seed,
                param_dict=param_dict,
                cell_dict=cell_dict,
                verbosity=VERBOSITY,
                timings=(5, time.time()),
                polltime=2,
                executable=executable,
                compute_dir="/tmp/compute_test",
                start=True,
            )

        base_file_exists = [
            isfile("_LiC.res"),
            isfile("_LiC.castep"),
        ]

        self.assertFalse(isfile("_LiC.res.lock"), "failed to clean up lock")
        self.assertTrue(all(base_file_exists), "failed to keep valid files!")

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_convergence_runner(self):
        """Check that convergence tests run to completion."""
        shutil.copy(REAL_PATH + "data/structures/Li.res", "_LiAs_testcase.res")
        shutil.copy(REAL_PATH + "data/LiAs_tests/LiAs_scf.cell", ".")
        shutil.copy(REAL_PATH + "data/LiAs_tests/LiAs_scf.param", ".")
        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/As_00PBE.usp", ".")

        with open("kpt.conv", "w") as f:
            f.write("0.08\n")
            f.write("0.07")

        with open("cutoff.conv", "w") as f:
            f.write("300\n")
            f.write("400")

        runner = BatchRun(
            seed=["LiAs_scf"],
            debug=False,
            conv_cutoff=True,
            conv_kpt=True,
            verbosity=VERBOSITY,
            ncores=1,
            nprocesses=min(2, ACTUAL_NCORES),
            executable=EXECUTABLE,
        )
        runner.spawn(join=False)

        dirs_exist = [isdir(_dir) for _dir in ["completed_kpts", "completed_cutoff"]]
        bad_castep_exist = isdir("bad_castep")

        files_that_should_not_exist = [
            "_LiAs_testcase_300eV.res",
            "_LiAs_testcase_400eV.res",
            "_LiAs_testcase_0.08A.res",
            "_LiAs_testcase_0.07A.res",
        ]
        do_bad_files_exist = [isfile(_file) for _file in files_that_should_not_exist]

        results = [
            "completed_cutoff/_LiAs_testcase_{}eV.castep".format(cutoff)
            for cutoff in [300, 400]
        ]
        results += [
            "completed_kpts/_LiAs_testcase_{}A.castep".format(kpts)
            for kpts in [0.08, 0.07]
        ]
        files_exist = [isfile(_file) for _file in results]

        self.assertTrue(all(dirs_exist))
        self.assertFalse(bad_castep_exist)
        self.assertTrue(all(files_exist))
        self.assertFalse(any(do_bad_files_exist))

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_batch_failed_scf(self):
        """Check that SCF failures don't kill everything..."""

        shutil.copy(REAL_PATH + "data/structures/Li.res", "_Li.res")
        shutil.copy(REAL_PATH + "data/structures/LiC.res", "_LiC.res")
        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/C_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/fail_scf/LiC_scf.cell", ".")
        shutil.copy(REAL_PATH + "data/fail_scf/LiC_scf.param", ".")

        runner = BatchRun(
            seed=["LiC_scf"],
            debug=False,
            no_reopt=True,
            verbosity=VERBOSITY,
            ncores=NCORES,
            executable=EXECUTABLE,
        )
        runner.spawn()

        completed_folder_exists = isdir("completed")
        bad_castep_folder_exists = isdir("bad_castep")

        seeds = ["_Li.res", "_LiC.res"]
        output_files_exist = all(
            [isfile("bad_castep/{}".format(seed)) for seed in seeds]
        )

        cruft = glob.glob("_Li*")
        cruft_doesnt_exist = bool(len(cruft))

        compute_dir_exist = isdir(HOSTNAME)

        num = reset_job_folder()

        self.assertEqual(num, 0)
        self.assertFalse(completed_folder_exists, "couldn't find output file!")
        self.assertTrue(bad_castep_folder_exists, "couldn't find bad_castep")
        self.assertTrue(output_files_exist, "couldn't find both outputs")
        self.assertFalse(cruft_doesnt_exist, "found some cruft {}".format(cruft))
        self.assertFalse(compute_dir_exist, "Compute dir not cleaned up!")

    def test_failed_compute_dir_scf(self):
        """Check that using a garbage path for compute dir causes a safe crash."""
        shutil.copy(REAL_PATH + "data/structures/Li.res", "_Li.res")
        shutil.copy(REAL_PATH + "data/structures/LiC.res", "_LiC.res")
        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/C_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/fail_scf/LiC_scf.cell", ".")
        shutil.copy(REAL_PATH + "data/fail_scf/LiC_scf.param", ".")

        runner = BatchRun(
            seed=["LiC_scf"],
            debug=False,
            no_reopt=True,
            scratch_prefix="/this/drive/doesnt/exist",
            verbosity=VERBOSITY,
            ncores=NCORES,
            executable=EXECUTABLE,
        )
        with self.assertRaises(CriticalError):
            runner.spawn()

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_batch_max_walltime_threaded(self):
        """Check that WallTimeErrors do kill everything..."""

        shutil.copy(REAL_PATH + "data/structures/LiAs_testcase.res", ".")
        shutil.copy(REAL_PATH + "data/structures/LiAs_testcase_bad.res", ".")
        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/As_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/max_walltime/LiAs.cell", ".")
        shutil.copy(REAL_PATH + "data/max_walltime/LiAs.param", ".")

        runner = BatchRun(
            seed=["LiAs"],
            debug=False,
            no_reopt=True,
            verbosity=VERBOSITY,
            ncores=min(2, ACTUAL_NCORES),
            nprocesses=min(2, ACTUAL_NCORES // 2),
            executable=EXECUTABLE,
            max_walltime=15,
            polltime=10,
        )
        with self.assertRaises(WalltimeError):
            runner.spawn()

        castep_exists = isfile("LiAs_testcase.castep")
        bad_castep_exists = isfile("LiAs_testcase_bad.castep")
        res_exists = isfile("LiAs_testcase.res") and isfile("LiAs_testcase_bad.res")
        lock_exists = isfile("LiAs_testcase.res.lock") or isfile(
            "LiAs_testcase_bad.res.lock"
        )
        compute_dir_exist = isdir(HOSTNAME)

        self.assertTrue(castep_exists, "Could not find castep file!")
        self.assertTrue(bad_castep_exists, "Could not find bad castep file!")
        self.assertTrue(res_exists, "Could not find res file!")
        self.assertFalse(compute_dir_exist, "Compute dir not cleaned up!")
        self.assertFalse(lock_exists, "Lock file was not deleted!")

    def test_generic_batch(self):
        """Run a calculation in generic mode, check that res files
        are cycled over.

        """

        executable = "echo"
        files = glob.glob(REAL_PATH + "/data/structures/*.res")
        for _file in files:
            shutil.copy(_file, ".")

        runner = BatchRun(
            seed="*.res",
            debug=False,
            mode="generic",
            verbosity=4,
            ncores=1,
            executable=executable,
        )

        runner.spawn(join=True)

        completed_files_exist = all(
            [isfile("completed/" + _file.split("/")[-1]) for _file in files]
        )
        txt_files_exist = all(
            [isfile(_file) for _file in ["jobs.txt", "finished_cleanly.txt"]]
        )
        dirs = ["completed", "input", "logs"]
        dirs_exist = all([isdir(_dir) for _dir in dirs])

        logs = glob.glob("logs/*.log")
        num_logs = len(logs)
        log_lines = []
        for log in logs:
            with open(log, "r") as f:
                log_lines.append(len(f.readlines()))

        self.assertEqual(num_logs, len(files), msg="Not enough log files!")
        self.assertTrue(
            all(lines > 5 for lines in log_lines), msg="Log files were too short!"
        )
        self.assertTrue(completed_files_exist)
        self.assertTrue(dirs_exist)
        self.assertTrue(txt_files_exist)

    def test_generic_batch_threaded(self):
        """Run a calculation in generic mode, check that res files
        are cycled over and not repeated by multiple threads. Multiple
        threads competing for the same job.

        """

        executable = "echo"
        files = glob.glob(REAL_PATH + "/data/structures/*.res")[0:2]
        for _file in files:
            shutil.copy(_file, ".")

        runner = BatchRun(
            seed="*.res",
            debug=False,
            mode="generic",
            verbosity=4,
            ncores=1,
            nprocesses=2,
            executable=executable,
        )

        runner.spawn(join=True)

        completed_files_exist = isfile("completed/" + _file.split("/")[-1])
        txt_files_exist = all(
            [isfile(_file) for _file in ["jobs.txt", "finished_cleanly.txt"]]
        )
        dirs = ["completed", "input", "logs"]
        dirs_exist = all([isdir(_dir) for _dir in dirs])

        logs = glob.glob("logs/*.log")
        num_logs = len(logs)
        log_lines = []
        for log in logs:
            with open(log, "r") as f:
                log_lines.append(len(f.readlines()))

        with open("jobs.txt", "r") as f:
            jobs_len = len(f.readlines())

        self.assertEqual(num_logs, len(files), msg="Not enough log files!")
        for f in files:
            self.assertTrue(os.path.isfile("input/{}".format(f.split("/")[-1])))
        self.assertTrue(
            all(lines > 5 for lines in log_lines), msg="Log files were too short!"
        )
        self.assertEqual(jobs_len, len(files))
        self.assertTrue(dirs_exist)
        self.assertTrue(completed_files_exist)
        self.assertTrue(txt_files_exist)

    def test_batch_nothing_todo(self):
        """Check that nothing is done when there's nothing to do..."""
        for file in glob.glob(REAL_PATH + "data/nothing_to_do/*.*"):
            shutil.copy(file, ".")

        runner = BatchRun(
            seed=["LiAs"],
            debug=False,
            no_reopt=True,
            verbosity=VERBOSITY,
            ncores=min(2, ACTUAL_NCORES),
            nprocesses=min(2, ACTUAL_NCORES // 2),
            executable=EXECUTABLE,
        )
        start = time.time()
        runner.spawn()
        elapsed = time.time() - start
        self.assertTrue(elapsed < 10, "Sluggish to quit!")

    def test_res_name_collision(self):
        """Check that run3 safely falls over if there is a file called <seed>.res."""
        shutil.copy(REAL_PATH + "data/file_collision/LiAs.cell", ".")
        shutil.copy(REAL_PATH + "data/file_collision/LiAs.res", ".")
        shutil.copy(REAL_PATH + "data/file_collision/LiAs.param", ".")

        failed_safely = False
        try:
            runner = BatchRun(
                seed=["LiAs"],
                debug=False,
                no_reopt=True,
                verbosity=VERBOSITY,
                ncores=min(2, ACTUAL_NCORES),
                nprocesses=min(2, ACTUAL_NCORES // 2),
                executable=EXECUTABLE,
            )
            runner.spawn()
        except InputError:
            failed_safely = True

        self.assertTrue(failed_safely)

    def test_missing_basics(self):
        """ " Check that run3 falls over when e.g. xc_functional is missing."""
        for file in glob.glob(REAL_PATH + "data/misisng_basics/*.*"):
            shutil.copy(file, ".")

        tests = ["missing_cutoff", "missing_kpts", "missing_xc", "missing_pspot"]
        errors = [False for test in tests]
        for ind, test in enumerate(tests):
            try:
                runner = BatchRun(
                    seed=["LiC_" + test],
                    debug=False,
                    no_reopt=True,
                    verbosity=VERBOSITY,
                    executable=EXECUTABLE,
                )
                runner.spawn()
            except InputError:
                errors[ind] = True

        self.assertTrue(all(errors))


class BenchmarkCastep(MatadorUnitTest):
    """Run some short CASTEP calculations and compare the timings
    to single core & multicore references.

    """

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_benchmark_dual_core_scf(self):
        """Test the time taken to perform a set number of SCF steps
        on 2 cores. CASTEP prints no total timing data for single core jobs.

        """
        from os import makedirs

        seed = "_LiC.res"
        shutil.copy(REAL_PATH + "data/structures/LiC.res", "_LiC.res")

        cell_dict, s = cell2dict(
            REAL_PATH + "/data/benchmark/LiC_scf/LiC_scf.cell",
            verbosity=VERBOSITY,
            db=False,
        )
        self.assertTrue(s)
        param_dict, s = param2dict(
            REAL_PATH + "/data/benchmark/LiC_scf/LiC_scf.param",
            verbosity=VERBOSITY,
            db=False,
        )
        self.assertTrue(s)

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/C_00PBE.usp", ".")
        with self.assertRaises(CalculationError):
            ComputeTask(
                ncores=min(2, ACTUAL_NCORES),
                nnodes=None,
                node=None,
                res=seed,
                param_dict=param_dict,
                cell_dict=cell_dict,
                verbosity=0,
                executable=EXECUTABLE,
                start=True,
            )

        outputs_exist = [
            isfile("bad_castep/_LiC.res"),
            isfile("bad_castep/_LiC.castep"),
        ]

        results, s = castep2dict("bad_castep/_LiC.castep", db=False)
        makedirs(REAL_PATH + "/data/benchmark/results", exist_ok=True)
        shutil.copy(
            "bad_castep/_LiC.castep",
            REAL_PATH
            + "/data/benchmark/results/_LiC_2core_castep{}.castep".format(
                results.get("castep_version", "xxx")
            ),
        )

        self.assertTrue(all(outputs_exist), "couldn't find output files!")
        self.assertTrue(s, "couldn't read output files!")
        self.assertLess(results["_time_estimated"], 8)

    @unittest.skipIf(
        (not CASTEP_PRESENT or not MPI_PRESENT),
        "castep or mpirun executable not found in PATH",
    )
    def test_benchmark_manycore_scf(self):
        """Test the time taken to perform a set number of SCF steps
        on many cores.

        """
        from os import makedirs

        seed = "_LiC.res"
        shutil.copy(REAL_PATH + "data/structures/LiC.res", "_LiC.res")

        cell_dict, s = cell2dict(
            REAL_PATH + "/data/benchmark/LiC_scf/LiC_scf.cell",
            verbosity=VERBOSITY,
            db=False,
        )
        self.assertTrue(s)
        param_dict, s = param2dict(
            REAL_PATH + "/data/benchmark/LiC_scf/LiC_scf.param",
            verbosity=VERBOSITY,
            db=False,
        )
        self.assertTrue(s)

        shutil.copy(REAL_PATH + "data/pspots/Li_00PBE.usp", ".")
        shutil.copy(REAL_PATH + "data/pspots/C_00PBE.usp", ".")
        with self.assertRaises(CalculationError):
            ComputeTask(
                ncores=NCORES,
                nnodes=None,
                node=None,
                res=seed,
                param_dict=param_dict,
                cell_dict=cell_dict,
                verbosity=0,
                executable=EXECUTABLE,
                start=True,
            )

        outputs_exist = [
            isfile("bad_castep/_LiC.res"),
            isfile("bad_castep/_LiC.castep"),
        ]

        results, s = castep2dict("bad_castep/_LiC.castep", db=False)
        makedirs(REAL_PATH + "/data/benchmark/results", exist_ok=True)
        shutil.copy(
            "bad_castep/_LiC.castep",
            REAL_PATH
            + "/data/benchmark/results/_LiC_{}core_castep{}.castep".format(
                results.get("num_mpi_processes", 0),
                results.get("castep_version", "xxx"),
            ),
        )

        self.assertTrue(all(outputs_exist), "couldn't find output files!")
        self.assertTrue(s, "couldn't read output files!")
        print(results["_time_estimated"])
        benchmark_data = {2: 2 * 7.4, 4: 4 * 4.0, 12: 16.8, 14: 22.4, 18: 23.4}
        warnings.warn(
            RuntimeWarning(
                "Run took {} s with {} MPI processes, with cumulative CPU time of {:.1f} s. Benchmark data\n = {}".format(
                    results["_time_estimated"],
                    results["num_mpi_processes"],
                    results["_time_estimated"] * results["num_mpi_processes"],
                    benchmark_data,
                )
            )
        )


if __name__ == "__main__":
    unittest.main()
