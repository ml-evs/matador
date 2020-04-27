""" Useful base classes for tests. """

import unittest
import os
import subprocess as sp
import shutil

REAL_PATH = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"
TMP_DIR = "tmp_test"
ROOT_DIR = os.getcwd()


def detect_program(executable):
    try:
        with open("/dev/null", "w") as devnull:
            out, errs = sp.Popen(
                [executable, "--version"], stdout=devnull, stderr=devnull
            ).communicate()
        if errs:
            raise RuntimeError
        return True
    except Exception:
        return False


class MatadorUnitTest(unittest.TestCase):
    def tearDown(self):
        os.chdir(REAL_PATH)
        if os.path.isdir(TMP_DIR):
            shutil.rmtree(TMP_DIR)

        os.chdir(ROOT_DIR)

    def setUp(self):
        os.chdir(REAL_PATH)
        os.makedirs(TMP_DIR, exist_ok=True)
        os.chdir(TMP_DIR)
