#!/usr/bin/env python

import unittest

try:
    import ase  # noqa

    ASE_IMPORTED = True
except ImportError:
    ASE_IMPORTED = False

try:
    import pymatgen  # noqa

    PMG_IMPORTED = True
except ImportError:
    PMG_IMPORTED = False


@unittest.skipIf(~(PMG_IMPORTED and ASE_IMPORTED), "Unable to import ASE/pymtagen")
class PMGUtilTest(unittest.TestCase):
    """ Tests cursor util functions. """

    def test_pm2dict(self):
        pass


if __name__ == "__main__":
    unittest.main()
