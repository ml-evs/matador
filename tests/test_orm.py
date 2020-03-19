#!/usr/bin/env python
# standard library
import unittest
from os.path import realpath

# matador modules
from matador.orm.orm import DataContainer

# grab abs path for accessing test data
REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"


class ORMTest(unittest.TestCase):
    def test_data_container(self):
        data_dict = {"a": "b", "lattice_abc": [[10, 10, 10], [90, 90, 90]]}
        container = DataContainer(data_dict)
        self.assertEqual(container["a"], "b")
        data_dict["a"] = 2
        self.assertEqual(container["a"], "b")

        self.assertEqual(container["lattice_abc"], [[10, 10, 10], [90, 90, 90]])
        data_dict["lattice_abc"][0][0] = 100
        self.assertEqual(container["lattice_abc"], [[10, 10, 10], [90, 90, 90]])

        with self.assertRaises(AttributeError):
            container["lattice_abc"] = [[1, 2, 3], [4, 5, 6]]

        container["extra_data"] = [1, 2, 3]
