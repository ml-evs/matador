#!/usr/bin/env python

import unittest
import numpy as np
import json


class MiscUtilTest(unittest.TestCase):
    """Tests cursor util functions."""

    def test_serialize_numpy(self):
        from matador.utils.print_utils import dumps

        doc = {"a": "b", "1": [1, 2, 3], "2": np.zeros((3, 3))}
        doc_list = {"a": "b", "1": [1, 2, 3], "2": [[0, 0, 0], [0, 0, 0], [0, 0, 0]]}
        string = dumps(doc)
        self.assertEqual(json.loads(string), doc_list)


if __name__ == "__main__":
    unittest.main()
