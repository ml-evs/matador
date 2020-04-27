#!/usr/bin/env python
import unittest
from os.path import realpath

# grab abs path for accessing test data
REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"


class CursorUtilTest(unittest.TestCase):
    """ Tests cursor util functions. """

    def test_guess_prov(self):
        from matador.utils.cursor_utils import get_guess_doc_provenance

        sources = []
        answers = []

        sources.append("/my/stupid/long/path/KSnP-asdf123.castep")
        answers.append("AIRSS")

        sources.append(["/my/stupid/long/path/KSnP-asdf123.castep"])
        answers.append("AIRSS")

        sources.append(["/my/ICSD-CollCode678321/stupid/long/path/KSnP-asdf123.castep"])
        answers.append("AIRSS")

        sources.append(
            ["KSnP-asdf123.res", "KSnP-ICSD-CollCode123123.cell", "OQMD_11111.param"]
        )
        answers.append("AIRSS")

        sources.append(["Na-Collo.res"])
        answers.append("ICSD")

        sources.append(["OQMD_64572"])
        answers.append("OQMD")

        sources.append(["PK-BiK-OQMD_647548-CollCode55065.castep"])
        answers.append("SWAPS")

        sources.append(["KSnP-GA-abcdef-5x101.res"])
        answers.append("GA")

        sources.append(["/u/fs1/swaps+known/LiP-KSnP-GA-abcdef-5x101.res"])
        answers.append("SWAPS")

        sources.append(
            ["/u/fs1/swaps+known/LiP-CollCode10101-swap-KSnP-GA-abcdef-5x101.res"]
        )
        answers.append("SWAPS")

        sources.append(["Ag2Bi2I8-MP-35909_300eV.castep"])
        answers.append("MP")

        sources.append(
            ["AgBiI4-spinel-Config5-DOI-10.17638__datacat.liverpool.ac.uk__240.castep"]
        )
        answers.append("DOI")

        for source, answer in zip(sources, answers):
            self.assertEqual(
                get_guess_doc_provenance(source), answer, msg="failed {}".format(source)
            )

    def test_filter_cursor(self):
        from matador.utils.cursor_utils import filter_cursor

        cursor = [{"field": int} for int in range(100)]

        filtered = filter_cursor(cursor, "field", 98)
        self.assertEqual(len(filtered), 2)

        filtered = filter_cursor(cursor, "field", [-2, 100])
        self.assertEqual(len(filtered), len(cursor))

        cursor.append({"no_field": 1})
        filtered = filter_cursor(cursor, "field", [-2, 100])
        self.assertEqual(len(filtered), len(cursor) - 1)

        cursor.append({"field": [0]})
        filtered = filter_cursor(cursor, "field", [-2, 100])
        self.assertEqual(len(filtered), len(cursor) - 2)

    def test_recursive_get_set(self):
        from matador.utils.cursor_utils import recursive_get, recursive_set

        nested_dict = {
            "source": ["blah", "foo"],
            "lattice_cart": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "_beef": {
                "total_energy": [1, 2, 3, 4],
                "thetas": [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                "foo": {"blah": "bloop"},
            },
        }
        self.assertEqual(recursive_get(nested_dict, ["_beef", "thetas", -1]), [7, 8, 9])
        self.assertEqual(recursive_get(nested_dict, ["_beef", "total_energy", 2]), 3)
        self.assertEqual(recursive_get(nested_dict, ["_beef", "foo", "blah"]), "bloop")
        recursive_set(nested_dict, ["_beef", "thetas", -1], [1, 2, 3])
        self.assertEqual(recursive_get(nested_dict, ["_beef", "thetas", -1]), [1, 2, 3])
        recursive_set(nested_dict, ["_beef", "total_energy", 2], 3.5)
        self.assertEqual(recursive_get(nested_dict, ["_beef", "total_energy", 2]), 3.5)
        nested_dict["_beef"]["total_energy"][2] = 3
        self.assertEqual(recursive_get(nested_dict, ["_beef", "total_energy", 2]), 3)
        nested_dict["_beef"]["foo"]["blah"] = (1, 2)
        self.assertEqual(recursive_get(nested_dict, ["_beef", "foo", "blah"]), (1, 2))

        with self.assertRaises(IndexError):
            recursive_get(nested_dict, ["_beef", "thetas", 4])
        with self.assertRaises(KeyError):
            recursive_get(nested_dict, ["_beef", "thetaz", 4])
        with self.assertRaises(KeyError):
            recursive_get(nested_dict, ["_beef", "foo", "blahp"])


if __name__ == "__main__":
    unittest.main()
