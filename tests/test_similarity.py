#!/usr/bin/env python
import unittest
from os.path import realpath
import numpy as np
from matador.fingerprints.similarity import get_uniq_cursor
from matador.fingerprints import Fingerprint
from matador.scrapers.castep_scrapers import res2dict

REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"


class SimilarityFilterTest(unittest.TestCase):
    """Test similarity filter."""

    def test_icsd_priority(self):
        test_docs = []
        i = 0
        while i < 10:
            test_doc, _ = res2dict(REAL_PATH + "data/KP_primitive.res", db=False)
            test_doc["text_id"] = ["primitive", "cell"]
            test_docs.append(test_doc)
            i += 1

        uniq_inds, _, _, _ = get_uniq_cursor(test_docs)
        self.assertEqual(uniq_inds, [0])

        test_docs[6]["source"] = ["KP-CollCode999999.res"]
        test_docs[6]["icsd"] = 999999
        test_docs[6]["text_id"] = ["keep", "this"]

        uniq_inds, _, _, _ = get_uniq_cursor(
            test_docs, **{"dr": 0.1, "gaussian_width": 0.1}
        )
        self.assertEqual(uniq_inds, [6])

    def test_k3p_uniq_default(self):
        cursor, _ = res2dict(REAL_PATH + "data/K3P_uniq/*.res")
        cursor = sorted(cursor, key=lambda x: x["enthalpy_per_atom"])
        uniq_inds, _, _, _ = get_uniq_cursor(cursor)
        filtered_cursor = [cursor[ind] for ind in uniq_inds]
        self.assertEqual(len(cursor), 10)
        self.assertEqual(len(filtered_cursor), 5)
        found = []
        correct_structures = [
            "K3P-OQMD_4786-CollCode25550",
            "K3P-mode-follow-swap-Na3N-OQMD_21100-CollCode165992",
            "KP-fvsqdf",
            "PK-NNa3-OQMD_21100-CollCode165992",
            "KP-yzcni8",
        ]
        for struct in correct_structures:
            for doc in filtered_cursor:
                if struct in doc["source"][0]:
                    found.append(True)
                    break
            else:
                found.append(False)

        if not all(found):
            print([doc["source"][0] for doc in filtered_cursor])

        self.assertTrue(all(found))

    def test_volume_rescale(self):
        import numpy as np

        test_doc, success = res2dict(REAL_PATH + "data/KP_primitive.res", db=False)
        self.assertTrue(success)
        test_docs = []
        rescale = np.linspace(0.1, 10, 8)
        lattice = np.asarray(test_doc["lattice_abc"])
        for val in rescale:
            test_docs.append(test_doc)
            test_docs[-1]["lattice_abc"] = lattice
            test_docs[-1]["lattice_abc"][0] *= val
            test_docs[-1]["lattice_abc"] = test_docs[-1]["lattice_abc"].tolist()
        uniq_inds, _, _, _ = get_uniq_cursor(test_docs)
        self.assertEqual(uniq_inds, [0])

    def test_uniq_filter_with_hierarchy(self):
        import glob

        files = glob.glob(REAL_PATH + "data/uniqueness_hierarchy/*.res")
        cursor = [res2dict(f)[0] for f in files]
        cursor = sorted(cursor, key=lambda x: x["enthalpy_per_atom"])[0:10]
        uniq_inds, _, _, _ = get_uniq_cursor(
            cursor,
            sim_tol=0.08,
            energy_tol=0.05,
            projected=True,
            **{"dr": 0.01, "gaussian_width": 0.1}
        )
        filtered_cursor = [cursor[ind] for ind in uniq_inds]
        self.assertEqual(len(uniq_inds), 2)
        self.assertEqual(len(filtered_cursor), 2)
        self.assertTrue(
            "KP-NaP-OQMD_2817-CollCode14009" in filtered_cursor[0]["source"][0]
        )
        self.assertTrue("KP-NaP-CollCode421420" in filtered_cursor[1]["source"][0])

    def test_uniq_filter_with_hierarchy_2(self):
        cursor, f_ = res2dict(REAL_PATH + "data/hull-LLZO/*LLZO*.res")
        cursor = sorted(cursor, key=lambda x: x["enthalpy_per_atom"])[0:10]
        uniq_inds, _, _, _ = get_uniq_cursor(
            cursor,
            sim_tol=0.1,
            energy_tol=1e10,
            projected=True,
            **{"dr": 0.01, "gaussian_width": 0.1}
        )
        filtered_cursor = [cursor[ind] for ind in uniq_inds]
        self.assertEqual(len(uniq_inds), 1)
        self.assertEqual(len(filtered_cursor), 1)
        self.assertTrue("cubic-LLZO-CollCode999999" in filtered_cursor[0]["source"][0])

    def test_double_uniqueness_hierarchy(self):
        import glob

        files = glob.glob(REAL_PATH + "data/uniqueness_hierarchy/*.res")
        files += glob.glob(REAL_PATH + "data/hull-LLZO/*LLZO*.res")
        cursor = [res2dict(f)[0] for f in files]
        cursor = sorted(cursor, key=lambda x: x["enthalpy_per_atom"])[0:10]
        uniq_inds, _, _, _ = get_uniq_cursor(
            cursor,
            sim_tol=0.1,
            energy_tol=1e20,
            projected=True,
            **{"dr": 0.01, "gaussian_width": 0.1}
        )
        filtered_cursor = [cursor[ind] for ind in uniq_inds]
        self.assertEqual(len(uniq_inds), 3)
        self.assertEqual(len(filtered_cursor), 3)
        print([doc["source"] for doc in filtered_cursor])
        self.assertTrue("cubic-LLZO-CollCode999999" in filtered_cursor[0]["source"][0])
        self.assertTrue(
            "KP-NaP-OQMD_2817-CollCode14009" in filtered_cursor[1]["source"][0]
        )
        self.assertTrue("KP-NaP-CollCode421420" in filtered_cursor[2]["source"][0])

    def test_no_overlap_retains_all_structures(self):
        import glob

        files = glob.glob(REAL_PATH + "data/uniqueness_hierarchy/*.res")
        cursor = [res2dict(f)[0] for f in files]
        uniq_inds, _, _, _ = get_uniq_cursor(
            cursor,
            sim_tol=0,
            energy_tol=1e20,
            projected=True,
            debug=True,
            **{"dr": 0.1, "gaussian_width": 0.1}
        )
        filtered_cursor = [cursor[ind] for ind in uniq_inds]
        self.assertEqual(len(filtered_cursor), len(cursor))

    def test_with_crystals(self):
        from matador.crystal import Crystal
        import glob

        files = glob.glob(REAL_PATH + "data/uniqueness_hierarchy/*.res")
        cursor = [Crystal(res2dict(f)[0]) for f in files]
        uniq_inds, _, _, _ = get_uniq_cursor(
            cursor,
            sim_tol=0,
            energy_tol=1e20,
            projected=True,
            debug=True,
            **{"dr": 0.1, "gaussian_width": 0.1}
        )
        filtered_cursor = [cursor[ind] for ind in uniq_inds]
        self.assertEqual(len(filtered_cursor), len(cursor))

    def test_with_skips(self):
        from matador.crystal import Crystal
        from matador.utils.cursor_utils import filter_unique_structures
        import glob

        files = glob.glob(REAL_PATH + "data/uniqueness_hierarchy/*.res")
        cursor = [Crystal(res2dict(f)[0]) for f in files]
        filtered_cursor = filter_unique_structures(cursor, energy_tol=0)
        self.assertEqual(len(filtered_cursor), len(cursor))

        cursor = sorted(
            [res2dict(f)[0] for f in files], key=lambda doc: doc["enthalpy_per_atom"]
        )[0:10]
        for ind, doc in enumerate(cursor):
            doc["enthalpy_per_atom"] = float(-ind)

        cursor[8]["enthalpy_per_atom"] = -5.0
        cursor[9]["enthalpy_per_atom"] = -5.0001

        filtered_cursor = filter_unique_structures(cursor, energy_tol=0.003)
        self.assertEqual(len(filtered_cursor), 8)


class TestBroadening(unittest.TestCase):
    def test_broadening_agreement(self):

        hist = np.zeros((1000))
        hist[100] = 10
        hist[800] = 20
        x = np.linspace(0, 10, 1000)
        dists = [x[100]] * 10 + [x[800]] * 20

        for btype in ["gaussian", "lorentzian"]:
            gr1 = Fingerprint._broadening_distance_dominated(hist, x, 1, btype)
            gr2 = Fingerprint._broadening_unrolled(hist, x, 1, btype)
            gr3 = Fingerprint._broadening_space_dominated(dists, x, 1, btype)

            np.testing.assert_almost_equal(gr1, gr2, decimal=10)
            np.testing.assert_almost_equal(gr2, gr3, decimal=10)
