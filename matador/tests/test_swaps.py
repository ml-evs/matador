#!/usr/bin/env python
import unittest
import sys
import os
from matador.swaps import AtomicSwapper
from matador.utils.chem_utils import get_periodic_table


class SwapTest(unittest.TestCase):
    """ Test atomic swap functions. """

    def testSingleSimpleSwap(self):
        # spoof AtomicSwapper __init__
        swap_args = {"swap": ["AsP"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        with open(os.devnull, "w") as sys.stdout:
            bare_swap.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_swap.swap_pairs, [[["As"], ["P"]]])
        # set up test data for real swap
        doc = dict()
        doc["atom_types"] = ["Li", "Li", "As", "As"]
        swapped_docs, num_swapped = bare_swap.atomic_swaps(doc)
        self.assertEqual(num_swapped, 1)
        self.assertEqual(swapped_docs[0]["atom_types"], ["Li", "Li", "P", "P"])

    def testNullSelfSwap(self):
        # spoof AtomicSwapper __init__
        swap_args = {"swap": ["KK:PP"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        with open(os.devnull, "w") as sys.stdout:
            bare_swap.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_swap.swap_pairs, [[["K"], ["K"]], [["P"], ["P"]]])
        # set up test data for real swap
        doc = dict()
        doc["atom_types"] = ["K", "K", "P", "P"]
        swapped_docs, num_swapped = bare_swap.atomic_swaps(doc)
        self.assertEqual(num_swapped, 0)

    def testMultipleSimpleSwap(self):
        # spoof AtomicSwapper __init__
        swap_args = {"swap": ["AsP:LiNa"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        with open(os.devnull, "w") as sys.stdout:
            bare_swap.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_swap.swap_pairs, [[["As"], ["P"]], [["Li"], ["Na"]]])
        # set up test data for real swap
        doc = dict()
        doc["atom_types"] = ["Li", "Li", "As", "As"]
        swapped_docs, num_swapped = bare_swap.atomic_swaps(doc)
        self.assertEqual(num_swapped, 1)
        self.assertEqual(swapped_docs[0]["atom_types"], ["Na", "Na", "P", "P"])

    def testOneToManySwap(self):
        # spoof AtomicSwapper __init__
        swap_args = {"swap": ["As[P,Sb,Zn,Cu]"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        # with open(os.devnull, 'w') as sys.stdout:
        bare_swap.parse_swaps()
        # sys.stdout = sys.__stdout__
        self.assertEqual(bare_swap.swap_pairs, [[["As"], ["P", "Sb", "Zn", "Cu"]]])
        # set up test data for real swap
        doc = dict()
        doc["atom_types"] = ["Li", "Li", "As", "As"]
        swapped_docs, num_swapped = bare_swap.atomic_swaps(doc)
        self.assertEqual(num_swapped, 4)
        P_found = False
        Sb_found = False
        Zn_found = False
        Cu_found = False
        for new_doc in swapped_docs:
            self.assertTrue("As" not in new_doc["atom_types"])
            if "P" in new_doc["atom_types"]:
                self.assertTrue(x not in new_doc["atom_types"] for x in ["Sb", "Zn", "Cu"])
                self.assertEqual(new_doc["atom_types"], ["Li", "Li", "P", "P"])
                P_found = True
            if "Sb" in new_doc["atom_types"]:
                self.assertTrue(x not in new_doc["atom_types"] for x in ["P", "Zn", "Cu"])
                self.assertEqual(new_doc["atom_types"], ["Li", "Li", "Sb", "Sb"])
                Sb_found = True
            if "Zn" in new_doc["atom_types"]:
                self.assertTrue(x not in new_doc["atom_types"] for x in ["P", "Sb", "Cu"])
                self.assertEqual(new_doc["atom_types"], ["Li", "Li", "Zn", "Zn"])
                Zn_found = True
            if "Cu" in new_doc["atom_types"]:
                self.assertTrue(x not in new_doc["atom_types"] for x in ["P", "Sb", "Zn"])
                self.assertEqual(new_doc["atom_types"], ["Li", "Li", "Cu", "Cu"])
                Cu_found = True
        self.assertTrue(P_found)
        self.assertTrue(Sb_found)
        self.assertTrue(Zn_found)
        self.assertTrue(Cu_found)

    def testMistakenMacro(self):
        swap_args = {"swap": ["VP"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        with open(os.devnull, "w") as sys.stdout:
            bare_swap.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_swap.swap_pairs, [[["V"], ["P"]]])

    def testManyToOneMacroSwap(self):
        swap_args = {"swap": ["[V]P"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        with open(os.devnull, "w") as sys.stdout:
            bare_swap.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_swap.swap_pairs, [[["N", "P", "As", "Sb", "Bi"], ["P"]]])
        # set up test data for real swap
        doc = dict()
        doc["atom_types"] = ["P", "Sb", "As", "As"]
        swapped_docs, num_swapped = bare_swap.atomic_swaps(doc)
        self.assertEqual(num_swapped, 1)
        self.assertEqual(swapped_docs[0]["atom_types"], ["P", "P", "P", "P"])

    def testManyToManyMacroSwap(self):
        swap_args = {"swap": ["[V][Tc,Mo]"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        with open(os.devnull, "w") as sys.stdout:
            bare_swap.parse_swaps()
        sys.stdout = sys.__stdout__
        # set up test data for real swap
        doc = dict()
        doc["atom_types"] = ["P", "Sb", "As", "As", "Bi"]
        self.assertEqual(bare_swap.swap_pairs, [[["N", "P", "As", "Sb", "Bi"], ["Tc", "Mo"]]])
        swapped_docs, num_swapped = bare_swap.atomic_swaps(doc)
        self.assertEqual(num_swapped, 2)
        self.assertEqual(swapped_docs[0]["atom_types"], ["Tc", "Tc", "Tc", "Tc", "Tc"])
        self.assertEqual(swapped_docs[1]["atom_types"], ["Mo", "Mo", "Mo", "Mo", "Mo"])

    def testMultipleManyToOneSwap(self):
        # spoof AtomicSwapper __init__
        swap_args = {"swap": ["[Li,Na]K:[Ru, Rh]La"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        with open(os.devnull, "w") as sys.stdout:
            bare_swap.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_swap.swap_pairs, [[["Li", "Na"], ["K"]], [["Ru", "Rh"], ["La"]]])
        # set up test data for real swap
        doc = dict()
        doc["atom_types"] = ["Li", "Na", "Ru", "Rh"]
        swapped_docs, num_swapped = bare_swap.atomic_swaps(doc)
        self.assertEqual(num_swapped, 1)
        self.assertEqual(swapped_docs[0]["atom_types"], ["K", "K", "La", "La"])

    def testMultipleManyToManySwapsWithProblematicElementNames(self):
        # spoof AtomicSwapper __init__
        swap_args = {"swap": ["[Li,Na]K:[V,I][I]"], "debug": True}
        bare_swap = AtomicSwapper.__new__(AtomicSwapper)
        bare_swap.periodic_table = get_periodic_table()
        bare_swap.args = swap_args
        # try to parse swaps
        with open(os.devnull, "w") as sys.stdout:
            bare_swap.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_swap.swap_pairs, [[["Li", "Na"], ["K"]], [["V", "I"], ["Li", "Na", "K", "Rb", "Cs", "Fr"]]])
        # set up test data for real swap
        doc = dict()
        doc["atom_types"] = ["Li", "Na", "V", "I"]
        swapped_docs, num_swapped = bare_swap.atomic_swaps(doc)
        self.assertEqual(num_swapped, 6)
        self.assertEqual(swapped_docs[0]["atom_types"], ["K", "K", "Li", "Li"])
        self.assertEqual(swapped_docs[5]["atom_types"], ["K", "K", "Fr", "Fr"])


if __name__ == "__main__":
    unittest.main()