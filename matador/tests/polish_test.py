#!/usr/bin/env python
import unittest
from matador.polish import Polisher
from matador.utils.chem_utils import get_periodic_table
import sys
import os


class SwapTest(unittest.TestCase):
    """ Test atomic swap functions. """
    def testOneToManySwap(self):
        # spoof Polisher __init__
        swap_args = {'swap': ['As[P,Sb,Zn,Cu]'], 'debug': False}
        bare_polish = Polisher.__new__(Polisher)
        bare_polish.periodic_table = get_periodic_table()
        bare_polish.args = swap_args
        # try to parse swaps
        with open(os.devnull, 'w') as sys.stdout:
            bare_polish.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_polish.swap_pairs, [[['As'], ['P', 'Sb', 'Zn', 'Cu']]])
        # set up test data for real swap
        doc = dict()
        doc['atom_types'] = ['Li', 'Li', 'As', 'As']
        swapped_docs, num_swapped = bare_polish.atomic_swaps(doc)
        self.assertEqual(num_swapped, 4)
        P_found = False
        Sb_found = False
        Zn_found = False
        Cu_found = False
        for new_doc in swapped_docs:
            self.assertTrue('As' not in new_doc['atom_types'])
            if 'P' in new_doc['atom_types']:
                self.assertTrue(x not in new_doc['atom_types'] for x in ['Sb', 'Zn', 'Cu'])
                self.assertEqual(new_doc['atom_types'], ['Li', 'Li', 'P', 'P'])
                P_found = True
            if 'Sb' in new_doc['atom_types']:
                self.assertTrue(x not in new_doc['atom_types'] for x in ['P', 'Zn', 'Cu'])
                self.assertEqual(new_doc['atom_types'], ['Li', 'Li', 'Sb', 'Sb'])
                Sb_found = True
            if 'Zn' in new_doc['atom_types']:
                self.assertTrue(x not in new_doc['atom_types'] for x in ['P', 'Sb', 'Cu'])
                self.assertEqual(new_doc['atom_types'], ['Li', 'Li', 'Zn', 'Zn'])
                Zn_found = True
            if 'Cu' in new_doc['atom_types']:
                self.assertTrue(x not in new_doc['atom_types'] for x in ['P', 'Sb', 'Zn'])
                self.assertEqual(new_doc['atom_types'], ['Li', 'Li', 'Cu', 'Cu'])
                Cu_found = True
        self.assertTrue(P_found)
        self.assertTrue(Sb_found)
        self.assertTrue(Zn_found)
        self.assertTrue(Cu_found)

    def testManyToOneMacroSwap(self):
        swap_args = {'swap': ['[V]P'], 'debug': False}
        bare_polish = Polisher.__new__(Polisher)
        bare_polish.periodic_table = get_periodic_table()
        bare_polish.args = swap_args
        # try to parse swaps
        with open(os.devnull, 'w') as sys.stdout:
            bare_polish.parse_swaps()
        sys.stdout = sys.__stdout__
        self.assertEqual(bare_polish.swap_pairs, [[['N', 'P', 'As', 'Sb', 'Bi'], ['P']]])
        # set up test data for real swap
        doc = dict()
        doc['atom_types'] = ['P', 'Sb', 'As', 'As']
        swapped_docs, num_swapped = bare_polish.atomic_swaps(doc)
        self.assertEqual(num_swapped, 1)
        self.assertEqual(swapped_docs[0]['atom_types'], ['P', 'P', 'P', 'P'])

if __name__ == '__main__':
    unittest.main()
