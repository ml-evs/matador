#!/usr/bin/env python
import unittest
from matador.polish import Polisher
from matador.utils.chem_utils import get_periodic_table


class SwapTest(unittest.TestCase):
    """ Test swap functions. """
    def testOneToManySwap(self):
        """ Test one-to-many swaps. """
        print('\n', 20*'-', sep='')
        # spoof Polisher __init__
        swap_args = {'swap': ['As[P,Sb]'], 'debug': True}
        bare_polish = Polisher.__new__(Polisher)
        bare_polish.periodic_table = get_periodic_table()
        bare_polish.args = swap_args
        # try to parse swaps
        bare_polish.parse_swaps()
        self.assertEqual(bare_polish.swap_pairs, [[['As'], ['P', 'Sb']]])
        # set up test data for real swap
        doc = dict()
        doc['atom_types'] = ['Li', 'Li', 'As', 'As']
        swapped_docs, num_swapped = bare_polish.atomic_swaps(doc)
        self.assertEqual(num_swapped, 2)
        for new_doc in swapped_docs:
            self.assertTrue('As' not in new_doc['atom_types'])
            if 'P' in new_doc['atom_types']:
                self.assertTrue('Sb' not in new_doc['atom_types'])
            if 'Sb' in new_doc['atom_types']:
                self.assertTrue('Sb' not in new_doc['atom_types'])

    def testMacroSwap(self):
        """ Test atomic swaps from a group macro. """
        print('\n', 20*'-', sep='')
        swap_args = {'swap': ['[V]P'], 'debug': True}
        bare_polish = Polisher.__new__(Polisher)
        bare_polish.periodic_table = get_periodic_table()
        bare_polish.args = swap_args
        # try to parse swaps
        bare_polish.parse_swaps()
        self.assertEqual(bare_polish.swap_pairs, [[['N', 'P', 'As', 'Sb', 'Bi'], ['P']]])
        # set up test data for real swap
        doc = dict()
        doc['atom_types'] = ['P', 'Sb', 'As', 'As']
        swapped_docs, num_swapped = bare_polish.atomic_swaps(doc)
        self.assertEqual(num_swapped, 1)
        self.assertEqual(swapped_docs[0]['atom_types'], ['P', 'P', 'P', 'P'])

if __name__ == '__main__':
    unittest.main()
