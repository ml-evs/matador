#!/usr/bin/env python
import unittest
from matador.similarity.similarity import get_uniq_cursor
from matador.utils.cell_utils import cart2volume, abc2cart
from matador.scrapers.castep_scrapers import res2dict
from os.path import realpath

REAL_PATH = '/'.join(realpath(__file__).split('/')[:-1]) + '/'


class SimilarityFilterTest(unittest.TestCase):
    """ Test similarity filter. """
    def testICSDPriority(self):
        test_docs = []
        i = 0
        while i < 10:
            test_doc, success = res2dict(REAL_PATH + 'data/KP_primitive.res', db=False)
            test_doc['text_id'] = ['primitive', 'cell']
            test_doc['lattice_cart'] = abc2cart(test_doc['lattice_abc'])
            test_doc['cell_volume'] = cart2volume(test_doc['lattice_cart'])
            test_doc['enthalpy_per_atom'] = 0
            test_docs.append(test_doc)
            i += 1

        uniq_inds, dupe_dict, _, _ = get_uniq_cursor(test_docs)
        self.assertEqual(uniq_inds, {0})

        test_docs[6]['source'] = ['KP-CollCode999999.res']
        test_docs[6]['icsd'] = 999999
        test_docs[6]['text_id'] = ['keep', 'this']

        uniq_inds, dupe_dict, _, _ = get_uniq_cursor(test_docs)
        self.assertEqual(uniq_inds, {6})


if __name__ == '__main__':
    unittest.main()
