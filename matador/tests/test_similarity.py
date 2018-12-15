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

        uniq_inds, dupe_dict, _, _ = get_uniq_cursor(test_docs,
                                                     **{'dr': 0.1, 'gaussian_width': 0.1})
        self.assertEqual(uniq_inds, {6})

    def testUniqFilterWithHierarchy(self):
        import glob
        files = glob.glob(REAL_PATH + 'data/uniqueness_hierarchy/*.res')
        cursor = [res2dict(f)[0] for f in files]
        cursor = sorted(cursor, key=lambda x: x['enthalpy_per_atom'])[0:10]
        uniq_inds, _, _, _ = get_uniq_cursor(cursor, sim_tol=0.08, energy_tol=0.05, projected=True,
                                             **{'dr': 0.01, 'gaussian_width': 0.1})
        filtered_cursor = [cursor[ind] for ind in uniq_inds]
        self.assertEqual(len(uniq_inds), 2)
        self.assertEqual(len(filtered_cursor), 2)
        self.assertTrue('KP-NaP-OQMD_2817-CollCode14009' in filtered_cursor[0]['source'][0])
        self.assertTrue('KP-NaP-CollCode421420' in filtered_cursor[1]['source'][0])

    def testNoUniquenessRetainsAllStructures(self):
        import glob
        files = glob.glob(REAL_PATH + 'data/uniqueness_hierarchy/*.res')
        cursor = [res2dict(f)[0] for f in files]
        uniq_inds, _, _, _ = get_uniq_cursor(cursor, sim_tol=0, energy_tol=1e20, projected=True,
                                             **{'dr': 0.1, 'gaussian_width': 0.1})
        filtered_cursor = [cursor[ind] for ind in uniq_inds]
        self.assertEqual(len(filtered_cursor), len(cursor))


if __name__ == '__main__':
    unittest.main()
