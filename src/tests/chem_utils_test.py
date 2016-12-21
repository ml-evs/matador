#!/usr/bin/env python
import unittest
from utils.chem_utils import get_concentration, get_num_intercalated
from utils.chem_utils import get_binary_grav_capacities, get_molar_mass
from utils.chem_utils import get_generic_grav_capacity


class CapacityTest(unittest.TestCase):
    """ Test capacity functions. """
    def test(self):
        print('Testing capacity...')
        test_docs = []
        test_elements = []
        Q = []

        # Li3P and symmetry test
        doc = dict()
        doc['stoichiometry'] = [['Li', 3], ['P', 1]]
        test_docs.append(doc)
        test_elements.append(['Li', 'P'])
        doc = dict()
        doc['stoichiometry'] = [['P', 1], ['Li', 3]]
        test_docs.append(doc)
        test_elements.append(['Li', 'P'])
        # ternary test
        doc = dict()
        doc['stoichiometry'] = [['Li', 1], ['Mo', 1], ['S', 2]]
        test_docs.append(doc)
        test_elements.append(['Li', 'Mo', 'S'])
        doc = dict()
        doc['stoichiometry'] = [['Li', 8], ['Mo', 1], ['S', 2]]
        test_docs.append(doc)
        test_elements.append(['Li', 'Mo', 'S'])
        doc = dict()
        doc['stoichiometry'] = [['Li', 1], ['Mo', 2], ['S', 4]]
        test_docs.append(doc)
        test_elements.append(['Li', 'Mo', 'S'])
        for doc, elem in zip(test_docs, test_elements):
            doc['concentration'] = get_concentration(doc, elem)
            temp_conc = list(doc['concentration'])
            temp_conc.append(1.0)
            for conc in doc['concentration']:
                temp_conc[-1] -= conc

            Q.append(get_generic_grav_capacity(temp_conc, elem))
            print(doc['stoichiometry'], Q[-1], 'mAh/g')

        assert(Q[0] == Q[1])
        assert(round(8*Q[2], 3) == round(Q[3], 3))
        assert(round(Q[2], 3) == round(2*Q[4], 3))

if __name__ == '__main__':
    unittest.main()
