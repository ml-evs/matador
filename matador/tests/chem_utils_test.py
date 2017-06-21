#!/usr/bin/env python
import unittest
from matador.utils.chem_utils import get_concentration
from matador.utils.chem_utils import get_generic_grav_capacity
from matador.utils.chem_utils import get_stoich, get_formula_from_stoich


class ChemUtilsTest(unittest.TestCase):
    """ Test chem utils functionality. """
    def testGravimetricCapacity(self):
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
        self.assertAlmostEqual(Q[0], 2596.09660218)
        self.assertAlmostEqual(Q[2], 167.449398573)
        self.assertEqual(Q[0], Q[1])
        self.assertEqual(round(8*Q[2], 3), round(Q[3], 3))
        self.assertEqual(round(Q[2], 3), round(2*Q[4], 3))

    def testAtoms2Stoich(self):
        atoms = 5*['Li']
        atoms.extend(5*['P'])
        stoich = [['Li', 1], ['P', 1]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 99*['Li']
        atoms.extend(1*['P'])
        stoich = [['Li', 99], ['P', 1]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 4*['Li']
        atoms.extend(36*['P'])
        stoich = [['Li', 1], ['P', 9]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 3*['Li']
        atoms.extend(2*['P'])
        stoich = [['Li', 3], ['P', 2]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 9*['Li']
        atoms.extend(6*['P'])
        stoich = [['Li', 3], ['P', 2]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 36*['P']
        atoms.extend(4*['Li'])
        stoich = [['Li', 1], ['P', 9]]
        self.assertEqual(stoich, get_stoich(atoms))

    def testStoich2Form(self):
        stoich = [['Li', 1], ['P', 9]]
        form = 'LiP9'
        self.assertEqual(form, get_formula_from_stoich(stoich))
        stoich = [['P', 9], ['Li', 1]]
        form = 'LiP9'
        self.assertEqual(form, get_formula_from_stoich(stoich))
        stoich = [['Li', 1], ['P', 9]]
        form = 'LiP\\textsubscript{9}'
        self.assertEqual(form, get_formula_from_stoich(stoich, tex=True))
        stoich = [['Li', 1], ['P', 9]]
        form = 'P9Li'
        self.assertEqual(form, get_formula_from_stoich(stoich, elements=['P', 'Li']))


if __name__ == '__main__':
    unittest.main()
