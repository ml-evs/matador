#!/usr/bin/env python
import unittest
from matador.query import DBQuery, parse_element_string
import json


class QueryTest(unittest.TestCase):
    """ Test query functionality. """
    def testBasicQueries(self):
        kwargs = {'composition': 'KP', 'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'atom_types': {'$in': ['K']}},
                    {'atom_types': {'$in': ['P']}},
                    {'stoichiometry': {'$size': 2}}
                    ]},
                {'$or': [
                    {'quality': {'$gt': 0}},
                    {'quality': {'$exists': False}}
                    ]},
                ]
        })
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'formula': 'K3P4', 'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'stoichiometry': {'$in': [['K', 3.0]]}},
                    {'stoichiometry': {'$in': [['P', 4.0]]}},
                    {'stoichiometry': {'$size': 2}}
                    ]},
                {'$or': [
                    {'quality': {'$gt': 0}},
                    {'quality': {'$exists': False}}
                    ]},
                ]
        })
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'formula': 'K3P6', 'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'stoichiometry': {'$in': [['K', 1.0]]}},
                    {'stoichiometry': {'$in': [['P', 2.0]]}},
                    {'stoichiometry': {'$size': 2}}
                    ]},
                {'$or': [
                    {'quality': {'$gt': 0}},
                    {'quality': {'$exists': False}}
                    ]},
                ]
        })
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'formula': 'K27P9', 'partial_formula': True, 'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'stoichiometry': {'$in': [['K', 3]]}},
                    {'stoichiometry': {'$in': [['P', 1]]}}
                    ]},
                {'$or': [
                    {'quality': {'$gt': 0}},
                    {'quality': {'$exists': False}}
                    ]},
                ]
        })
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'formula': ['K27P9'], 'partial_formula': True, 'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'stoichiometry': {'$in': [['K', 3]]}},
                    {'stoichiometry': {'$in': [['P', 1]]}}
                    ]},
                {'$or': [
                    {'quality': {'$gt': 0}},
                    {'quality': {'$exists': False}}
                    ]},
                ]
        })
        self.assertDictEqual(test_dict, query.query_dict)

    def testComplexQueries(self):
        """ Test long queries with multiple mismatching of lists
        (emulating argparse) and values.
        """
        self.maxDiff = None
        kwargs = {'formula': 'K27P9', 'space_group': ['Fd-3m'], 'num_species': 3,
                  'num_fu': [4], 'tags': ['foo', 'Bar'], 'doi': ['1001/4001'],
                  'icsd': 100020, 'cutoff': 300, 'spin': 1,
                  'partial_formula': True, 'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'stoichiometry': {'$in': [['K', 3]]}},
                    {'stoichiometry': {'$in': [['P', 1]]}},
                    ]},
                {'stoichiometry': {'$size': 3}},
                {'space_group': 'Fd-3m'},
                {'num_fu': {'$gte': 4}},
                {'$and': [
                    {'tags': {'$in': ['foo']}},
                    {'tags': {'$in': ['Bar']}}
                ]},
                {'doi': {'$in': ['1001/4001']}},
                {'icsd': {'$eq': '100020'}},
                {'cut_off_energy': {'$gte': 300}},
                {'spin_polarized': True},
                {'$or': [
                    {'quality': {'$gt': 0}},
                    {'quality': {'$exists': False}}
                    ]},
                ]
        })
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': 'LiFeBe', 'icsd': 0, 'ignore_warnings': True,
                  'src_str': '/Foo/bar/foo/Bar.res', 'pressure': 5,
                  'cutoff': [300, 400], 'encapsulated': True, 'cnt_radius': 5.21,
                  'sedc': 'null', 'mp_spacing': [0.05], 'spin': 0,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'atom_types': {'$in': ['Li']}},
                    {'atom_types': {'$in': ['Fe']}},
                    {'atom_types': {'$in': ['Be']}},
                    {'stoichiometry': {'$size': 3}}
                    ]},
                {'icsd': {'$exists': True}},
                {'cut_off_energy': {'$gte': 300, '$lte': 400}},
                {'source': {'$in': ['/Foo/bar/foo/Bar.res']}},
                {'pressure': {'$lt': 5.55, '$gt': 4.45}},
                {'encapsulated': {'$exists': True}},
                {'cnt_radius': {'$gt': 5.20, '$lt': 5.22}},
                {'sedc_scheme': {'$exists': False}},
                {'kpoints_mp_spacing': {'$gte': 0.045000000000000005, '$lte': 0.055}},
                {'spin_polarized': {'$ne': True}},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

    def testParseElementStr(self):
        arg = '[VII][Fe,Ru,Os][I]'
        elements = parse_element_string(arg)
        self.assertEqual(elements, ['[VII]', '[Fe,Ru,Os]', '[I]'])

        arg = '[VII][Fe,Ru,Os][I][V][VIII][ASDASD]'
        elements = parse_element_string(arg)
        self.assertEqual(elements, ['[VII]', '[Fe,Ru,Os]', '[I]', '[V]', '[VIII]', '[ASDASD]'])

        arg = '[VII]5[Fe,Ru,Os]2[I][V]6[VIII]2[ASDASD]'
        elements = parse_element_string(arg, stoich=True)
        self.assertEqual(elements, ['[VII]', '5', '[Fe,Ru,Os]', '2', '[I]', '[V]', '6', '[VIII]', '2', '[ASDASD]'])

    def testHarderCompositions(self):
        """ Test awkward compositions with macros etc. """
        kwargs = {'composition': ['[I]FeBe'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$or': [
                        {'atom_types': {'$in': ['Li']}},
                        {'atom_types': {'$in': ['Na']}},
                        {'atom_types': {'$in': ['K']}},
                        {'atom_types': {'$in': ['Rb']}},
                        {'atom_types': {'$in': ['Cs']}},
                        {'atom_types': {'$in': ['Fr']}},
                    ]},
                    {'atom_types': {'$in': ['Fe']}},
                    {'atom_types': {'$in': ['Be']}},
                    {'stoichiometry': {'$size': 3}}
                    ]},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': ['LiFeBe'], 'ignore_warnings': True, 'intersection': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$or': [
                    {'$and': [
                        {'stoichiometry': {'$size': 1}},
                        {'atom_types': {'$in': ['Li']}},
                    ]},
                    {'$and': [
                        {'stoichiometry': {'$size': 1}},
                        {'atom_types': {'$in': ['Fe']}},
                    ]},
                    {'$and': [
                        {'stoichiometry': {'$size': 1}},
                        {'atom_types': {'$in': ['Be']}},
                    ]},
                    {'$and': [
                        {'stoichiometry': {'$size': 2}},
                        {'atom_types': {'$in': ['Li']}},
                        {'atom_types': {'$in': ['Fe']}},
                    ]},
                    {'$and': [
                        {'stoichiometry': {'$size': 2}},
                        {'atom_types': {'$in': ['Li']}},
                        {'atom_types': {'$in': ['Be']}},
                    ]},
                    {'$and': [
                        {'stoichiometry': {'$size': 2}},
                        {'atom_types': {'$in': ['Fe']}},
                        {'atom_types': {'$in': ['Be']}},
                    ]},
                    {'$and': [
                        {'stoichiometry': {'$size': 3}},
                        {'atom_types': {'$in': ['Li']}},
                        {'atom_types': {'$in': ['Fe']}},
                        {'atom_types': {'$in': ['Be']}},
                    ]},
                ]}
            ]})
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': ['Fe[I]Be'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'atom_types': {'$in': ['Fe']}},
                    {'$or': [
                        {'atom_types': {'$in': ['Li']}},
                        {'atom_types': {'$in': ['Na']}},
                        {'atom_types': {'$in': ['K']}},
                        {'atom_types': {'$in': ['Rb']}},
                        {'atom_types': {'$in': ['Cs']}},
                        {'atom_types': {'$in': ['Fr']}},
                    ]},
                    {'atom_types': {'$in': ['Be']}},
                    {'stoichiometry': {'$size': 3}}
                    ]},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': ['[Fe,Ru,Os][I]Be'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$or': [
                        {'atom_types': {'$in': ['Fe']}},
                        {'atom_types': {'$in': ['Ru']}},
                        {'atom_types': {'$in': ['Os']}}
                    ]},
                    {'$or': [
                        {'atom_types': {'$in': ['Li']}},
                        {'atom_types': {'$in': ['Na']}},
                        {'atom_types': {'$in': ['K']}},
                        {'atom_types': {'$in': ['Rb']}},
                        {'atom_types': {'$in': ['Cs']}},
                        {'atom_types': {'$in': ['Fr']}}
                    ]},
                    {'atom_types': {'$in': ['Be']}},
                    {'stoichiometry': {'$size': 3}}
                    ]},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': ['Be[Fe,Ru,Os][I]'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'atom_types': {'$in': ['Be']}},
                    {'$or': [
                        {'atom_types': {'$in': ['Fe']}},
                        {'atom_types': {'$in': ['Ru']}},
                        {'atom_types': {'$in': ['Os']}}
                    ]},
                    {'$or': [
                        {'atom_types': {'$in': ['Li']}},
                        {'atom_types': {'$in': ['Na']}},
                        {'atom_types': {'$in': ['K']}},
                        {'atom_types': {'$in': ['Rb']}},
                        {'atom_types': {'$in': ['Cs']}},
                        {'atom_types': {'$in': ['Fr']}}
                    ]},
                    {'stoichiometry': {'$size': 3}}
                    ]},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': ['[VII][Fe,Ru,Os][I]'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$or': [
                        {'atom_types': {'$in': ['F']}},
                        {'atom_types': {'$in': ['Cl']}},
                        {'atom_types': {'$in': ['Br']}},
                        {'atom_types': {'$in': ['I']}},
                        {'atom_types': {'$in': ['At']}},
                    ]},
                    {'$or': [
                        {'atom_types': {'$in': ['Fe']}},
                        {'atom_types': {'$in': ['Ru']}},
                        {'atom_types': {'$in': ['Os']}}
                    ]},
                    {'$or': [
                        {'atom_types': {'$in': ['Li']}},
                        {'atom_types': {'$in': ['Na']}},
                        {'atom_types': {'$in': ['K']}},
                        {'atom_types': {'$in': ['Rb']}},
                        {'atom_types': {'$in': ['Cs']}},
                        {'atom_types': {'$in': ['Fr']}}
                    ]},
                    {'stoichiometry': {'$size': 3}}
                    ]},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': ['[Si,Ge,Sn][Fe,Ru,Os]'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$or': [
                        {'atom_types': {'$in': ['Si']}},
                        {'atom_types': {'$in': ['Ge']}},
                        {'atom_types': {'$in': ['Sn']}}
                    ]},
                    {'$or': [
                        {'atom_types': {'$in': ['Fe']}},
                        {'atom_types': {'$in': ['Ru']}},
                        {'atom_types': {'$in': ['Os']}}
                    ]},
                    {'stoichiometry': {'$size': 2}}
                    ]},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))

        self.assertDictEqual(test_dict, query.query_dict)

    def testTrickyStoichs(self):
        kwargs = {'formula': ['[VII]2[Fe,Ru,Os]3[I]'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$or': [
                        {'stoichiometry': {'$in': [['F', 2.0]]}},
                        {'stoichiometry': {'$in': [['Cl', 2.0]]}},
                        {'stoichiometry': {'$in': [['Br', 2.0]]}},
                        {'stoichiometry': {'$in': [['I', 2.0]]}},
                        {'stoichiometry': {'$in': [['At', 2.0]]}},
                    ]},
                    {'$or': [
                        {'stoichiometry': {'$in': [['Fe', 3.0]]}},
                        {'stoichiometry': {'$in': [['Ru', 3.0]]}},
                        {'stoichiometry': {'$in': [['Os', 3.0]]}}
                    ]},
                    {'$or': [
                        {'stoichiometry': {'$in': [['Li', 1.0]]}},
                        {'stoichiometry': {'$in': [['Na', 1.0]]}},
                        {'stoichiometry': {'$in': [['K', 1.0]]}},
                        {'stoichiometry': {'$in': [['Rb', 1.0]]}},
                        {'stoichiometry': {'$in': [['Cs', 1.0]]}},
                        {'stoichiometry': {'$in': [['Fr', 1.0]]}}
                    ]},
                    {'stoichiometry': {'$size': 3}}
                    ]},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'formula': ['[Ag,Cd,In]2[Fe,Ru,Os]3[I]'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$or': [
                        {'stoichiometry': {'$in': [['Ag', 2.0]]}},
                        {'stoichiometry': {'$in': [['Cd', 2.0]]}},
                        {'stoichiometry': {'$in': [['In', 2.0]]}}
                    ]},
                    {'$or': [
                        {'stoichiometry': {'$in': [['Fe', 3.0]]}},
                        {'stoichiometry': {'$in': [['Ru', 3.0]]}},
                        {'stoichiometry': {'$in': [['Os', 3.0]]}}
                    ]},
                    {'$or': [
                        {'stoichiometry': {'$in': [['Li', 1.0]]}},
                        {'stoichiometry': {'$in': [['Na', 1.0]]}},
                        {'stoichiometry': {'$in': [['K', 1.0]]}},
                        {'stoichiometry': {'$in': [['Rb', 1.0]]}},
                        {'stoichiometry': {'$in': [['Cs', 1.0]]}},
                        {'stoichiometry': {'$in': [['Fr', 1.0]]}}
                    ]},
                    {'stoichiometry': {'$size': 3}}
                    ]},
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))

    def testRatioQuery(self):
        kwargs = {'composition': ['Li:TiP4'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$and': [
                        {'atom_types': {'$in': ['Li']}},
                    ]},
                    {'ratios.TiP': 0.25},
                    {'$and': [
                        {'atom_types': {'$in': ['Ti']}},
                    ]},
                    {'$and': [
                        {'atom_types': {'$in': ['P']}},
                    ]},
                    {'stoichiometry': {'$size': 3}}
                ]}
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': ['LiMn:Mo2S3'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$and': [
                        {'atom_types': {'$in': ['Li']}},
                    ]},
                    {'$and': [
                        {'atom_types': {'$in': ['Mn']}},
                    ]},
                    {'ratios.MoS': 0.667},
                    {'$and': [
                        {'atom_types': {'$in': ['Mo']}},
                    ]},
                    {'$and': [
                        {'atom_types': {'$in': ['S']}},
                    ]},
                    {'stoichiometry': {'$size': 4}}
                ]}
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {'composition': ['LiMn:Mo2S3B5'], 'ignore_warnings': True,
                  'testing': True}
        query = DBQuery(**kwargs)
        test_dict = ({
            '$and': [
                {'$and': [
                    {'$and': [
                        {'atom_types': {'$in': ['Li']}},
                    ]},
                    {'$and': [
                        {'atom_types': {'$in': ['Mn']}},
                    ]},
                    {'ratios.MoS': 0.667, 'ratios.SB': 0.6, 'ratios.MoB': 0.4},
                    {'$and': [
                        {'atom_types': {'$in': ['Mo']}},
                    ]},
                    {'$and': [
                        {'atom_types': {'$in': ['S']}},
                    ]},
                    {'$and': [
                        {'atom_types': {'$in': ['B']}},
                    ]},
                    {'stoichiometry': {'$size': 5}}
                ]}
            ]
        })
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)


if __name__ == '__main__':
    unittest.main()
