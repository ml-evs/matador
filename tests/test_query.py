#!/usr/bin/env python
import unittest
import re

from matador.query import DBQuery
from matador.utils.chem_utils import parse_element_string


class QueryParseTest(unittest.TestCase):
    """ Test query functionality. """

    def test_basic_queries(self):
        kwargs = {"composition": "KP", "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["K"]}},
                        {"elems": {"$in": ["P"]}},
                        {"stoichiometry": {"$size": 2}},
                    ]
                },
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {"formula": "K3P4", "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"stoichiometry": {"$in": [["K", 3]]}},
                        {"stoichiometry": {"$in": [["P", 4]]}},
                        {"stoichiometry": {"$size": 2}},
                    ]
                },
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {"formula": "K3P6", "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"stoichiometry": {"$in": [["K", 1]]}},
                        {"stoichiometry": {"$in": [["P", 2]]}},
                        {"stoichiometry": {"$size": 2}},
                    ]
                },
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {"formula": "K27P9", "partial_formula": True, "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"stoichiometry": {"$in": [["K", 3]]}},
                        {"stoichiometry": {"$in": [["P", 1]]}},
                    ]
                },
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {"formula": ["K27P9"], "partial_formula": True, "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"stoichiometry": {"$in": [["K", 3]]}},
                        {"stoichiometry": {"$in": [["P", 1]]}},
                    ]
                },
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {
            "field": ["xc_functional", "cut_off_energy"],
            "filter": [["PBE"], ["301.2312", 400.0]],
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {"xc_functional": "PBE"},
                {"cut_off_energy": {"$gte": 301.2312, "$lte": 400.0}},
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {
            "field": ["cut_off_energy", "xc_functional"],
            "filter": [["301.2312", 400.0], ["PBE", "LDA"]],
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {"cut_off_energy": {"$gte": 301.2312, "$lte": 400.0}},
                {"$or": [{"xc_functional": "PBE"}, {"xc_functional": "LDA"}]},
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_complex_queries(self):
        """ Test long queries with multiple mismatching of lists
        (emulating argparse) and values.
        """
        self.maxDiff = None
        kwargs = {
            "formula": "K27P9",
            "space_group": ["Fd-3m"],
            "num_species": 3,
            "num_fu": [4],
            "tags": ["foo", "Bar"],
            "doi": ["1001/4001"],
            "icsd": 100020,
            "cutoff": 300,
            "spin": 1,
            "partial_formula": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"stoichiometry": {"$in": [["K", 3]]}},
                        {"stoichiometry": {"$in": [["P", 1]]}},
                    ]
                },
                {"stoichiometry": {"$size": 3}},
                {"space_group": "Fd-3m"},
                {"num_fu": {"$gte": 4}},
                {"$and": [{"tags": {"$in": ["foo"]}}, {"tags": {"$in": ["Bar"]}}]},
                {"doi": {"$in": ["1001/4001"]}},
                {"$or": [{'icsd': {"$eq": "100020"}}, {'icsd': {'$eq': 100020}}]},
                {"cut_off_energy": {"$eq": 300}},
                {"spin_polarized": True},
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {
            "composition": "LiFeBe",
            "icsd": True,
            "ignore_warnings": True,
            "src_str": "/Foo/bar/foo/Bar.res",
            "pressure": 5,
            "pressure_tolerance": 1,
            "cutoff": [300, 400],
            "encapsulated": True,
            "cnt_radius": 5.21,
            "sedc": "null",
            "mp_spacing": [0.05],
            "spin": 0,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["Li"]}},
                        {"elems": {"$in": ["Fe"]}},
                        {"elems": {"$in": ["Be"]}},
                        {"stoichiometry": {"$size": 3}},
                    ]
                },
                {"icsd": {"$exists": True}},
                {"cut_off_energy": {"$gte": 300, "$lte": 400}},
                {"source": {"$in": [re.compile("/Foo/bar/foo/Bar.res")]}},
                {"pressure": {"$lte": 6.0, "$gte": 4.0}},
                {"encapsulated": {"$exists": True}},
                {"cnt_radius": {"$gte": 5.20, "$lte": 5.22}},
                {"sedc_scheme": {"$exists": False}},
                {"kpoints_mp_spacing": {"$gte": 0.04, "$lte": 0.06}},
                {"spin_polarized": {"$ne": True}},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {
            "composition": "LiFeBe",
            "icsd": False,
            "ignore_warnings": True,
            "src_str": "/Foo/bar/foo/Bar.res",
            "pressure": 5,
            "pressure_tolerance": None,
            "kpoint_tolerance": None,
            "cutoff": [300, 400],
            "encapsulated": True,
            "cnt_radius": 5.21,
            "sedc": "null",
            "mp_spacing": [0.05],
            "spin": "any",
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["Li"]}},
                        {"elems": {"$in": ["Fe"]}},
                        {"elems": {"$in": ["Be"]}},
                        {"stoichiometry": {"$size": 3}},
                    ]
                },
                {"icsd": {"$exists": False}},
                {"cut_off_energy": {"$gte": 300, "$lte": 400}},
                {"source": {"$in": [re.compile("/Foo/bar/foo/Bar.res")]}},
                {"pressure": {"$lte": 5.5, "$gte": 4.5}},
                {"encapsulated": {"$exists": True}},
                {"cnt_radius": {"$gte": 5.20, "$lte": 5.22}},
                {"sedc_scheme": {"$exists": False}},
                {"kpoints_mp_spacing": {"$lte": 0.06, "$gte": 0.04}},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {
            "composition": "LiFeBe",
            "icsd": True,
            "ignore_warnings": True,
            "src_str": "/Foo/bar/foo/Bar.res",
            "pressure": 5,
            "pressure_tolerance": 10,
            "kpoint_tolerance": 0.02,
            "cutoff": [300, 400],
            "encapsulated": True,
            "cnt_radius": 5.21,
            "sedc": "null",
            "mp_spacing": [0.04, 0.05],
            "spin": "any",
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["Li"]}},
                        {"elems": {"$in": ["Fe"]}},
                        {"elems": {"$in": ["Be"]}},
                        {"stoichiometry": {"$size": 3}},
                    ]
                },
                {"icsd": {"$exists": True}},
                {"cut_off_energy": {"$gte": 300, "$lte": 400}},
                {"source": {"$in": [re.compile("/Foo/bar/foo/Bar.res")]}},
                {"pressure": {"$lte": 15, "$gte": -5}},
                {"encapsulated": {"$exists": True}},
                {"cnt_radius": {"$gte": 5.20, "$lte": 5.22}},
                {"sedc_scheme": {"$exists": False}},
                {"kpoints_mp_spacing": {"$lte": 0.05, "$gte": 0.04}},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_non_elemental_query(self):
        kwargs = {"composition": "Li4Ge:Fe2Be", "testing": True}
        query = DBQuery(**kwargs)

        kwargs = {"composition": "LiGeFeBe", "intersection": True, "testing": True}
        query2 = DBQuery(**kwargs)
        print(query.query_dict)
        print(query2.query_dict)

        self.assertDictEqual(query.query_dict, query2.query_dict)
        self.assertEqual(query._chempots, ["Li4Ge", "Fe2Be"])

    def test_duplicate_elements(self):
        kwargs = {"composition": "PPPPP", "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {"$and": [{"elems": {"$in": ["P"]}}, {"stoichiometry": {"$size": 1}}]},
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }

        self.assertDictEqual(query.query_dict, test_dict)

    def test_parse_element_strings(self):
        arg = "[VII][Fe,Ru,Os][I]"
        elements = parse_element_string(arg)
        self.assertEqual(elements, ["[VII]", "[Fe,Ru,Os]", "[I]"])

        arg = "[VII][Fe,Ru,Os][I][V][VIII][ASDASD]"
        elements = parse_element_string(arg)
        self.assertEqual(
            elements, ["[VII]", "[Fe,Ru,Os]", "[I]", "[V]", "[VIII]", "[ASDASD]"]
        )

        arg = "{VII}[Fe,Ru,Os][I]{V}[VIII][ASDASD]"
        elements = parse_element_string(arg)
        self.assertEqual(
            elements, ["{VII}", "[Fe,Ru,Os]", "[I]", "{V}", "[VIII]", "[ASDASD]"]
        )

        arg = "{VII}[Fe,Ru,Os][I][V][VIII][ASDASD]"
        elements = parse_element_string(arg)
        self.assertEqual(
            elements, ["{VII}", "[Fe,Ru,Os]", "[I]", "[V]", "[VIII]", "[ASDASD]"]
        )

        arg = "[VII]5[Fe,Ru,Os]2[I][V]6[VIII]2[ASDASD]"
        elements = parse_element_string(arg, stoich=True)
        self.assertEqual(
            elements,
            [
                "[VII]",
                "5",
                "[Fe,Ru,Os]",
                "2",
                "[I]",
                "[V]",
                "6",
                "[VIII]",
                "2",
                "[ASDASD]",
            ],
        )

        raised = False
        arg = "{VII}[Fe,Ru,Os}[I][V][VIII][ASDASD]"
        try:
            elements = parse_element_string(arg)
        except RuntimeError:
            raised = True
        self.assertTrue(raised, msg="Failed to raise error for unmatched brackets")

    def test_macro_composition(self):
        kwargs = {"composition": ["[I]FeBe"], "ignore_warnings": True, "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Na"]}},
                                {"elems": {"$in": ["K"]}},
                                {"elems": {"$in": ["Rb"]}},
                                {"elems": {"$in": ["Cs"]}},
                                {"elems": {"$in": ["Fr"]}},
                            ]
                        },
                        {"elems": {"$in": ["Fe"]}},
                        {"elems": {"$in": ["Be"]}},
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_intersection(self):
        kwargs = {
            "composition": ["LiFeBe"],
            "ignore_warnings": True,
            "intersection": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$or": [
                        {
                            "$and": [
                                {"stoichiometry": {"$size": 1}},
                                {"elems": {"$in": ["Li"]}},
                            ]
                        },
                        {
                            "$and": [
                                {"stoichiometry": {"$size": 1}},
                                {"elems": {"$in": ["Fe"]}},
                            ]
                        },
                        {
                            "$and": [
                                {"stoichiometry": {"$size": 1}},
                                {"elems": {"$in": ["Be"]}},
                            ]
                        },
                        {
                            "$and": [
                                {"stoichiometry": {"$size": 2}},
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Fe"]}},
                            ]
                        },
                        {
                            "$and": [
                                {"stoichiometry": {"$size": 2}},
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Be"]}},
                            ]
                        },
                        {
                            "$and": [
                                {"stoichiometry": {"$size": 2}},
                                {"elems": {"$in": ["Fe"]}},
                                {"elems": {"$in": ["Be"]}},
                            ]
                        },
                        {
                            "$and": [
                                {"stoichiometry": {"$size": 3}},
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Fe"]}},
                                {"elems": {"$in": ["Be"]}},
                            ]
                        },
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_middle_macro_composition(self):
        kwargs = {"composition": ["Fe[I]Be"], "ignore_warnings": True, "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["Fe"]}},
                        {
                            "$or": [
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Na"]}},
                                {"elems": {"$in": ["K"]}},
                                {"elems": {"$in": ["Rb"]}},
                                {"elems": {"$in": ["Cs"]}},
                                {"elems": {"$in": ["Fr"]}},
                            ]
                        },
                        {"elems": {"$in": ["Be"]}},
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_double_macro_composition(self):
        kwargs = {
            "composition": ["[Fe,Ru,Os][I]Be"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {"elems": {"$in": ["Fe"]}},
                                {"elems": {"$in": ["Ru"]}},
                                {"elems": {"$in": ["Os"]}},
                            ]
                        },
                        {
                            "$or": [
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Na"]}},
                                {"elems": {"$in": ["K"]}},
                                {"elems": {"$in": ["Rb"]}},
                                {"elems": {"$in": ["Cs"]}},
                                {"elems": {"$in": ["Fr"]}},
                            ]
                        },
                        {"elems": {"$in": ["Be"]}},
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_double_end_macro_composition(self):
        kwargs = {
            "composition": ["Be[Fe,Ru,Os][I]"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["Be"]}},
                        {
                            "$or": [
                                {"elems": {"$in": ["Fe"]}},
                                {"elems": {"$in": ["Ru"]}},
                                {"elems": {"$in": ["Os"]}},
                            ]
                        },
                        {
                            "$or": [
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Na"]}},
                                {"elems": {"$in": ["K"]}},
                                {"elems": {"$in": ["Rb"]}},
                                {"elems": {"$in": ["Cs"]}},
                                {"elems": {"$in": ["Fr"]}},
                            ]
                        },
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_triple_macro_composition(self):
        kwargs = {
            "composition": ["[VII][Fe,Ru,Os][I]"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {"elems": {"$in": ["F"]}},
                                {"elems": {"$in": ["Cl"]}},
                                {"elems": {"$in": ["Br"]}},
                                {"elems": {"$in": ["I"]}},
                                {"elems": {"$in": ["At"]}},
                            ]
                        },
                        {
                            "$or": [
                                {"elems": {"$in": ["Fe"]}},
                                {"elems": {"$in": ["Ru"]}},
                                {"elems": {"$in": ["Os"]}},
                            ]
                        },
                        {
                            "$or": [
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Na"]}},
                                {"elems": {"$in": ["K"]}},
                                {"elems": {"$in": ["Rb"]}},
                                {"elems": {"$in": ["Cs"]}},
                                {"elems": {"$in": ["Fr"]}},
                            ]
                        },
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        # print(json.dumps(test_dict, indent=2))
        # print(json.dumps(query.query_dict, indent=2))
        self.assertDictEqual(test_dict, query.query_dict)

    def test_double_list_composition(self):
        kwargs = {
            "composition": ["[Si,Ge,Sn][Fe,Ru,Os]"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {"elems": {"$in": ["Si"]}},
                                {"elems": {"$in": ["Ge"]}},
                                {"elems": {"$in": ["Sn"]}},
                            ]
                        },
                        {
                            "$or": [
                                {"elems": {"$in": ["Fe"]}},
                                {"elems": {"$in": ["Ru"]}},
                                {"elems": {"$in": ["Os"]}},
                            ]
                        },
                        {"stoichiometry": {"$size": 2}},
                    ]
                }
            ]
        }

        self.assertDictEqual(test_dict, query.query_dict)

    def test_triple_list_composition(self):
        kwargs = {
            "composition": ["[Li,Na,K][Fe,Ru,Os][F,P,S]"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {"elems": {"$in": ["Li"]}},
                                {"elems": {"$in": ["Na"]}},
                                {"elems": {"$in": ["K"]}},
                            ]
                        },
                        {
                            "$or": [
                                {"elems": {"$in": ["Fe"]}},
                                {"elems": {"$in": ["Ru"]}},
                                {"elems": {"$in": ["Os"]}},
                            ]
                        },
                        {
                            "$or": [
                                {"elems": {"$in": ["F"]}},
                                {"elems": {"$in": ["P"]}},
                                {"elems": {"$in": ["S"]}},
                            ]
                        },
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        print(query.query_dict)
        self.assertDictEqual(test_dict, query.query_dict)

    def test_single_set_composition(self):
        kwargs = {"composition": ["{Li,Na}"], "ignore_warnings": True, "testing": True}
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {"$and": [{"elems": {"$in": ["Li"]}}]},
                                {"$and": [{"elems": {"$in": ["Na"]}}]},
                            ]
                        },
                        {"stoichiometry": {"$size": 1}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_double_set_composition(self):
        kwargs = {
            "composition": ["{Li,Na}{Be,Fe}"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Li"]}},
                                        {"elems": {"$in": ["Be"]}},
                                    ]
                                },
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Li"]}},
                                        {"elems": {"$in": ["Fe"]}},
                                    ]
                                },
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Na"]}},
                                        {"elems": {"$in": ["Be"]}},
                                    ]
                                },
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Na"]}},
                                        {"elems": {"$in": ["Fe"]}},
                                    ]
                                },
                            ]
                        },
                        {"stoichiometry": {"$size": 2}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_ternary_set_composition(self):
        kwargs = {
            "composition": ["LiRu{Ru,S}"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Li"]}},
                                        {"elems": {"$in": ["Ru"]}},
                                        {"elems": {"$in": ["S"]}},
                                    ]
                                }
                            ]
                        },
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_solo_set_composition(self):
        kwargs = {
            "composition": ["{Li}Ru{Ru,S}"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Li"]}},
                                        {"elems": {"$in": ["Ru"]}},
                                        {"elems": {"$in": ["S"]}},
                                    ]
                                }
                            ]
                        },
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_triple_set_composition(self):
        kwargs = {
            "composition": ["{Li,Na}{Ru,Os}{Ru,S}"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Li"]}},
                                        {"elems": {"$in": ["Ru"]}},
                                        {"elems": {"$in": ["S"]}},
                                    ]
                                },
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Li"]}},
                                        {"elems": {"$in": ["Os"]}},
                                        {"elems": {"$in": ["Ru"]}},
                                    ]
                                },
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Li"]}},
                                        {"elems": {"$in": ["Os"]}},
                                        {"elems": {"$in": ["S"]}},
                                    ]
                                },
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Na"]}},
                                        {"elems": {"$in": ["Ru"]}},
                                        {"elems": {"$in": ["S"]}},
                                    ]
                                },
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Na"]}},
                                        {"elems": {"$in": ["Os"]}},
                                        {"elems": {"$in": ["Ru"]}},
                                    ]
                                },
                                {
                                    "$and": [
                                        {"elems": {"$in": ["Na"]}},
                                        {"elems": {"$in": ["Os"]}},
                                        {"elems": {"$in": ["S"]}},
                                    ]
                                },
                            ]
                        },
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_illegal_character(self):
        kwargs = {
            "composition": ["{Li,Na,K}.{Fe,Ru,Os>}{F,P,S}"],
            "ignore_warnings": True,
            "testing": True,
        }
        with self.assertRaises(RuntimeError):
            DBQuery(**kwargs)

    def test_triple_macro_stoich(self):
        kwargs = {
            "formula": ["[VII]2[Fe,Ru,Os]3[I]"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {"stoichiometry": {"$in": [["F", 2]]}},
                                {"stoichiometry": {"$in": [["Cl", 2]]}},
                                {"stoichiometry": {"$in": [["Br", 2]]}},
                                {"stoichiometry": {"$in": [["I", 2]]}},
                                {"stoichiometry": {"$in": [["At", 2]]}},
                            ]
                        },
                        {
                            "$or": [
                                {"stoichiometry": {"$in": [["Fe", 3]]}},
                                {"stoichiometry": {"$in": [["Ru", 3]]}},
                                {"stoichiometry": {"$in": [["Os", 3]]}},
                            ]
                        },
                        {
                            "$or": [
                                {"stoichiometry": {"$in": [["Li", 1]]}},
                                {"stoichiometry": {"$in": [["Na", 1]]}},
                                {"stoichiometry": {"$in": [["K", 1]]}},
                                {"stoichiometry": {"$in": [["Rb", 1]]}},
                                {"stoichiometry": {"$in": [["Cs", 1]]}},
                                {"stoichiometry": {"$in": [["Fr", 1]]}},
                            ]
                        },
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        print(test_dict)
        print(query.query_dict)
        self.assertDictEqual(test_dict, query.query_dict)

        kwargs = {
            "formula": ["[Ag,Cd,In]2[Fe,Ru,Os]3[I]"],
            "ignore_warnings": True,
            "testing": True,
        }
        query = DBQuery(**kwargs)
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {
                            "$or": [
                                {"stoichiometry": {"$in": [["Ag", 2]]}},
                                {"stoichiometry": {"$in": [["Cd", 2]]}},
                                {"stoichiometry": {"$in": [["In", 2]]}},
                            ]
                        },
                        {
                            "$or": [
                                {"stoichiometry": {"$in": [["Fe", 3]]}},
                                {"stoichiometry": {"$in": [["Ru", 3]]}},
                                {"stoichiometry": {"$in": [["Os", 3]]}},
                            ]
                        },
                        {
                            "$or": [
                                {"stoichiometry": {"$in": [["Li", 1]]}},
                                {"stoichiometry": {"$in": [["Na", 1]]}},
                                {"stoichiometry": {"$in": [["K", 1]]}},
                                {"stoichiometry": {"$in": [["Rb", 1]]}},
                                {"stoichiometry": {"$in": [["Cs", 1]]}},
                                {"stoichiometry": {"$in": [["Fr", 1]]}},
                            ]
                        },
                        {"stoichiometry": {"$size": 3}},
                    ]
                }
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_query_time(self):
        from bson.objectid import ObjectId
        from datetime import datetime, timedelta
        from time import mktime

        num_days = 5
        kwargs = {"composition": "KP", "testing": True, "time": num_days}
        query = DBQuery(**kwargs)
        days_ago = (datetime.today() - timedelta(days=num_days)).timetuple()
        time_str = str(hex(int(mktime(days_ago))))[2:]
        time_str += (24 - len(time_str)) * "0"
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["K"]}},
                        {"elems": {"$in": ["P"]}},
                        {"stoichiometry": {"$size": 2}},
                    ]
                },
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
                {"_id": {"$lte": ObjectId(time_str)}},
            ]
        }
        print(test_dict["$and"][0])
        print(query.query_dict["$and"][0])
        self.assertDictEqual(test_dict, query.query_dict)

        num_days = 3001
        kwargs = {"composition": "KP", "testing": True, "time": num_days}
        query = DBQuery(**kwargs)
        days_ago = (datetime.today() - timedelta(days=num_days)).timetuple()
        time_str = str(hex(int(mktime(days_ago))))[2:]
        time_str += (24 - len(time_str)) * "0"
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["K"]}},
                        {"elems": {"$in": ["P"]}},
                        {"stoichiometry": {"$size": 2}},
                    ]
                },
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
                {"_id": {"$lte": ObjectId(time_str)}},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

        num_days = -1
        kwargs = {"composition": "KP", "testing": True, "time": num_days}
        query = DBQuery(**kwargs)
        days_ago = (datetime.today() - timedelta(days=num_days)).timetuple()
        time_str = str(hex(int(mktime(days_ago))))[2:]
        time_str += (24 - len(time_str)) * "0"
        test_dict = {
            "$and": [
                {
                    "$and": [
                        {"elems": {"$in": ["K"]}},
                        {"elems": {"$in": ["P"]}},
                        {"stoichiometry": {"$size": 2}},
                    ]
                },
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
                {"_id": {"$lte": ObjectId(time_str)}},
            ]
        }
        self.assertDictEqual(test_dict, query.query_dict)

    def test_id_query(self):
        kwargs = {"id": "testing testing", "ignore_warnings": False, "testing": True}
        query = DBQuery(**kwargs)

        test_dict = {
            "$and": [
                {"text_id": ["testing", "testing"]},
                {"$or": [{"quality": {"$gt": 0}}, {"quality": {"$exists": False}}]},
            ]
        }

        print(test_dict)
        print(query.query_dict)
        self.assertDictEqual(test_dict, query.query_dict)


if __name__ == "__main__":
    unittest.main()
