#!/usr/bin/env python

import unittest
import os
import numpy as np
from matador.utils.chem_utils import get_concentration
from matador.utils.chem_utils import (
    get_generic_grav_capacity,
    get_binary_volumetric_capacity,
)
from matador.utils.chem_utils import (
    get_stoich,
    get_formula_from_stoich,
    get_stoich_from_formula,
)
from matador.utils.chem_utils import get_ratios_from_stoichiometry, get_root_source
from matador.utils.chem_utils import get_formation_energy
from matador.utils.chem_utils import get_padded_composition

REAL_PATH = "/".join(os.path.realpath(__file__).split("/")[:-1]) + "/"


class ChemUtilsTest(unittest.TestCase):
    """Test chem utils functionality."""

    def test_gravimetric_capacity(self):
        test_docs = []
        test_elements = []
        Q = []

        # Li3P and symmetry test
        doc = dict()
        doc["stoichiometry"] = [["Li", 3], ["P", 1]]
        test_docs.append(doc)
        test_elements.append(["Li", "P"])
        doc = dict()
        doc["stoichiometry"] = [["P", 1], ["Li", 3]]
        test_docs.append(doc)
        test_elements.append(["Li", "P"])
        # ternary test
        doc = dict()
        doc["stoichiometry"] = [["Li", 1], ["Mo", 1], ["S", 2]]
        test_docs.append(doc)
        test_elements.append(["Li", "Mo", "S"])
        doc = dict()
        doc["stoichiometry"] = [["Li", 8], ["Mo", 1], ["S", 2]]
        test_docs.append(doc)
        test_elements.append(["Li", "Mo", "S"])
        doc = dict()
        doc["stoichiometry"] = [["Li", 1], ["Mo", 2], ["S", 4]]
        test_docs.append(doc)
        test_elements.append(["Li", "Mo", "S"])
        for doc, elem in zip(test_docs, test_elements):
            doc["concentration"] = get_concentration(doc, elem)
            temp_conc = list(doc["concentration"])
            temp_conc.append(1.0)
            for conc in doc["concentration"]:
                temp_conc[-1] -= conc

            Q.append(get_generic_grav_capacity(temp_conc, elem))
        self.assertAlmostEqual(Q[0], 2596.096626125, places=3)
        self.assertAlmostEqual(Q[2], 167.449, places=3)
        self.assertEqual(Q[0], Q[1])
        self.assertEqual(round(8 * Q[2], 3), round(Q[3], 3))
        self.assertEqual(round(Q[2], 3), round(2 * Q[4], 3))

    def test_get_concentration(self):
        doc = [["Na", 3], ["P", 1]]
        conc = get_concentration(doc, ["Na", "P"])
        self.assertEqual(conc, [0.75])
        conc = get_concentration(doc, ["Na", "P"], include_end=True)
        self.assertEqual(conc, [0.75, 0.25])
        conc = get_concentration(doc, ["K", "P"], include_end=True)
        self.assertEqual(conc, [0.0, 0.25])

        doc = {"stoichiometry": [["Na", 3], ["P", 1]]}
        conc = get_concentration(doc, ["Na", "P"])
        self.assertEqual(conc, [0.75])
        conc = get_concentration(doc, ["Na", "P"], include_end=True)
        self.assertEqual(conc, [0.75, 0.25])

    def test_volumetric_capacity(self):
        initial_doc = dict()
        final_doc = dict()
        initial_doc["stoichiometry"] = [["P", 1]]
        initial_doc["cell_volume"] = 84.965349
        initial_doc["num_fu"] = 4

        final_doc["stoichiometry"] = sorted([["Li", 3], ["P", 1]])
        vol_cap = get_binary_volumetric_capacity(initial_doc, final_doc)
        self.assertAlmostEqual(vol_cap, 6286, places=0)

    def test_atoms_to_stoich(self):
        atoms = 5 * ["Li"]
        atoms.extend(5 * ["P"])
        stoich = [["Li", 1], ["P", 1]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 99 * ["Li"]
        atoms.extend(1 * ["P"])
        stoich = [["Li", 99], ["P", 1]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 4 * ["Li"]
        atoms.extend(36 * ["P"])
        stoich = [["Li", 1], ["P", 9]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 3 * ["Li"]
        atoms.extend(2 * ["P"])
        stoich = [["Li", 3], ["P", 2]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 9 * ["Li"]
        atoms.extend(6 * ["P"])
        stoich = [["Li", 3], ["P", 2]]
        self.assertEqual(stoich, get_stoich(atoms))
        atoms = 36 * ["P"]
        atoms.extend(4 * ["Li"])
        stoich = [["Li", 1], ["P", 9]]
        self.assertEqual(stoich, get_stoich(atoms))

    def test_stoich_to_form(self):
        stoich = [["Li", 1], ["P", 9]]
        form = "LiP9"
        self.assertEqual(form, get_formula_from_stoich(stoich))
        stoich = [["P", 9], ["Li", 1]]
        form = "LiP9"
        self.assertEqual(form, get_formula_from_stoich(stoich))
        stoich = [["Li", 1], ["P", 9]]
        form = "LiP$_\\mathrm{9}$"
        self.assertEqual(
            form, get_formula_from_stoich(stoich, tex=True, latex_sub_style="\\mathrm")
        )
        stoich = [["Li", 1], ["P", 9]]
        form = "P9Li"
        self.assertEqual(form, get_formula_from_stoich(stoich, elements=["P", "Li"]))

    def test_form_to_stoich(self):
        formula = "Li12P1N18"
        stoich = [["Li", 12], ["P", 1], ["N", 18]]
        self.assertEqual(stoich, get_stoich_from_formula(formula, sort=False))

    def test_padded_composition(self):
        stoich = [["Li", 2], ["O", 1]]
        elements = ["O", "Li", "Ba"]
        self.assertEqual(get_padded_composition(stoich, elements), [1, 2, 0])
        elements = ["O", "Ba", "Li"]
        self.assertEqual(get_padded_composition(stoich, elements), [1, 0, 2])
        elements = ["Ba", "Li", "O"]
        self.assertEqual(get_padded_composition(stoich, elements), [0, 2, 1])
        stoich = [["F", 8]]
        elements = ["F"]
        self.assertEqual(get_padded_composition(stoich, elements), [8])
        elements = ["Fe", "Ba", "Ca"]
        with self.assertRaises(RuntimeError):
            get_padded_composition(stoich, elements)

    def test_formation_energy_binary(self):
        chempots = []
        chempots.append(
            {
                "stoichiometry": [["Li", 1]],
                "enthalpy": -200,
                "enthalpy_per_atom": -100,
                "temperature": {"free_energy_per_atom": {"300": -100, "400": -120}},
                "num_fu": 2,
                "atom_types": ["Li", "Li"],
            }
        )
        chempots.append(
            {
                "stoichiometry": [["P", 1]],
                "enthalpy": -600,
                "enthalpy_per_atom": -200,
                "num_fu": 3,
                "num_atoms": 3,
                "temperature": {
                    "free_energy_per_atom": {"300": -200, "400": -210, "450": -300}
                },
                "atom_types": ["P", "P", "P"],
            }
        )
        cursor = []
        cursor.append(
            {
                "stoichiometry": [["Li", 1], ["P", 1]],
                "enthalpy": -1200,
                "enthalpy_per_atom": -300,
                "num_fu": 2,
                "num_atoms": 4,
                "temperature": {"free_energy_per_atom": {"300": -300, "400": -315}},
                "atom_types": ["Li", "Li", "P", "P"],
            }
        )

        cursor.append(
            {
                "stoichiometry": [["Li", 3], ["P", 1]],
                "enthalpy": -1600,
                "enthalpy_per_atom": -400,
                "num_fu": 1,
                "num_atoms": 4,
                "atom_types": ["Li", "Li", "Li", "P"],
            }
        )

        ef = []
        for doc in cursor:
            ef.append(get_formation_energy(chempots, doc))

        self.assertEqual(ef[0], -600 / 4)
        self.assertEqual(ef[1], -1100 / 4)

        ef_300 = get_formation_energy(
            chempots,
            cursor[0],
            energy_key=["temperature", "free_energy_per_atom", "300"],
        )
        self.assertEqual(ef_300, -150)
        ef_400 = get_formation_energy(
            chempots,
            cursor[0],
            energy_key=["temperature", "free_energy_per_atom", "400"],
        )
        self.assertEqual(ef_400, -150)

        with self.assertRaises(KeyError):
            get_formation_energy(
                chempots,
                cursor[0],
                energy_key=["temperature", "free_energy_per_atom", "450"],
            )

    def test_formation_energy_nonbinary(self):
        from matador.utils.chem_utils import get_formation_energy

        chempots = []
        chempots.append(
            {
                "stoichiometry": [["Li", 2], ["O", 1]],
                "enthalpy": -1200,
                "enthalpy_per_atom": -200,
                "num_fu": 2,
                "num_atoms": 6,
                "atom_types": ["Li", "Li", "Li", "Li", "O", "O"],
            }
        )
        chempots.append(
            {
                "stoichiometry": [["La", 2], ["O", 3]],
                "enthalpy": -500,
                "enthalpy_per_atom": -100,
                "num_fu": 1,
                "num_atoms": 5,
                "atom_types": ["La", "La", "O", "O", "O"],
            }
        )
        chempots.append(
            {
                "stoichiometry": [["Zr", 1], ["O", 2]],
                "enthalpy": -4500,
                "enthalpy_per_atom": -150,
                "num_fu": 10,
                "num_atoms": 30,
            }
        )
        cursor = []
        cursor.append(
            {
                "stoichiometry": [["La", 2], ["Li", 4], ["O", 5]],
                "enthalpy": -1700,
                "enthalpy_per_atom": -1700 / 11,
                "num_fu": 1,
                "num_atoms": 11,
                "atom_types": [
                    "La",
                    "La",
                    "Li",
                    "Li",
                    "Li",
                    "Li",
                    "O",
                    "O",
                    "O",
                    "O",
                    "O",
                ],
            }
        )
        cursor.append(
            {
                "stoichiometry": [["La", 2], ["Li", 4], ["O", 5]],
                "enthalpy": -1711,
                "enthalpy_per_atom": -1711 / 11,
                "num_fu": 1,
                "num_atoms": 11,
                "atom_types": [
                    "La",
                    "La",
                    "Li",
                    "Li",
                    "Li",
                    "Li",
                    "O",
                    "O",
                    "O",
                    "O",
                    "O",
                ],
            }
        )
        cursor.append(
            {
                "stoichiometry": [["Li", 7], ["La", 3], ["Zr", 2], ["O", 12]],
                "enthalpy": -3.5 * 3 * 200 - 1.5 * 5 * 100 - 2 * 3 * -150,
                "enthalpy_per_atom": (-3.5 * 3 * 200 - 1.5 * 5 * 100 - 2 * 3 * 150)
                / 24,
            }
        )
        cursor.append(
            {
                "stoichiometry": [["Li", 7], ["La", 3], ["Zr", 2], ["O", 12]],
                "enthalpy": -3.5 * 3 * 200 - 1.5 * 5 * 100 - 2 * 3 * 150 - 48,
                "enthalpy_per_atom": (-3.5 * 3 * 200 - 1.5 * 5 * 100 - 2 * 3 * 150 - 48)
                / 24,
            }
        )

        ef = []
        for doc in cursor:
            ef.append(get_formation_energy(chempots, doc))

        np.testing.assert_array_almost_equal(ef, [0, -1, 0, -2])

    def test_get_number_of_chempots(self):
        from matador.utils.chem_utils import get_number_of_chempots

        stoich = [["V", 2]]
        chempots = [[["Na", 1]], [["S", 1], ["O", 4]], [["V", 1]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [0, 0, 2])

        stoich = [["Li", 7], ["La", 3], ["Zr", 2], ["O", 12]]
        chempots = [[["Li", 2], ["O", 1]], [["Zr", 1], ["O", 2]], [["La", 2], ["O", 3]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [3.5, 2, 1.5])

        stoich = [["Li", 7], ["La", 3], ["Zr", 2], ["O", 12]]
        chempots = [[["Li", 2], ["O", 1]], [["Zr", 1], ["O", 2]], [["La", 2], ["O", 3]]]
        np.testing.assert_array_almost_equal(
            get_number_of_chempots(stoich, chempots, precision=None), [3.5, 2, 1.5]
        )

        stoich = [["Li", 8], ["Zr", 3], ["O", 10]]
        chempots = [[["Li", 2], ["O", 1]], [["Zr", 1], ["O", 2]], [["La", 2], ["O", 3]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [4, 3, 0])

        stoich = [["Li", 3], ["P", 1]]
        chempots = [[["Li", 1]], [["Zn", 1]], [["P", 1]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [3, 0, 1])

        stoich = [["Li", 1]]
        chempots = [[["Li", 1]], [["Zn", 1]], [["P", 1]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [1, 0, 0])

        stoich = [["La", 2], ["Li", 4], ["O", 5]]
        chempots = [[["La", 2], ["O", 3]], [["Li", 2], ["O", 1]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [1, 2])

        stoich = [["Li", 1], ["P", 1], ["O", 3]]
        chempots = [[["Li", 2], ["S", 1]], [["P", 2], ["O", 5]], [["P", 2], ["S", 5]]]
        with self.assertRaises(RuntimeError):
            get_number_of_chempots(stoich, chempots)

        stoich = [["Na", 1], ["B", 1], ["S", 1], ["O", 4]]
        chempots = [[["Na", 1]], [["B", 1]], [["S", 1], ["O", 4]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [1, 1, 1])

        stoich = [["Na", 1]]
        chempots = [[["Na", 1]], [["B", 1]], [["S", 1], ["O", 4]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [1, 0, 0])

        stoich = [["B", 1]]
        chempots = [[["Na", 1]], [["B", 1]], [["S", 1], ["O", 4]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [0, 1, 0])

        stoich = [["S", 2], ["O", 8]]
        chempots = [[["Na", 1]], [["B", 1]], [["S", 1], ["O", 4]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [0, 0, 2])

        stoich = [["B", 2], ["Na", 9]]
        chempots = [[["Na", 1]], [["S", 1], ["O", 4]], [["B", 1]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [9, 0, 2])

        stoich = [["B", 2], ["Na", 9]]
        chempots = [[["Na", 1]], [["S", 1], ["O", 4]], [["B", 2]]]
        self.assertEqual(get_number_of_chempots(stoich, chempots), [9, 0, 1])

        stoich = [["A", 1], ["B", 2]]
        chempots = [[["A", 1], ["B", 2], ["C", 3]]]
        with self.assertRaises(RuntimeError):
            get_number_of_chempots(stoich, chempots)

    def test_ratios_from_stoich(self):
        stoich = [["Li", 12], ["N", 18], ["P", 1]]
        ratios = {
            "LiN": round(12.0 / 18, 3),
            "LiP": 12,
            "NP": 18,
            "NLi": round(18.0 / 12, 3),
            "PLi": round(1.0 / 12, 3),
            "PN": round(1.0 / 18, 3),
        }
        self.assertEqual(ratios, get_ratios_from_stoichiometry(stoich))

        stoich = [["K", 8], ["Sn", 1], ["P", 4]]
        ratios = {
            "KSn": 8,
            "KP": 2,
            "SnP": 0.25,
            "SnK": round(1.0 / 8, 3),
            "PK": 0.5,
            "PSn": 4,
        }
        self.assertEqual(ratios, get_ratios_from_stoichiometry(stoich))

    def test_root_source(self):
        source = ["KP.cell", "KP.param", "KP.castep"]
        src = "KP"
        self.assertEqual(src, get_root_source(source))

        source = ["KP.cell", "KP.param", "KP-1234-abcd.castep"]
        src = "KP-1234-abcd"
        self.assertEqual(src, get_root_source(source))

        source = [
            "KP.cell",
            "KP.param",
            "abcd-123.fdasf/efgf/KP-0.02.-1234-abcd.castep",
            "KP-0.02.-1234-abcd.res",
        ]
        src = "KP-0.02.-1234-abcd"
        self.assertEqual(src, get_root_source(source))

        source = [
            "KP.cell",
            "KP.param",
            "abcd-123.fdasf/efgf/KP-0.02.-1234-abcd.history",
        ]
        src = "KP-0.02.-1234-abcd"
        self.assertEqual(src, get_root_source(source))

        source = ["KP.cell", "KP.param", "PK-OQMD_12345.history"]
        src = "PK-OQMD_12345"
        self.assertEqual(src, get_root_source(source))

        source = ["OQMD 12345"]
        src = "OQMD 12345"
        self.assertEqual(src, get_root_source(source))

        source = ["OQMD-12345"]
        src = "OQMD-12345"
        self.assertEqual(src, get_root_source(source))

        source = [
            "KP.cell",
            "KP.param",
            "abcd-123.fdasf/efgf/KP-0.02.-1234-abcd.castep",
            "KP-1234-abcde.res",
        ]
        with self.assertRaises(RuntimeError):
            src = get_root_source(source)

        source = ["not a file name"]
        src = "not a file name"
        self.assertEqual(src, get_root_source(source))

        source = "not even a list"
        src = "not even a list"
        self.assertEqual(src, get_root_source(source))

    def test_magres_reference_shifts(self):
        from matador.utils.chem_utils import magres_reference_shifts

        magres = {
            "chemical_shielding_isos": [100, 75, -50, -25, 0, 100, 75, -50, -25, 0],
            "atom_types": ["H", "H", "H", "H", "H", "Li", "Li", "Li", "Li", "Li"],
        }

        magres_with_shifts = {
            "chemical_shielding_isos": [100, 75, -50, -25, 0, 100, 75, -50, -25, 0],
            "atom_types": ["H", "H", "H", "H", "H", "Li", "Li", "Li", "Li", "Li"],
            "chemical_shift_isos": [
                -49.0,
                -36.5,
                26.0,
                13.5,
                1.0,
                None,
                None,
                None,
                None,
                None,
            ],
        }

        magres_reference_shifts(magres=magres, reference={"H": [-0.5, 1]})

        self.assertEqual(magres, magres_with_shifts)


if __name__ == "__main__":
    unittest.main()
