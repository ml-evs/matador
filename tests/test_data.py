import unittest


class ScatteringDataLoaderTests(unittest.TestCase):
    """Test whether the atomic scattering factors can be loaded."""

    def test_load_gsas(self):
        from matador.data.atomic_scattering import GSAS_ATOMIC_SCATTERING_COEFFS

        self.assertEqual(len(GSAS_ATOMIC_SCATTERING_COEFFS["Cu"]), 3)
        for key in GSAS_ATOMIC_SCATTERING_COEFFS:
            self.assertEqual(len(GSAS_ATOMIC_SCATTERING_COEFFS[key]), 3)

    def test_load_raspa(self):
        from matador.data.atomic_scattering import RASPA_ATOMIC_SCATTERING_COEFFS

        self.assertEqual(len(RASPA_ATOMIC_SCATTERING_COEFFS["Cu"]), 3)
        for key in RASPA_ATOMIC_SCATTERING_COEFFS:
            self.assertEqual(len(RASPA_ATOMIC_SCATTERING_COEFFS[key]), 3)


class MagresDataLoaderTests(unittest.TestCase):
    """Test whether magres data can be loaded correctly."""

    def test_load_efg(self):
        from matador.data.magres import ELECTRIC_QUADRUPOLE_MOMENTS

        self.assertEqual(len(ELECTRIC_QUADRUPOLE_MOMENTS), 93)


class PeriodicTableTests(unittest.TestCase):
    """Test that the periodic table is written correctly."""

    def test_numbers(self):
        from matador.data.periodic_table import PERIODIC_TABLE

        for i in range(92):
            assert list(PERIODIC_TABLE.values())[i].number == i

    def test_keys_and_symbols(self):
        from matador.data.periodic_table import PERIODIC_TABLE

        for element in PERIODIC_TABLE:
            assert PERIODIC_TABLE[element].symbol == element
