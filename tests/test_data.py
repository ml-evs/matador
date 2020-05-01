import unittest


class ScatteringDataLoaderTests(unittest.TestCase):
    """ Test whether the atomic scattering factors can be loaded. """

    def test_load_gsas(self):
        from matador.data import GSAS_ATOMIC_SCATTERING_COEFFS

        self.assertEqual(len(GSAS_ATOMIC_SCATTERING_COEFFS["Cu"]), 3)
        for key in GSAS_ATOMIC_SCATTERING_COEFFS:
            self.assertEqual(len(GSAS_ATOMIC_SCATTERING_COEFFS[key]), 3)

    def test_load_raspa(self):
        from matador.data import RASPA_ATOMIC_SCATTERING_COEFFS

        self.assertEqual(len(RASPA_ATOMIC_SCATTERING_COEFFS["Cu"]), 3)
        for key in RASPA_ATOMIC_SCATTERING_COEFFS:
            self.assertEqual(len(RASPA_ATOMIC_SCATTERING_COEFFS[key]), 3)
