#!/usr/bin/env python
import unittest


class ImportTest(unittest.TestCase):
    """ Test importing of all matador submodules. """
    def test_all_imports(self):
        import matador.battery
        import matador.calculators
        import matador.compute
        import matador.config
        import matador.crystal
        import matador.db
        import matador.export
        import matador.fingerprints
        import matador.hull
        import matador.orm
        import matador.plotting
        import matador.scrapers
        import matador.query
        import matador.swaps


if __name__ == '__main__':
    unittest.main()
