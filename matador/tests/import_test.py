#!/usr/bin/env python
import unittest


class ImportTest(unittest.TestCase):
    """ Test importing of all matador submodules. """
    def testQuery(self):
        from matador.query import DBQuery
        self.assertTrue(DBQuery)

    def testHull(self):
        from matador.hull import QueryConvexHull
        self.assertTrue(QueryConvexHull)

if __name__ == '__main__':
    unittest.main()
