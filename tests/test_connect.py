#!/usr/bin/env python
import unittest


class ConnectTest(unittest.TestCase):
    """Test simple database connection methods."""

    def testFuzzyCollnames(self):
        """Test the fuzzy collection name matcher."""
        from matador.db.connect import fuzzy_collname_match

        trial = "me388_proto_test"
        targets = [
            "me388_prototype_test",
            "proto_test_me388",
            "me388_fun",
            "ajm255_test",
            "ajm255_no",
        ]
        options = fuzzy_collname_match(trial, targets)

        results = [targets[0], targets[1], targets[2], targets[3]]
        self.assertEqual(results, options)


if __name__ == "__main__":
    unittest.main()
