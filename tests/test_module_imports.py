#!/usr/bin/env python
import unittest


class ImportTest(unittest.TestCase):
    """ Test importing of all matador submodules. """

    def test_battery(self):
        from matador.battery import Electrode

        self.assertTrue(Electrode)

    def test_calculators(self):
        from matador.calculators import CastepCalculator

        self.assertTrue(CastepCalculator)

    def test_compute(self):
        from matador.compute import ComputeTask, BatchRun, reset_job_folder

        self.assertTrue(ComputeTask)
        self.assertTrue(BatchRun)
        self.assertTrue(reset_job_folder)

    def test_config(self):
        from matador.config import load_custom_settings

        self.assertTrue(load_custom_settings)

    def test_crystal(self):
        from matador.crystal import Crystal

        self.assertTrue(Crystal)

    def test_db(self):
        from matador.db import (
            make_connection_to_collection,
            Spatula,
            DatabaseChanges,
            Refiner,
        )

        self.assertTrue(make_connection_to_collection)
        self.assertTrue(Spatula)
        self.assertTrue(DatabaseChanges)
        self.assertTrue(Refiner)

    def test_export(self):
        from matador.export import (
            doc2param,
            doc2cell,
            doc2pdb,
            doc2pwscf,
            doc2res,
            doc2xsf,
            query2files,
            doc2arbitrary,
        )

        self.assertTrue(doc2param)
        self.assertTrue(doc2cell)
        self.assertTrue(doc2pdb)
        self.assertTrue(doc2pwscf)
        self.assertTrue(doc2res)
        self.assertTrue(doc2xsf)
        self.assertTrue(query2files)
        self.assertTrue(doc2arbitrary)

    def test_fingerprints(self):
        from matador.fingerprints import (
            PDF,
            PDFOverlap,
            CombinedProjectedPDF,
            PXRD,
            Fingerprint,
            FingerprintFactory,
            get_uniq_cursor,
        )

        self.assertTrue(
            all(
                [
                    PDF,
                    PDFOverlap,
                    CombinedProjectedPDF,
                    PXRD,
                    Fingerprint,
                    FingerprintFactory,
                    get_uniq_cursor,
                ]
            )
        )

    def test_hull(self):
        from matador.hull import (
            QueryConvexHull,
            EnsembleHull,
            PhaseDiagram,
            HullDiff,
            diff_hulls,
        )

        self.assertTrue(QueryConvexHull)
        self.assertTrue(EnsembleHull)
        self.assertTrue(PhaseDiagram)
        self.assertTrue(HullDiff)
        self.assertTrue(diff_hulls)

    def test_query(self):
        from matador.query import DBQuery

        self.assertTrue(DBQuery)

    def test_swaps(self):
        from matador.swaps import AtomicSwapper

        self.assertTrue(AtomicSwapper)


if __name__ == "__main__":
    unittest.main()
