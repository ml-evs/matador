#!/usr/bin/env python
import unittest
from matador.fingerprints.pdf import PDF, PDFOverlap, PDFFactory
from matador.scrapers.castep_scrapers import res2dict
from matador.utils.cell_utils import abc2cart, cart2volume
import numpy as np
from os.path import realpath

REAL_PATH = "/".join(realpath(__file__).split("/")[:-1]) + "/"
DEBUG = True


class PDFCalculatorTest(unittest.TestCase):
    """ Test PDF calculator. """

    def test_auto_images_vs_large(self):
        doc, success = res2dict(REAL_PATH + "data/LiPZn-r57des.res")
        doc["lattice_cart"] = abc2cart(doc["lattice_abc"])
        doc["text_id"] = ["pdf", "test"]
        doc["pdf_num_images"] = PDF(
            doc,
            low_mem=True,
            num_images=5,
            rmax=15,
            dr=0.1,
            **{"debug": True, "projected": False}
        )
        doc["pdf_auto_images"] = PDF(
            doc,
            low_mem=True,
            num_images="auto",
            rmax=15,
            dr=0.1,
            **{"debug": True, "projected": False}
        )
        np.testing.assert_array_almost_equal(
            doc["pdf_num_images"].gr, doc["pdf_auto_images"].gr
        )

    def test_pdf_from_projected(self):
        doc, success = res2dict(REAL_PATH + "data/LiPZn-r57des.res")
        doc["lattice_cart"] = abc2cart(doc["lattice_abc"])
        doc["text_id"] = ["unprojected", "test"]
        doc["pdf_unprojected"] = PDF(doc, dr=0.01, **{"debug": False})
        doc["text_id"] = ["projected", "test"]
        doc["pdf_projected"] = PDF(doc, dr=0.01, **{"debug": False})
        np.testing.assert_array_almost_equal(
            doc["pdf_unprojected"].gr, doc["pdf_projected"].gr
        )

    def test_identity_overlap(self):
        doc, success = res2dict(REAL_PATH + "data/LiPZn-r57des.res")
        doc["lattice_cart"] = abc2cart(doc["lattice_abc"])
        doc["text_id"] = ["pdf", "test"]
        doc["pdf_smear"] = PDF(
            doc,
            num_images=3,
            dr=0.001,
            gaussian_width=0.1,
            style="smear",
            debug=False,
            low_mem=True,
        )
        overlap = PDFOverlap(doc["pdf_smear"], doc["pdf_smear"])
        self.assertEqual(overlap.similarity_distance, 0.0)

    @unittest.skip
    def test_overlap_smear_vs_hist(self):
        doc, success = res2dict(REAL_PATH + "data/LiPZn-r57des.res")
        doc["lattice_cart"] = abc2cart(doc["lattice_abc"])
        doc["text_id"] = ["smear", "test"]
        doc["pdf_smear"] = PDF(
            doc,
            num_images=3,
            dr=0.01,
            gaussian_width=0.1,
            projected=False,
            style="smear",
            low_mem=True,
        )
        doc["text_id"] = ["hist", "test"]
        doc["pdf_hist"] = PDF(
            doc, num_images=3, dr=0.1, projected=False, style="histogram"
        )
        overlap = PDFOverlap(doc["pdf_smear"], doc["pdf_hist"])
        self.assertLessEqual(overlap.similarity_distance, 0.02)
        self.assertGreater(overlap.similarity_distance, 0.0)

    def test_pdf_primitive_vs_supercell(self):
        test_doc, success = res2dict(REAL_PATH + "data/KP_primitive.res", db=False)
        test_doc["text_id"] = ["primitive", "cell"]
        test_doc["lattice_cart"] = abc2cart(test_doc["lattice_abc"])
        test_doc["cell_volume"] = cart2volume(test_doc["lattice_cart"])
        supercell_doc, success = res2dict(REAL_PATH + "data/KP_supercell.res", db=False)
        supercell_doc["text_id"] = ["supercell", "cell"]
        supercell_doc["lattice_cart"] = abc2cart(supercell_doc["lattice_abc"])
        supercell_doc["cell_volume"] = cart2volume(supercell_doc["lattice_cart"])
        test_doc["pdf"] = PDF(
            test_doc, dr=0.01, low_mem=True, rmax=10, num_images="auto", debug=DEBUG
        )
        supercell_doc["pdf"] = PDF(
            supercell_doc,
            dr=0.01,
            low_mem=True,
            rmax=10,
            num_images="auto",
            debug=DEBUG,
        )
        overlap = PDFOverlap(test_doc["pdf"], supercell_doc["pdf"])
        self.assertLessEqual(overlap.similarity_distance, 1e-3)
        self.assertGreaterEqual(overlap.similarity_distance, 0.0)

    @unittest.skip("Not sure this test is worthwhile anymore...")
    def test_ideal_gas_pdf(self, retry=0):
        """ DEPRECATED.

        Slow, and not very useful.

        """
        # create fake matador doc
        doc = dict()
        max_retries = 1
        self.assertLess(
            retry, max_retries, msg="After {} attempts, PDF still failed.".format(retry)
        )
        num_atoms = 100
        box_size = 10
        num_samples = 1
        rmax = 15
        dr = 0.01
        num_images = "auto"
        i = 0
        doc["atom_types"] = num_atoms * ["C"]
        doc["lattice_cart"] = np.asarray(
            [[box_size, 0, 0], [0, box_size, 0], [0, 0, box_size]]
        )
        doc["cell_volume"] = box_size ** 3
        doc["text_id"] = ["ideal", "gas"]
        while i < num_samples:
            doc["positions_frac"] = np.random.rand(num_atoms, 3)
            doc["text_id"] = "hist"
            doc["pdf"] = PDF(
                doc,
                num_images=num_images,
                dr=dr,
                rmax=rmax,
                lazy=True,
                style="histogram",
                debug=True,
            )
            doc["pdf"].calc_pdf()
            self.assertAlmostEqual(np.mean(doc["pdf"].gr[50:]), 1.0, places=1)

    def test_single_atom_pdf(self):
        from math import ceil

        doc = dict()
        box_size = 20
        rmax = 41
        dr = 0.1
        num_images = "auto"
        doc["positions_frac"] = [[0.5, 0.5, 0.5]]
        doc["atom_types"] = ["C"]
        doc["lattice_cart"] = np.asarray(
            [[box_size, 0, 0], [0, box_size, 0], [0, 0, box_size]]
        )
        doc["cell_volume"] = box_size ** 3
        doc["text_id"] = ["hist", "ogram"]
        doc["pdf"] = PDF(
            doc,
            num_images=num_images,
            dr=dr,
            rmax=rmax,
            lazy=True,
            style="histogram",
            debug=True,
        )
        doc["pdf"].calc_pdf()
        doc["text_id"] = ["smear"]
        doc["pdf_smear"] = PDF(
            doc,
            num_images=num_images,
            gaussian_width=0.01,
            dr=0.1,
            rmax=rmax,
            lazy=True,
            style="smear",
            debug=True,
        )
        doc["pdf_smear"].calc_pdf()
        doc["text_id"] = ["low"]
        doc["pdf_low"] = PDF(
            doc,
            low_mem=True,
            num_images=num_images,
            gaussian_width=0.01,
            dr=0.1,
            rmax=rmax,
            lazy=True,
            style="smear",
            debug=True,
        )
        doc["pdf_low"].calc_pdf()
        peaks = [20, np.sqrt(2) * 20, np.sqrt(3) * 20, 40]
        indices = [ceil(peak / dr) for peak in peaks]
        self.assertListEqual(np.where(doc["pdf_low"].gr > 1e-8)[0].tolist(), indices)
        self.assertListEqual(np.where(doc["pdf"].gr > 1e-8)[0].tolist(), indices)
        self.assertListEqual(np.where(doc["pdf_smear"].gr > 1e-8)[0].tolist(), indices)


class PDFFactoryCalculatorTest(unittest.TestCase):
    """ Test PDFFactory. """

    def test_concurrent_pdf(self):
        import glob
        import numpy as np
        import time
        from copy import deepcopy

        files = glob.glob(REAL_PATH + "data/hull-KPSn-KP/*.res")[0:24]
        cursor = [res2dict(file, db=False)[0] for file in files]
        serial_cursor = deepcopy(cursor)
        pdf_args = {
            "dr": 0.1,
            "num_images": "auto",
            "gaussian_width": 0.1,
            "lazy": False,
            "projected": False,
        }
        start = time.time()
        pdf_factory = PDFFactory(cursor, **pdf_args)
        factory_elapsed = time.time() - start
        start = time.time()
        for doc in serial_cursor:
            doc["pdf"] = PDF(doc, **pdf_args, timing=False)
        serial_elapsed = time.time() - start
        print(
            "{:.2f} s over {} processes vs {:.2f} s in serial".format(
                factory_elapsed, pdf_factory.nprocs, serial_elapsed
            )
        )
        print(
            "Corresponding to a speedup of {:.1f} vs ideal {:.1f}".format(
                serial_elapsed / factory_elapsed, pdf_factory.nprocs
            )
        )
        for ind, doc in enumerate(serial_cursor):
            np.testing.assert_array_almost_equal(
                doc["pdf"].gr, cursor[ind]["pdf"].gr, decimal=6
            )

    def test_concurrent_pdf_stoichs(self):
        import glob
        import numpy as np
        import time
        from copy import deepcopy
        from matador.hull import QueryConvexHull

        files = glob.glob(REAL_PATH + "data/hull-KPSn-KP/*.res")
        cursor = [res2dict(file, db=True)[0] for file in files]
        hull = QueryConvexHull(
            cursor=cursor,
            no_plot=True,
            hull_cutoff=0.5,
            summary=True,
            elements=["K", "Sn", "P"],
            quiet=True,
        )
        serial_cursor = deepcopy(hull.cursor)

        pdf_args = {
            "dr": 0.1,
            "num_images": "auto",
            "gaussian_width": 0.1,
            "lazy": False,
            "projected": False,
        }
        start = time.time()
        pdf_factory = PDFFactory(hull.cursor, **pdf_args)
        factory_elapsed = time.time() - start
        start = time.time()
        for doc in serial_cursor:
            doc["pdf"] = PDF(doc, **pdf_args, timing=False)
        serial_elapsed = time.time() - start
        print(
            "{:.2f} s over {} processes vs {:.2f} s in serial".format(
                factory_elapsed, pdf_factory.nprocs, serial_elapsed
            )
        )
        print(
            "Corresponding to a speedup of {:.1f} vs ideal {:.1f}".format(
                serial_elapsed / factory_elapsed, pdf_factory.nprocs
            )
        )
        for ind, doc in enumerate(serial_cursor):
            np.testing.assert_array_almost_equal(
                doc["pdf"].gr, hull.cursor[ind]["pdf"].gr, decimal=6
            )

    def test_concurrent_projected_pdf(self):
        import glob
        import numpy as np
        import time
        from copy import deepcopy

        files = glob.glob(REAL_PATH + "data/hull-KPSn-KP/*.res")[0:24]
        cursor = [res2dict(file, db=False)[0] for file in files]
        serial_cursor = deepcopy(cursor)
        pdf_args = {
            "dr": 0.1,
            "num_images": "auto",
            "gaussian_width": 0.1,
            "lazy": False,
            "projected": True,
        }
        start = time.time()
        pdf_factory = PDFFactory(cursor, **pdf_args)
        factory_elapsed = time.time() - start
        start = time.time()
        for doc in serial_cursor:
            doc["pdf"] = PDF(doc, **pdf_args, timing=False)
        serial_elapsed = time.time() - start
        print(
            "{:.2f} s over {} processes vs {:.2f} s in serial".format(
                factory_elapsed, pdf_factory.nprocs, serial_elapsed
            )
        )
        print(
            "Corresponding to a speedup of {:.1f} vs ideal {:.1f}".format(
                serial_elapsed / factory_elapsed, pdf_factory.nprocs
            )
        )
        for ind, doc in enumerate(serial_cursor):
            np.testing.assert_array_almost_equal(
                doc["pdf"].gr, cursor[ind]["pdf"].gr, decimal=6
            )
            for key in doc["pdf"].elem_gr:
                np.testing.assert_array_almost_equal(
                    doc["pdf"].elem_gr[key], cursor[ind]["pdf"].elem_gr[key], decimal=6
                )


if __name__ == "__main__":
    unittest.main()
