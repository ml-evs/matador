# coding: utf-8
# Distributed under the terms of the MIT License.

""" The export module provides interfaces to write matador-created
dictionaries/Crystal objects to several file formats.

Currently supported formats:

    - CASTEP's .cell and .param
    - the Protein Data Bank's .pdb
    - input files for Quantum Espresso (written as .in)
    - XCrysden's .xsf
    - the custom .res file format based on SHELX used first by AIRSS
        (https://www.mtg.msm.cam.ac.uk/Codes/AIRSS)

"""


__all__ = [
    "doc2param",
    "doc2cell",
    "doc2pdb",
    "doc2pwscf",
    "doc2res",
    "doc2xsf",
    "query2files",
    "doc2arbitrary",
]
__author__ = "Matthew Evans"
__maintainer__ = "Matthew Evans"


from .export import doc2param, doc2cell, doc2pdb
from .export import doc2pwscf, doc2res, doc2xsf, doc2arbitrary
from .export import query2files
