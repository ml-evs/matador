# coding: utf-8
# Distributed under the terms of the MIT License.

""" The fingerprints module provides functionality for calculating, in
parallel, structural fingerprints like PDF and PXRD.

"""


__all__ = ['get_uniq_cursor',
           'Fingerprint', 'FingerprintFactory',
           'PDF', 'PDFOverlap', 'CombinedProjectedPDF',
           'PXRD']

__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.fingerprints.similarity import get_uniq_cursor
from matador.fingerprints.fingerprint import Fingerprint, FingerprintFactory
from matador.fingerprints.pdf import PDF, PDFOverlap, CombinedProjectedPDF
from matador.fingerprints.pxrd import PXRD
