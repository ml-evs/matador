# coding: utf-8
# Distributed under the terms of the MIT License.

""" `matador` performs, aggregates, and analyses high-throughput
electronic structure calculations, primarily for crystal structure
prediction, leveraging CASTEP and Quantum Epsresso as density-functional
theory compute engines.
"""


__all__ = ['__version__', 'DBQuery', 'QueryConvexHull', 'Crystal']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from pkg_resources import require, DistributionNotFound
try:
    __version__ = require('matador')[0].version
except DistributionNotFound:
    __version__ = 'xxx'

from matador.query import DBQuery
from matador.hull import QueryConvexHull
from matador.crystal import Crystal
