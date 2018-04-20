""" Define matador version. """
from pkg_resources import require, DistributionNotFound
try:
    __version__ = require('matador')[0].version
except DistributionNotFound:
    __version__ = 'xxx'
