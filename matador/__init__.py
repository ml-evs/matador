# coding: utf-8
# Distributed under the terms of the MIT License.

""" `matador` performs, aggregates, and analyses high-throughput
electronic structure calculations, primarily for crystal structure
prediction, leveraging CASTEP and Quantum Epsresso as density-functional
theory compute engines.
"""

import sys


__all__ = ['__version__']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'
__version__ = "0.9.11"

script_epilog = f"Written and maintained by Matthew Evans (me388@cam.ac.uk) 2016-2021, version {__version__}."

if sys.version_info.minor == 6:
    # Python 3.6
    import warnings

    warnings.filterwarnings(
        action="once",
        message=r"v0\.9.x of the `matador-db` package.*",
        category=DeprecationWarning,
        append=False,
    )
    warnings.warn(
        "v0.9.x of the `matador-db` package will be the last to support Python 3.6. "
        "Please upgrade to Python 3.7+ to use v0.10 and later versions of `matador-db`.",
        DeprecationWarning,
    )
