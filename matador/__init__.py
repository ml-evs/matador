# coding: utf-8
# Distributed under the terms of the MIT License.

""" `matador` performs, aggregates, and analyses high-throughput
electronic structure calculations, primarily for crystal structure
prediction, leveraging CASTEP and Quantum Epsresso as density-functional
theory compute engines.
"""

import sys

__all__ = ["__version__"]
__author__ = "Matthew Evans"
__maintainer__ = "Matthew Evans"
__version__ = "0.10.1"

script_epilog = f"Written and maintained by Matthew Evans (me388@cam.ac.uk) 2016-2022, version {__version__}."


if sys.version_info.minor == 7:
    import warnings

    warnings.warn(
        "v0.10 of the `matador` library will be the last to support Python 3.7. "
        "Please upgrade to Python 3.8+ to use v0.11 and later versions of `matador`.",
        DeprecationWarning,
    )
