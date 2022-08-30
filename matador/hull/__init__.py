# coding: utf-8
# Distributed under the terms of the MIT License.

""" The hull module provides functionality for creating phase diagrams,
voltage curves and volume curves, either directly from a database, or
from files.

"""


__all__ = [
    "QueryConvexHull",
    "PhaseDiagram",
    "EnsembleHull",
    "TemperatureDependentHull",
    "HullDiff",
    "diff_hulls",
]
__author__ = "Matthew Evans"
__maintainer__ = "Matthew Evans"


from .hull import QueryConvexHull
from .phase_diagram import PhaseDiagram
from .hull_ensemble import EnsembleHull
from .hull_temperature import TemperatureDependentHull
from .hull_diff import HullDiff, diff_hulls
