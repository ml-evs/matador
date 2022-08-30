# coding: utf-8
# Distributed under the terms of the MIT License.

""" The query module provides the DBQuery class that implements all
database queries, either for individual compositions/formulae, or for
matching calculations at the same accuracy to create phase diagrams.

"""


__all__ = ["DBQuery"]
__author__ = "Matthew Evans"
__maintainer__ = "Matthew Evans"


from .query import DBQuery
