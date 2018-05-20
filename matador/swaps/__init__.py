# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the AtomicSwapper class which
takes a Query object and a desired "swap", e.g. swap all group [I]
elements to Li, ready to be re-relaxed.

"""


__all__ = ['AtomicSwapper']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.swaps.swaps import AtomicSwapper
