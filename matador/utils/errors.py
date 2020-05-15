# coding: utf-8
# Distributed under the terms of the MIT license.

""" This submodule module implements some useful exception types,
mostly for use in the :mod:`compute` and :mod:`calculator` submodules.

"""


class CalculationError(Exception):
    """ Raised when a particular calculation fails, for non-fatal reasons. """


class MaxMemoryEstimateExceeded(Exception):
    """ Raised when a structure is estimated to exceed the max memory. """


class CriticalError(RuntimeError):
    """ Raise this when you don't want any more jobs to run because something
    uncorrectable has happened! Plays more nicely with multiprocessing than
    SystemExit.

    """


class InputError(RuntimeError):
    """ Raise this when there is an issue with the input files. """


class WalltimeError(RuntimeError):
    """ Raise this when you don't want any more jobs to run
    because they're about to exceed the max walltime.

    """


class NodeCollisionError(CalculationError):
    """ Dummy exception to raise when one node has tried to run
    a calculation that another node is performing.

    """
