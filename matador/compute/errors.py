# coding: utf-8
# Distributed under the terms of the MIT license.

""" This module implements some useful exceptions. """


class CalculationError(Exception):
    """ Raised when a particular calculation fails, for non-fatal reasons. """
    pass


class MaxMemoryEstimateExceeded(Exception):
    """ Raised when a structure is estimated to exceed the max memory. """
    pass


class CriticalError(RuntimeError):
    """ Raise this when you don't want any more jobs to run because something
    uncorrectable has happened! Plays more nicely with multiprocessing than
    SystemExit.

    """
    pass


class InputError(RuntimeError):
    """ Raise this when there is an issue with the input files. """
    pass


class WalltimeError(RuntimeError):
    """ Raise this when you don't want any more jobs to run
    because they're about to exceed the max walltime.

    """
    pass
