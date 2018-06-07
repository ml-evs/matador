# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file defines some useful scraper functionality,
like custom errors and a scraper function wrapper.

"""


from traceback import print_exc
import json


def scraper_function(function):
    """ Wrapper for scraper functions to handle exceptions. """

    from functools import wraps

    @wraps(function)
    def wrapped_scraper_function(*args, **kwargs):
        """ Wrap and return the plotting function. """
        result = None
        seed = args[0]
        try:
            result = function(*args, **kwargs)
        except Exception as oops:
            if kwargs.get('verbosity', 0) > 0:
                print_exc()
                print('Error in file', seed, 'skipping...')
            return '{} {}: {}'.format(seed, type(oops), str(oops)), False

        if kwargs.get('verbosity', 0) > 4:
            print(json.dumps(result[0], indent=2))

        return result

    return wrapped_scraper_function


class DFTError(Exception):
    """ Quick DFT exception class for unconverged or
    non-useful calculations.

    """
    pass


class CalculationError(Exception):
    """ Raised when the calculation fails to do the DFT.
    Distinct from DFTError as this is an issue of numerics
    or chemistry, where this is raised for technical issues,
    e.g. CASTEP crashes.

    """
    pass
