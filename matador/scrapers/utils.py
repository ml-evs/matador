# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file defines some useful scraper functionality,
like custom errors and a scraper function wrapper.

"""


import glob


def scraper_function(function):
    """ Wrapper for scraper functions to handle exceptions and
    template the scraper functions to work for multiples files
    at once.

    """

    from functools import wraps

    @wraps(function)
    def wrapped_scraper_function(*args, **kwargs):
        """ Wrap and return the scraper function, handling the
        multiplicity of file names.

        """

        result = None
        seed = args[0]
        if isinstance(seed, str):
            if '*' in seed:
                seed = glob.glob(seed)
            else:
                seed = [seed]

        failures = []
        cursor = []

        if not seed:
            raise RuntimeError('Found nothing to scrape.')
        else:
            for _seed in seed:
                # we can get away with this as each
                # scraper function only has one arg
                try:
                    result, success = function(_seed, **kwargs)
                except FileNotFoundError as oops:
                    raise oops
                except Exception as oops:
                    success = False
                    msg = '{}: {} {}'.format(_seed, type(oops), oops)
                    print(msg)
                    result = type(oops)('{}: {}'.format(_seed, oops))

                if not success:
                    failures += [_seed]
                else:
                    cursor.append(result)

            if len(seed) == 1:
                return result, success

            return cursor, failures

    return wrapped_scraper_function


def f90_float_parse(val):
    """ Wrapper to float that handles Fortran's horrible behaviour for
    float exponents <= 100, e.g. 1e-100 -> 1.0000000-100 rather than
    1.000000E-100.

    Parameters:
        val (str): the string to cast to a float.

    """
    try:
        return float(val)
    except ValueError as exc:
        # if the E is being displayed, then something else has gone wrong
        if 'E' in val:
            raise exc
        # if there's a minus sign after the first char, but no E...
        if len(val) > 1 and '-' in val[1:]:
            val = val[0] + val[1:].replace('-', 'E-')
        return float(val)


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
