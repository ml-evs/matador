# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file defines some useful scraper functionality,
like custom errors and a scraper function wrapper.
"""


import glob
import traceback as tb
import warnings

from matador.orm.spectral import VibrationalDOS

MODEL_REGISTRY = {
    "phonon_dos2dict": VibrationalDOS
}


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
            if '*' in seed and not kwargs.get('noglob'):
                seed = glob.glob(seed)
            else:
                seed = [seed]

        failures = []
        cursor = []

        if not seed:
            print('Nothing to scrape.')
            return

        for _seed in seed:
            # we can get away with this as each
            # scraper function only has one arg
            try:
                result, success = function(_seed, **kwargs)
            # UnicodeDecodeErrors require 5 arguments, so handle these separately
            except (FileNotFoundError, UnicodeError) as oops:
                raise oops
            except Exception as oops:
                success = False
                result = type(oops)('{}: {}\n'.format(_seed, oops))

                msg = '{}: {} {}'.format(_seed, type(oops), oops)
                if kwargs.get('verbosity', 1) > 0:
                    tb.print_exc()
                    print(msg)

            if len(seed) == 1:
                if success and not isinstance(result, dict):
                    raise AssertionError('Scraping succeeded, but dict not returned for {}'.format(seed))
                if not success and isinstance(result, dict):
                    raise AssertionError('Scraping failed, but dict returned for {}'.format(seed))

                return result, success

            if not success:
                failures += [_seed]
            else:

                if kwargs.get('as_model'):
                    model = MODEL_REGISTRY.get(function.__name__)
                    orm = None
                    if model is not None:
                        try:
                            orm = model(result)
                            cursor.append(orm)
                        except Exception:
                            print('Unable to convert scraped dict to model {}'.format(model.__name__))
                    else:
                        warnings.warn('as_model keyword not supported for {}'.format(function.__name__))

                if not kwargs.get('as_model') or orm is None:
                    cursor.append(result)

        return cursor, failures

    return wrapped_scraper_function


def f90_float_parse(val):
    """ Wrapper to float that handles Fortran's horrible behaviour for
    float exponents <= 100, e.g. 1e-100 -> 1.0000000-100 rather than
    1.000000E-100. Also handles "+" signs in front of numbers.

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
        if val.startswith('+'):
            val = val[1:]
        return float(val)


class DFTError(Exception):
    """ Quick DFT exception class for unconverged or
    non-useful calculations.

    """


class ComputationError(Exception):
    """ Raised when the calculation fails to do the DFT.
    Distinct from DFTError as this is an issue of numerics
    or chemistry, where this is raised for technical issues,
    e.g. CASTEP crashes.

    """
