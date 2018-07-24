# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file defines some useful scraper functionality,
like custom errors and a scraper function wrapper.

"""


from traceback import print_exc
import json


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
        if isinstance(seed, list):
            failures = []
            cursor = []
            for _seed in seed:
                # we can get away with this as each
                # scraper function only has one arg
                try:
                    result, success = function(_seed, **kwargs)
                except FileNotFoundError as oops:
                    raise oops
                except Exception as oops:
                    success = False
                    if kwargs.get('verbosity', 0) > 0:
                        print_exc()
                        print(oops)
                        print('Error in file', _seed, 'skipping...')

                if not success:
                    failures += [_seed]
                else:
                    cursor.append(result)
            if len(cursor) != len(seed):
                print('Scraped {}/{} structures.'.format(len(cursor), len(seed)))

            return cursor, failures

        try:
            result = function(*args, **kwargs)
        except FileNotFoundError as oops:
            raise oops
        except Exception as oops:
            if kwargs.get('verbosity', 0) > 0:
                print_exc()
                print('Error in file', seed, 'skipping...')
            if kwargs.get('debug'):
                raise oops
            return '{} {}: {}\n'.format(seed, type(oops), str(oops)), False

        if kwargs.get('verbosity', 0) > 4:
            try:
                print(json.dumps(result[0], indent=2))
            except TypeError:
                pass

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
