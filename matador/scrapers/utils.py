# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file defines some useful scraper functionality,
like custom errors and a scraper function wrapper.
"""


import glob
import os
import gzip
import traceback as tb

from matador.orm.spectral import VibrationalDOS, VibrationalDispersion
from matador.orm.spectral import ElectronicDOS, ElectronicDispersion
from matador.crystal import Crystal

MODEL_REGISTRY = {
    "phonon_dos2dict": VibrationalDOS,
    "phonon2dict": VibrationalDispersion,
    "optados2dict": ElectronicDOS,
    "bands2dict": ElectronicDispersion,
    "castep2dict": Crystal,
    "res2dict": Crystal,
    "cif2dict": Crystal,
    "magres2dict": Crystal,
    "pwout2dict": Crystal,
}


def get_flines_extension_agnostic(fname, ext):
    """Try to open and read the filename provided, if it doesn't exist
    then try adding the given file extension to it.

    Parameters:
        fname (str): the filename with or without extension.
        ext (list of str or str): the extension or list of file extensions to try,
            or None. Should not contain ".".

    Raises:
        FileNotFoundError: if the file was not found in either form.

    Returns:
        (list of str, str): the contents of the file and the filename.

    """
    if isinstance(ext, str):
        ext = [ext]

    if ext is not None and not os.path.isfile(fname):
        for exts in ext:
            if not fname.endswith(exts):
                _fname = f"{fname}.{exts}"
                if os.path.isfile(_fname):
                    fname = _fname
                    break

    try:
        if fname.endswith(".gz"):
            with gzip.open(fname, "r") as f:
                flines = [line.decode("utf-8") for line in f.readlines()]
        else:
            try:
                with open(fname, "r", encoding="utf-8") as f:
                    flines = f.readlines()
            except Exception:
                with open(fname, "r", encoding="latin1") as f:
                    flines = f.readlines()

    except FileNotFoundError as exc:
        if ext is not None:
            raise FileNotFoundError(f"Neither {fname} or {fname}.{ext} could be found.")
        raise exc

    return flines, fname


def scraper_function(function):
    """Wrapper for scraper functions to handle exceptions and
    template the scraper functions to work for multiples files
    at once.

    """

    from functools import wraps

    @wraps(function)
    def wrapped_scraper_function(*args, verbosity=1, fail_fast=False, **kwargs):
        """Wrap and return the scraper function, handling the
        multiplicity of file names.

        """

        if kwargs.get("no_wrap"):
            return function(*args, **kwargs)

        result = None
        seed = args[0]
        if isinstance(seed, str):
            if "*" in seed and not kwargs.get("noglob"):
                seed = sorted(glob.glob(seed))
            else:
                seed = [seed]

        failures = []
        cursor = []

        if not seed:
            print("Nothing to scrape.")
            return

        for _seed in seed:
            # we can get away with this as each
            # scraper function only has one arg
            try:
                result, success = function(_seed, verbosity=verbosity, **kwargs)
            # UnicodeDecodeErrors require 5 arguments, so handle these separately
            except (FileNotFoundError, UnicodeError) as oops:
                raise oops
            except Exception as oops:
                success = False
                result = type(oops)("{}: {}\n".format(_seed, oops))

                if verbosity >= 1:
                    msg = "{}: {} {}".format(_seed, type(oops), oops)
                    print(msg)
                if verbosity >= 2:
                    tb.print_exc()

                if fail_fast:
                    raise oops

            if len(seed) == 1:
                if success and kwargs.get("as_model"):
                    orm = _as_model(result, function)
                    if orm is not None:
                        result = orm
                if not success and verbosity >= 1:
                    print("Failed to scrape file {}".format(seed))

                return result, success

            if not success:
                failures += [_seed]
            else:
                if kwargs.get("as_model"):
                    orm = _as_model(result, function, debug=kwargs.get("debug"))
                    cursor.append(orm)
                if not kwargs.get("as_model") or orm is None:
                    cursor.append(result)

        if verbosity >= 1:
            print(
                "Successfully scraped {} out of {} files.".format(
                    len(cursor), len(cursor) + len(failures)
                )
            )

        return cursor, failures

    return wrapped_scraper_function


def _as_model(doc, function, debug=True):
    """Convert the document to the appropriate orm model."""
    model = MODEL_REGISTRY.get(function.__name__)
    orm = None
    if model is not None:
        try:
            orm = model(doc)
        except Exception as exc:
            if debug:
                tb.print_exc()
            print("Unable to convert scraped dict to model {}".format(model.__name__))
            raise exc
    else:
        print(
            "`as_model` keyword not supported for {}, not converting".format(
                function.__name__
            )
        )

    return orm


def f90_float_parse(val):
    """Wrapper to float that handles Fortran's horrible behaviour for
    float exponents <= 100, e.g. 1e-100 -> 1.0000000-100 rather than
    1.000000E-100. Also handles "+" signs in front of numbers.

    Parameters:
        val (str): the string to cast to a float.

    """
    try:
        return float(val)
    except ValueError as exc:
        # if the E is being displayed, then something else has gone wrong
        if "E" in val:
            raise exc
        # if there's a minus sign after the first char, but no E...
        if len(val) > 1 and "-" in val[1:]:
            val = val[0] + val[1:].replace("-", "E-")
        if val.startswith("+"):
            val = val[1:]
        return float(val)


class DFTError(Exception):
    """Quick DFT exception class for unconverged or
    non-useful calculations.

    """


class ComputationError(Exception):
    """Raised when the calculation fails to do the DFT.
    Distinct from DFTError as this is an issue of numerics
    or chemistry, where this is raised for technical issues,
    e.g. CASTEP crashes.

    """
