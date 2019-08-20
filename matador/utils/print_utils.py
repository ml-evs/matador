# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements some useful wrappers to the print function for
writing errors and warnings to stderr.

"""


import sys
import json
import numpy as np


def print_warning(string):
    """ Print but angry. """
    print('\033[93m', end='', file=sys.stderr)
    print(string, end='', file=sys.stderr)
    print('\033[0m', file=sys.stderr)


def print_failure(string):
    """ Print but sad. """
    print('\033[91m\033[4m', end='', file=sys.stderr)
    print(string, end='', file=sys.stderr)
    print('\033[0m', file=sys.stderr)


def print_success(string):
    """ Print but happy. """
    print('\033[92m\033[1m', end='')
    print(string, end='')
    print('\033[0m')


def print_notify(string):
    """ Print but aloof. """
    print('\033[94m', end='')
    print(string, end='')
    print('\033[0m')


def dumps(obj, **kwargs):
    """ Mirrors `json.dumps` whilst handling numpy arrays. """
    return json.dumps(obj, cls=NumpyEncoder, **kwargs)


class NumpyEncoder(json.JSONEncoder):
    """ This encoder handles NumPy arrays in JSON, and was
    taken from StackOverflow (where else)
    (https://stackoverflow.com/a/47626762).

    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
