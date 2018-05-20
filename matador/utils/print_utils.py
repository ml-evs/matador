# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements some useful wrappers to the print function for
writing errors and warnings to stderr.

"""


import sys


def print_warning(string):
    print('\033[93m', end='', file=sys.stderr)
    print(string, end='', file=sys.stderr)
    print('\033[0m', file=sys.stderr)


def print_failure(string):
    print('\033[91m\033[4m', end='', file=sys.stderr)
    print(string, end='', file=sys.stderr)
    print('\033[0m', file=sys.stderr)


def print_success(string):
    print('\033[92m\033[1m', end='')
    print(string, end='')
    print('\033[0m')


def print_notify(string):
    print('\033[94m', end='')
    print(string, end='')
    print('\033[0m')
