""" This file implements some useful
wrappers to the print function.
"""
from __future__ import print_function


def print_warning(string):
    print('\033[93m', end='')
    print(string, end='')
    print('\033[0m')


def print_failure(string):
    print('\033[91m\033[4m', end='')
    print(string, end='')
    print('\033[0m')


def print_success(string):
    print('\033[92m\033[1m', end='')
    print(string, end='')
    print('\033[0m')


def print_notify(string):
    print('\033[94m', end='')
    print(string, end='')
    print('\033[0m')
