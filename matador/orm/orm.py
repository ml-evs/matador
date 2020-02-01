# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the base DataContainer class which wraps raw
matador dictionaries and adds useful methods to be inherited by its
children.

"""

import copy
from abc import ABC


class DataContainer(ABC):
    """ Base class for matador data classes. This class is a read-only
    store of the underlying dictionary of raw data; its children can
    implement useful methods to inspect and analyse the underlying data.

    """
    def __init__(self, data):
        """ Initalise copy of raw data. """
        self._data = copy.deepcopy(data)

    def __getitem__(self, key):
        """ Allow properties to be used with key access.

        Parameters:
            key (str): name of key or attribute to get.

        Raises:
            AttributeError: if key or attribute can't be found.

        """
        try:
            return getattr(self, key)
        except AttributeError:
            pass

        try:
            return self._data[key]
        except KeyError:
            raise AttributeError('Object has no data or implementation for requested key: "{}"'
                                 .format(key))

    def __delitem__(self, key: str):
        raise AttributeError('Object does not support deletion of keys in `_data`.')

    def __setitem__(self, key: str, item):
        if key not in self._data:
            self._data[key] = item
        else:
            raise AttributeError('Cannot assign value to existing key {}'
                                 .format(key))

    def __contains__(self, key):
        if key in self._data:
            return True
        return False

    def get(self, *args):
        """ Overload dictionary.get() method.

        Parameters:
            key (str): key to try and obtain.

        Keyword arguments:
            default: return value raise if key is not present.

        Returns:
            :attr:`_data[key]` if it exists, else None.

        """
        key = args[0]
        if len(args) > 1:
            default = args[1]
        else:
            default = None
        if len(args) > 2:
            raise TypeError('get() takes up to 2 arguments, not {}'.format(len(args)))

        return self._data.get(key, default)
