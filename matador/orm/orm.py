# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the base DataContainer class which wraps raw
matador dictionaries and adds useful methods to be inherited by its
children.

"""

import copy


class DataContainer:
    """ Base class for matador data classes. This class is a read-only
    store of the underlying dictionary of raw data; its children can
    implement useful methods to inspect and analyse the underlying data.

    """
    def __init__(self, data):
        """ Initalise copy of raw data. """
        self._data = copy.deepcopy(data)

    def __getitem__(self, key):
        """ Allow properties to be used with key access.

        Raises a KeyError if property or key can't be found.

        """
        try:
            return getattr(self, key)
        except AttributeError:
            pass

        try:
            return self._data[key]
        except KeyError:
            raise AttributeError('Object has no data or implementation for requested {}'
                                 .format(key))

    def __setitem__(self, key, item):
        if key not in self._data:
            raise AttributeError('Cannot set value of {} inside DataContainer.'
                                 .format(key))
            # self._data[key] = item
        else:
            raise AttributeError('Cannot assign value to existing key {}'
                                 .format(key))

    def __contains__(self, key):
        if key in self._data:
            return True
        return False

    def get(self, key):
        """ Overload dictionary.get() method.

        Parameters:
            key (str): key to try and obtain.

        Returns:
            :attr:`_data[key]` if it exists, else None.

        """
        return self._data.get(key)