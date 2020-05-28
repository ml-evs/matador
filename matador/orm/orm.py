# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the base DataContainer class which wraps raw
matador dictionaries and adds useful methods to be inherited by its
children.

"""

import copy
import math
from abc import ABC


class DataContainer(ABC):
    """ Base class for matador data classes. This class is a read-only
    store of the underlying dictionary of raw data; its children can
    implement useful methods to inspect and analyse the underlying data.

    """

    required_keys = []

    def __init__(self, data=None, **kwargs):
        """ Initalise copy of raw data. """
        if isinstance(data, dict):
            self._data = copy.deepcopy(data)
        else:
            self._data = {key: kwargs[key] for key in kwargs}

        self._validate_inputs()

        # set root source to source filename
        from matador.utils.chem_utils import get_root_source
        self.root_source = 'unknown'
        try:
            if 'source' in self._data:
                self.root_source = get_root_source(self._data['source'])
        except RuntimeError:
            pass

    @property
    def source(self):
        """ Return the source of the data. """
        if 'source' not in self._data:
            return 'unknown'
        return self._data['source']

    def _validate_inputs(self):
        """ Validate the incoming data by checking the existence
        of the required keys for the subclass.

        """
        missing_keys = []
        for key in self.required_keys:
            if key not in self._data:
                missing_keys.append(key)

        if missing_keys:
            raise RuntimeError(
                "Unable to create object of type {} as the following keys are missing: {}"
                .format(type(self).__name__, missing_keys)
            )

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
            raise KeyError('Object has no data or implementation for requested key: "{}"'
                           .format(key))

    def __delitem__(self, key: str):
        raise AttributeError('Object does not support deletion of keys in `_data`.')

    def __setitem__(self, key: str, item):
        if key not in self._data or self._data[key] is None:
            self._data[key] = item
            return
        elif self._data[key] != item:
            try:
                if (math.isnan(item) and math.isnan(self._data[key])):
                    return
            except TypeError:
                pass
            raise AttributeError('Cannot assign value {} to existing key {} with value {}'
                                 .format(item, key, self._data[key]))

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
