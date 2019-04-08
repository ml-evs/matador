# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the base DataContainer class which wraps raw
matador dictionaries and adds useful methods to be inherited by its
children.

"""

import copy


class DataContainer:
    """ Base class for any data stored alongside a crystal. This class
    is a read-only store of the underlying dictionary of raw data; its
    children can implement useful methods to inspect and analyse the
    underlying data.

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
            self._data[key] = item
        else:
            raise AttributeError('Cannot assign value to existing key {}'.format(key))

    def __contains__(self, key):
        if key in self._data:
            return True
        return False

    def __str__(self):
        repr_string = "{root_source}: {formula}\n".format(root_source=self.root_source, formula=self.formula)
        repr_string += "{num_atoms:<3} atoms. {space_group:<8}\n".format(num_atoms=self.num_atoms,
                                                                         space_group=self.space_group)

        if 'formation_enthalpy_per_atom' in self._data:
            repr_string += ("Formation enthalpy = {:6.6f} eV/atom\n".format(self._data['formation_enthalpy_per_atom']))

        repr_string += (
            "(a, b, c) = {lattice[0][0]:4.4f} Å, {lattice[0][1]:4.4f} Å, {lattice[0][2]:4.4f} Å\n"
            "(α, β, γ) = {lattice[1][0]:4.4f}° {lattice[1][1]:4.4f}° {lattice[1][2]:4.4f}°\n"
            .format(lattice=self.lattice.abc))

        return repr_string

    def get(self, key):
        """ Overload dictionary.get() method.

        Parameters:
            key (str): key to try and obtain.

        Returns:
            doc[key] if it exists, else None.

        """
        return self._data.get(key)
