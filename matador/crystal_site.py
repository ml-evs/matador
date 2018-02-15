#!/usr/bin/env python
# coding: utf-8
""" This file implements the Site class. """

from matador.utils.cell_utils import cart2frac, frac2cart


class Site:
    def __init__(self, species, position, lattice_cart,
                 position_unit='fractional',
                 spin=None, voronoi_substructure=None,
                 chemical_shift=None, magnetic_shielding=None):
        self.lattice = lattice_cart
        self.set_position(position, position_unit)
        self.species = species
        self.spin = spin

        if voronoi_substructure is not None:
            assert self.species == voronoi_substructure[0]
            self.voronoi_substructure = voronoi_substructure[1]
        else:
            self.voronoi_substructure = None

        self.chemical_shift = chemical_shift
        self.magnetic_shielding = magnetic_shielding

    def __str__(self):
        return '{species:<3} {pos[0]:4.4f} {pos[1]:4.4f} {pos[2]:4.4f}'.format(species=self.species, pos=self.coords)

    def set_position(self, position, units):
        if '_coords' not in self.__dict__:
            self._coords = dict()
        if units == 'fractional':
            self._coords['fractional'] = position
            self._coords['cartesian'] = frac2cart(self.lattice, self.coords)
        elif units == 'cartesian':
            self._coords['cartesian'] = position
            self._coords['fractional'] = cart2frac(self.lattice, self.coords(units='cartesian'))
        else:
            raise RuntimeError('Unit system {} not understood, expecting `fractional`/`cartesian`'.format(units))

    @property
    def coords(self):
        return self._coords['fractional']

    def get_coords(self, units='fractional'):
        if units not in ['fractional', 'cartesian']:
            raise RuntimeError('Unit system {} not understood, expecting `fractional`/`cartesian`'.format(units))
        else:
            return self._coords[units]

    @property
    def coordination(self):
        if '_coordination' in self.__dict__:
            return self._coordination
        if self.voronoi_substructure is None:
            raise RuntimeError('Voronoi substructure not found.')

        coordination = {}
        eps = 0.05
        for atom, weight in self.voronoi_substructure:
            if weight >= 1 - eps:
                if atom not in coordination:
                    coordination[atom] = 1
                else:
                    coordination[atom] += 1
        self._coordination = coordination

        return coordination

    def displacement_between_sites(self, other_site):
        import numpy as np
        return np.asarray(self.get_coords(units='cartesian')) - np.asarray(other_site.get_coords(units='cartesian'))

    def distance_between_sites(self, other_site):
        import numpy as np
        return np.linalg.norm(self.displacement_between_sites(other_site))
