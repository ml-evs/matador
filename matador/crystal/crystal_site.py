# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements the Site class for handling
atomic sites.

"""

import numpy as np
import copy
from matador.utils.cell_utils import cart2frac, frac2cart, wrap_frac_coords


class Site:
    def __init__(self, species: str, position: list, lattice_cart,
                 position_unit='fractional', **site_data):

        self.lattice = copy.deepcopy(lattice_cart)
        self.set_position(position, position_unit)
        self.species = species
        self._occupancy = None
        self.site_data = {}

        if site_data.get('voronoi_substructure') is not None:
            assert self.species == site_data['voronoi_substructure'][0]
            self.site_data['voronoi_substructure'] = site_data['voronoi_substructure'][1]
            del site_data['voronoi_substructure']

        self.site_data.update(site_data)

    def __str__(self):
        site_str = '{species} {pos[0]:4.4f} {pos[1]:4.4f} {pos[2]:4.4f}'.format(species=self.species, pos=self.coords)
        for key in self.site_data:
            site_str += '\n{} = {}'.format(key, self.site_data[key])
        return site_str

    def __repr__(self):
        site_str = '{species} {pos[0]:4.4f} {pos[1]:4.4f} {pos[2]:4.4f}'.format(species=self.species, pos=self.coords)
        for key in self.site_data:
            site_str += '\n{} = {}'.format(key, self.site_data[key])
        return site_str

    def __deepcopy__(self, memo):
        from copy import deepcopy
        species, position, lattice = (deepcopy(x) for x in (self.species, self._coords['fractional'], self.lattice))
        site_data = deepcopy(self.site_data)
        return Site(species, position, lattice, position_unit='fractional', **site_data)

    def set_position(self, position, units):
        if len(position) != 3 or not all(isinstance(p, (float, int)) for p in position):
            raise RuntimeError('CrystalSite position has wrong shape: {}'.format(position))
        if '_coords' not in self.__dict__:
            self._coords = dict()
        if units == 'fractional':
            self._coords['fractional'] = wrap_frac_coords(
                [float(pos) for pos in position],
                remove=False
            )
            self._coords['cartesian'] = frac2cart(self.lattice, self.coords)
        elif units == 'cartesian':
            self._coords['cartesian'] = [float(pos) for pos in position]
            self._coords['fractional'] = wrap_frac_coords(
                cart2frac(self.lattice, self.coords(units='cartesian')),
                remove=False
            )
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
    def occupancy(self):
        if self._occupancy is None:
            self._occupancy = self.site_data.get("site_occupancy", 1.0)
        return self._occupancy

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
        return np.asarray(self.get_coords(units='cartesian')) - np.asarray(other_site.get_coords(units='cartesian'))

    def distance_between_sites(self, other_site):
        return np.linalg.norm(
            self.displacement_between_sites(other_site)
        )
