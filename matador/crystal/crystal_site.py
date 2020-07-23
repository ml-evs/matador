# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements the Site class for handling
atomic sites.

"""

import numpy as np
from matador.utils.cell_utils import cart2frac, frac2cart, wrap_frac_coords
from matador.orm.orm import DataContainer


class Site(DataContainer):

    def __init__(self, species: str, position: list, lattice_cart,
                 position_unit='fractional', **site_data):

        if site_data.get('voronoi_substructure') is not None:
            assert self.species == site_data['voronoi_substructure'][0]
            site_data['voronoi_substructure'] = site_data['voronoi_substructure'][1]

        super().__init__(
            species=species,
            position=position,
            lattice_cart=lattice_cart,
            **site_data
        )

        self.set_position(position, position_unit)
        self._occupancy = None

        self.site_data = {}
        self.site_data.update(site_data)

    def __getitem__(self, key):
        """ Add extra look-up in `self.site_data` to
        :class:`DataContainer`'s `__getitem__`.

        Parameters:
            key (str): name of key or attribute to get.

        Raises:
            AttributeError: if key or attribute can't be found.

        """
        try:
            super().__getitem__(key)
        except (AttributeError, KeyError):
            pass

        try:
            return self.site_data[key]
        except KeyError:
            raise KeyError('Site has no data/site_data or implementation for requested key: "{}"'
                           .format(key))

    def __str__(self):
        site_str = '{species} {pos[0]:4.4f} {pos[1]:4.4f} {pos[2]:4.4f}'.format(species=self.species, pos=self.coords)
        for key in self.site_data:
            try:
                site_str += '\n{} = {:4.4f}'.format(key, float(self.site_data[key]))
            except ValueError:
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

    @property
    def species(self):
        return self._data['species']

    @species.setter
    def species(self, value):
        self._data['species'] = value

    @property
    def lattice(self):
        return self._data['lattice_cart']

    def get_coords(self, units='fractional'):
        if units not in ['fractional', 'cartesian']:
            raise RuntimeError('Unit system {} not understood, expecting `fractional`/`cartesian`'.format(units))
        else:
            return self._coords[units]

    @property
    def occupancy(self):
        if "site_occupancy" not in self._data:
            self._data["site_occupancy"] = 1.0
        return self._data["site_occupancy"]

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
