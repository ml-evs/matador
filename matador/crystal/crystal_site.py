# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements the Site class for handling
atomic sites.

"""

import numpy as np
from matador.utils.cell_utils import cart2frac, frac2cart, wrap_frac_coords
from matador.orm.orm import DataContainer


class Site(DataContainer):
    """The Site class contains a description of an individual
    site within a 3D periodic Crystal.

    """

    # This dictionary defines the map between fields in :obj:`Crystal`
    # that correspond to arrays of site properties and between the
    # relevant keys the :obj:`Site` object
    _crystal_key_map = {
        "site_occupancy": "site_occupancy",
        "chemical_shielding_isos": "chemical_shielding_iso",
        "chemical_shift_isos": "chemical_shift_iso",
        "magnetic_shielding_tensors": "magnetic_shielding_tensor",
        "electric_field_gradients": "electric_field_gradient",
        "chemical_shift_anisos": "chemical_shift_aniso",
        "chemical_shift_asymmetries": "chemical_shift_asymmetry",
        "quadrupolar_couplings": "quadrupolar_coupling",
        "quadrupolar_asymmetries": "quadrupolar_asymmetry",
        "voronoi_substructure": "voronoi_substructure",
    }

    def __init__(
        self,
        species: str,
        position: list,
        lattice,
        position_unit="fractional",
        mutable: bool = False,
        **site_data,
    ):
        """Initialise a Site object from its species, position and
        a reference to the lattice it exists in. Any other keys will be made available
        as site-level values.

        """
        if site_data.get("voronoi_substructure") is not None:
            assert self.species == site_data["voronoi_substructure"][0]
            site_data["voronoi_substructure"] = site_data["voronoi_substructure"][1]

        # DataContainer will take a copy of all data passed to it, but lets keep
        # lattice as a reference so that it can change externally
        self._lattice = lattice

        super().__init__(
            species=species, position=position, site_data=site_data, mutable=mutable
        )

        self._data["lattice_cart"] = self._lattice

        self.set_position(position, position_unit)
        self._occupancy = None

        self.site_data = {}
        self.site_data.update(site_data)

    def get(self, key, default=None):
        try:
            return self[key]
        except (KeyError, AttributeError):
            return default

    def __getitem__(self, key):
        """Add extra look-up in `self.site_data` to
        :class:`DataContainer`'s `__getitem__`.

        Parameters:
            key (str): name of key or attribute to get.

        Raises:
            AttributeError: if key or attribute can't be found.

        """
        if isinstance(key, int):
            raise ValueError("Object does not support indexing")
        try:
            super().__getitem__(key)
        except (AttributeError, KeyError):
            pass

        try:
            return self.site_data[key]
        except KeyError:
            raise KeyError(
                'Site has no data/site_data or implementation for requested key: "{}"'.format(
                    key
                )
            )

    def __setitem__(self, key: str, item):
        if key not in self.site_data or self.site_data[key] is None:
            self.site_data[key] = item
            return

        elif self.site_data[key] != item:
            try:
                import math

                if math.isnan(item) and math.isnan(self.site_data[key]):
                    return
            except TypeError:
                pass
            raise AttributeError(
                "Cannot assign value {} to existing key {} with value {}".format(
                    item, key, self.site_data[key]
                )
            )

    def __str__(self):
        site_str = "{species} {pos[0]:4.4f} {pos[1]:4.4f} {pos[2]:4.4f}".format(
            species=self.species, pos=self.coords
        )
        for key in self.site_data:
            try:
                site_str += "\n{} = {:4.4f}".format(
                    key, np.asarray(self.site_data[key])
                )
            except (ValueError, TypeError):
                with np.printoptions(precision=2, threshold=6, edgeitems=2):
                    site_str += "\n{} = \n{}".format(
                        key,
                        "\n".join(
                            f"  {row}"
                            for row in np.asarray(self.site_data[key])
                            .__str__()
                            .split("\n")
                        ),
                    )

        site_str += "\n---"

        return site_str

    def __deepcopy__(self, memo):
        from copy import deepcopy

        species, position, lattice = (
            deepcopy(x) for x in (self.species, self.coords, self.lattice)
        )
        site_data = deepcopy(self.site_data)
        return Site(species, position, lattice, position_unit="fractional", **site_data)

    def set_position(self, position, units):
        if len(position) != 3 or not all(isinstance(p, (float, int)) for p in position):
            raise RuntimeError(
                "CrystalSite position has wrong shape: {}".format(position)
            )
        if not hasattr(self, "_coords"):
            self._coords = dict()
        if units == "fractional":
            self._coords["fractional"] = wrap_frac_coords(
                [float(pos) for pos in position], remove=False
            )
        elif units == "cartesian":
            self._coords["fractional"] = wrap_frac_coords(
                cart2frac(self.lattice, self.coords), remove=False
            )
        else:
            raise RuntimeError(
                "Unit system {} not understood, expecting `fractional`/`cartesian`".format(
                    units
                )
            )

    @property
    def coords(self):
        return np.asarray(self._coords["fractional"])

    @property
    def coords_cartesian(self):
        return np.asarray(frac2cart(self.lattice, self.coords))

    @property
    def species(self):
        return self._data["species"]

    @species.setter
    def species(self, value):
        self._data["species"] = value

    @property
    def lattice(self):
        try:
            return self._lattice.lattice_cart
        except AttributeError:
            return self._lattice

    @property
    def occupancy(self):
        if "site_occupancy" not in self._data:
            self._data["site_occupancy"] = 1.0
        return self._data["site_occupancy"]

    @property
    def coordination(self):
        if "_coordination" in self.__dict__:
            return self._coordination
        if self.voronoi_substructure is None:
            raise RuntimeError("Voronoi substructure not found.")

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
        return self.coords_cartesian - other_site.coords_cartesian

    def distance_between_sites(self, other_site):
        return np.linalg.norm(self.displacement_between_sites(other_site))
