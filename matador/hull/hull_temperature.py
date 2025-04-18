# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the `TemperatureDependentHull` class
for assessing phase stability from finite temperature free energies.

"""

import copy
import numpy as np

from matador.hull.hull_ensemble import EnsembleHull
from matador.orm.spectral import VibrationalDOS


class TemperatureDependentHull(EnsembleHull):
    """Leverages `EnsembleHull` to construct temperature dependent
    hulls from phonon calculations.

    """

    data_key = "temperature"
    energy_key = "free_energy_per_atom"

    def __init__(
        self, cursor, energy_key="enthalpy_per_atom", temperatures=None, use_castep_thermo=False, **kwargs
    ):

        self.temperatures = temperatures
        self.use_castep_thermo = use_castep_thermo
        print("use_castep_thermo", self.use_castep_thermo)
        if temperatures is None:
            self.temperatures = np.linspace(0, 800, 21)

        _cursor = copy.deepcopy(cursor)

        # prepare the cursor by computing free energies
        # and store it in the format expected by EnsembleHull
        for ind, doc in enumerate(cursor):
            if not isinstance(doc, VibrationalDOS):
                #Computing the vibrational free energy CAN be VERY SLOW
                #option to use the castep thermo data if it is present and the temperatures are within the range of the castep thermo data.
                vib_free_energies = None
                if self.use_castep_thermo and "thermo_temps" in doc:
                    castep_temps = np.array(doc["thermo_temps"])
                    castep_vib_free_energies = np.array([x for x in doc["thermo_free_energy"].values()])
                    #ONLY DO THIS if you are actually interpolating.. use 0.1 K tolerance.
                    if min(self.temperatures) > min(castep_temps)-0.1 and max(self.temperatures) < max(castep_temps) + 0.1:
                        vib_free_energies = [np.interp(T, castep_temps, castep_vib_free_energies) / doc["num_atoms"] for T in self.temperatures]
                if vib_free_energies is None:
                    _doc = VibrationalDOS(doc)
                    temps, vib_free_energies = _doc.vibrational_free_energy(
                        temperatures=self.temperatures
                    )

                _cursor[ind][self.data_key] = {}
                _cursor[ind][self.data_key][self.energy_key] = (
                    np.ones_like(self.temperatures) * _cursor[ind][energy_key]
                )
                _cursor[ind][self.data_key][self.energy_key] += vib_free_energies
                _cursor[ind][self.data_key]["temperatures"] = self.temperatures


        super().__init__(
            cursor=_cursor,
            data_key=self.data_key,
            energy_key=self.energy_key,
            chempot_energy_key=energy_key,
            parameter_key="temperatures",
            update_chem_pots=True,
            **kwargs
        )

    def plot_hull(self, **kwargs):
        """Hull plot helper function."""
        from matador.plotting.hull_plotting import plot_temperature_hull

        return plot_temperature_hull(self, **kwargs)
