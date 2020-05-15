# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the Electrode class used for calculating
relevant electrode properties from phase diagrams.

"""

import numpy as np
from typing import List, Tuple

from matador.utils.chem_utils import get_formula_from_stoich

EPS = 1e-12


class Electrode:
    """ The Electrode class wraps phase diagrams of battery electrode
    materials, and provides useful methods for electrochemical analysis.

    Note:
        The units for relevant quantities are as follows:

            gravimetric capacity --> mAh/g
            volumetric capacity  --> mAh/cm^3
            gravimetric energy density --> Wh/kg
            volumetric energy density --> Wh/l

    Attributes:
        kind (str): either 'anode' or 'cathode'. Effects conventions
            used when computing capacities.
        ion (str): species label for the active ion.
        valence (int): number of electrons transferred per ion inserted/
            removed. Currently assumes perfect Coulombic efficiency, and
            full oxidation at each stage (e.g. Li->Li+ and Mg->Mg2+).

    """

    _kinds = ['anode', 'cathode']
    valence_data = {'Li': 1, 'Na': 1, 'K': 1,
                    'Mg': 2, 'Ca': 2}

    def __init__(self, active_ion, base_material, phases, kind='anode', valence=None):
        """ Initialise Electrode material for the given ion from a convex hull.

        Parameters:
            active_ion (str): label of species to use as active ion.
            base_formula (str): formula of the starting material (either
                "empty" anode or "full" cathode).
            phases (:obj:`QueryConvexHull`/:obj:`list` of :obj:`dict`):
                either a QueryConvexHull object for the material system
                of interest, or a list of phases.

        Keyword arguments:
            kind (str): either 'anode' or 'cathode'. This decides what
                convention to use when computing capacities.
            valence (int): the valence of the active ion, if not included
                this will be looked up for common ions.

        """

        if kind not in self._kinds:
            raise TypeError('Electrode kind must be one of {}'.format(self._kinds))
        self.kind = kind
        self.ion = active_ion

        if valence is not None:
            self.valence = valence
        else:
            self.valence = self._valence_data.get('active_ion')
        if self.valence is None:
            raise RuntimeError('Unable to lookup valence of {}, please pass it manually'
                               .format(self.ion))

        # TODO: detect different materials inside convex hull
        # or directly compute if just a list of phases

    @property
    def gravimetric_capacities(self):
        pass

    @property
    def volumetric_capacities(self):
        pass

    @property
    def voltages(self):
        pass

    @property
    def average_voltage(self):
        pass

    @property
    def gravimetric_energy_density(self):
        pass

    @property
    def volumetric_energy_density(self):
        pass

    @property
    def max_gravimetric_capacity(self):
        pass

    @property
    def max_volumetric_capacity(self):
        pass

    @property
    def voltage_curve(self):
        pass

    @property
    def pdf_vs_voltage(self):
        pass

    @property
    def pdf_vs_capacity(self):
        pass

    @property
    def pxrd_vs_voltage(self):
        pass

    @property
    def pxrd_vs_capacity(self):
        pass

    @classmethod
    def calculate_average_voltage(cls, capacities, voltages):
        """ For a given set of capacities and voltages, compute
        the average voltage during charge/discharge.

        """
        trim = None
        if np.isnan(capacities[-1]):
            trim = -1
        return np.sum(voltages[1:trim] * np.diff(capacities[:trim])) / np.max(capacities[:trim])

    @classmethod
    def _find_starting_materials(self, points, stoichs):
        """ Iterate over an array of compositions and energies to find
        the starting points of the electrode, i.e. those with zero concentration
        of the active ion.

        """

        endpoints = []
        endstoichs = []
        for ind, point in enumerate(points):
            if point[0] == 0 and point[1] != 0 and point[1] != 1:
                if not any([point.tolist() == test_point.tolist() for test_point in endpoints]):
                    endpoints.append(point)
                    endstoichs.append(stoichs[ind])

        return endpoints, endstoichs


class VoltageProfile:
    """ Container class for data associated with a voltage profile.

    Attributes:
        starting_stoichiometry: the initial stoichiometry of the electrode.
        capacities: list of gravimetric capacities in mAh/g.
        voltage: list of voltages at each capacity step.
        average_voltage: the average voltage across the full cycle.
        reactions: a list of (coefficient, formula) tuples showing the
            progression of the balanced reaction, e.g.
            `[((1, "PSn")), ((1, "LiP"), (1, "LiSn"))]`
        active_ion: species label of the active ion.

    """
    def __init__(
        self,
        starting_stoichiometry: Tuple[Tuple[str, float], ...],
        capacities: List[float],
        voltages: List[float],
        average_voltage: float,
        active_ion: str,
        reactions: List[Tuple[Tuple[float, str], ...]],
    ):
        """ Initialise the voltage profile with the given data. """

        n_steps = len(voltages)
        if any(_len != n_steps for _len in [len(capacities), len(reactions)+1]):
            raise RuntimeError(
                "Invalid size of initial arrays, capacities and voltages must have same length."
                "reactions array must have length 1 smaller than voltages"
            )

        self.starting_stoichiometry = starting_stoichiometry
        self.starting_formula = get_formula_from_stoich(starting_stoichiometry)
        self.capacities = capacities
        self.average_voltage = average_voltage
        self.voltages = voltages
        self.reactions = reactions
        self.active_ion = active_ion

    def voltage_summary(self, csv=False):
        """ Prints a voltage data summary.

        Keyword arguments:
            csv (bool/str): whether or not to write a CSV file containing the data.
                If this contains a string use this as the filename.

        """

        data_str = '# {} into {}\n'.format(self.active_ion, self.starting_formula)
        data_str += "# Average voltage: {:4.2f} V\n".format(self.average_voltage)
        data_str += '# {:^10} \t{:^10}\n'.format('Q (mAh/g)', 'Voltage (V)')

        for idx, _ in enumerate(self.voltages):
            data_str += '{:>10.2f} \t{:>10.8f}'.format(self.capacities[idx], self.voltages[idx])
            if idx != len(self.voltages) - 1:
                data_str += '\n'

        if csv:
            if isinstance(csv, str):
                fname = csv
            else:
                fname = '{}_voltages.csv'.format(self.starting_formula)
            with open(fname, 'w') as f:
                f.write(data_str)

        return 'Voltage data:\n\n{}'.format(data_str)

    def __repr__(self):
        self.voltage_summary(csv=False)
