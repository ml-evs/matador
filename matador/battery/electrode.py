# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the Electrode class used for calculating
relevant electrode properties from phase diagrams.

"""


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
    _valence_data = {'Li': 1, 'Na': 1, 'K': 1,
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
        import numpy as np
        trim = None
        if np.isnan(capacities[-1]):
            trim = -1
        return np.sum(voltages[1:trim] * np.diff(capacities[:trim])) / np.max(capacities[:trim])
