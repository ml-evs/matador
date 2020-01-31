# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements classes to store and manipulate electronic and
vibrational DOS, with or without projection data.

"""

import warnings
import numpy as np
import scipy.integrate
import scipy.interpolate
from matador.orm.orm import DataContainer
from matador.utils.chem_utils import KELVIN_TO_EV
from matador.utils.chem_utils import INVERSE_CM_TO_EV

EPS = 1e-6


class DensityOfStates(DataContainer):
    """ Generic class for density of states. """

    def __init__(self, data):
        """ Initialise the DOS and trim the DOS data arrays.

        Parameters:
            data (dict): dictionary containing the phonon dos data.

        """

        super().__init__(data)
        self._trim_dos()

    def _trim_dos(self):
        """ Trim the density of states/frequencies to only include the non-zero
        section of the DOS.

        """
        first_index = np.argmax(self._data['dos'] > EPS)
        last_index = len(self._data['dos']) - np.argmax(self._data['dos'][::-1] > EPS)
        self._trimmed_dos = self._data['dos'][first_index:last_index]
        self._trimmed_energies = self._data['energies'][first_index:last_index]

    @property
    def sample_dos(self):
        """ Return the calculated density of states, trimmed at each end to
        only include non-zero values.

        """
        return self._trimmed_dos

    @property
    def sample_energies(self):
        """ Return the energies corresponding to the trimmed DOS. """
        return self._trimmed_energies

    def plot_dos(self, **kwargs):
        """ Plot the density of states. """
        from matador.plotting.spectral_plotting import plot_spectral
        plot_spectral(
            self,
            phonons='Vibrational' in self.__class__.__name__,
            plot_dos=True,
            plot_bandstructure=False,
            **kwargs
        )


class VibrationalDOS(DensityOfStates):
    """ Specific class for phonon DOS data, including free energy integration. """

    @property
    def zero_point_energy(self):
        """ Computes the zero-point energy for the vibrational DOS,
        using

        .. math::

            E_{zp} = \\frac{1}{2}\\int F(\\omega)\\hbar\\omega\\,\\mathrm{d}\\omega.

        where .. math:: F(\\omega)
        is the vibrational density of states.

        """
        warnings.warn(
            'Imaginary frequency phonons found in this structure, ZP energy '
            'calculation will be unreliable',
            Warning
        )

        def integrand(omega):
            return self.vdos_function(omega) * INVERSE_CM_TO_EV * omega

        result = scipy.integrate.quad(
            integrand,
            self.sample_energies[0],
            self.sample_energies[-1]
        )

        return 0.5 * result[0]

    def vibrational_free_energy(self, temperatures=None):
        """ Computes the vibrational contribution to the free energy
        at a given set of temperatures, using

        .. math::

            F_{\\text{vib}}(T) = kT \\int F(\\omega) \\ln{\\left[1 - \\exp{-\\frac{\\hbar\\omega}{kT}}\\right]\\, \\mathrm{d}\\omega,

        Keyword arguments:
            temperature (list): list, array or float of temperatures.

        """
        if temperatures is None:
            temperatures = np.linspace(0, 600, num=5)

        temperatures = np.asarray(temperatures)
        free_energy = np.zeros_like(temperatures)
        errs = np.zeros_like(free_energy)

        min_energy = self.sample_energies[0]
        max_energy = self.sample_energies[-1]
        if min_energy < 0:
            warnings.warn(
                'Imaginary frequency phonons found in this structure, free energy '
                'calculation will be unreliable, using 0 K as lower limit of integration.',
                Warning
            )
            min_energy = 1e-5

        for ind, temperature in enumerate(temperatures):

            # if 0 K is requested, return 0 and move on
            if temperature == 0:
                free_energy[ind] = 0.0
                errs[ind] = 0.0
                continue

            kT = KELVIN_TO_EV * temperature

            def integrand(omega):
                return self.vdos_function(omega) * np.log(1 - np.exp(- INVERSE_CM_TO_EV * omega/kT))

            result = scipy.integrate.quad(
                integrand,
                min_energy,
                max_energy
            )

            free_energy[ind] = INVERSE_CM_TO_EV * kT * result[0]
            errs[ind] = result[1]

        if len(temperatures) == 1:
            return free_energy[0]

        return temperatures, free_energy

    @property
    def vdos_function(self):
        """ From the data arrays :attr:`sample_energies` and :attr:`sample_dos`,
        return an interpolated function to integrate.

        """
        return scipy.interpolate.interp1d(
            self.sample_energies,
            self.sample_dos,
            fill_value=(0, 0),
            bounds_error=False,
            copy=False
        )

    def plot_free_energy(self, temperatures=None, ax=None, **kwargs):
        """ Plot G(T) on the array of given temperatures. Default T is [0, 800].

        Keyword arguments:
            temperatures (list/np.ndarray): list or array of temperatures to plot.
                If the array/list has length 2, use these as the start and endpoints
                with 21 plotting points.
            ax (matplotlib.pyplot.Axis): axis object to plot onto.

        """
        from matador.plotting.temperature_plotting import plot_free_energy
        plot_free_energy(self, temperatures=temperatures, ax=ax, **kwargs)


class ElectronicDOS(DensityOfStates):
    """ Specific class for electronic DOS data. """
    pass
