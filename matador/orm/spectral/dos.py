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

    def _from_dispersion(self, data, **kwargs):
        """ Convert a Dispersion instance to a DOS. """
        _data = {}
        dos, energies = self.bands_as_dos(data, **kwargs)

        for key in [
                'eigs_s_k', 'eigs_q',
                'kpoint_path', 'kpoint_weights',
                'num_kpoints', 'num_qpoints', 'num_modes']:
            if key in data:
                _data[key] = data[key]

        _data['dos'] = dos
        _data['energies'] = energies
        return _data

    @property
    def sample_dos(self):
        """ Return the calculated density of states, trimmed at each end to
        only include non-zero values.

        """
        return self._trimmed_dos

    @property
    def sample_energies(self):
        """ Return the energies corresponding to the trimmed DOS. """
        return self._trimmed_energies * INVERSE_CM_TO_EV

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

    @staticmethod
    def bands_as_dos(bands, gaussian_width=None):
        """ Convert bands data to DOS data. """
        if 'eigs_s_k' in bands:
            eigs_key = 'eigs_s_k'
            if gaussian_width is None:
                gaussian_width = 0.1
        elif 'eigs_q' in bands:
            eigs_key = 'eigs_q'
            if gaussian_width is None:
                gaussian_width = 10
        else:
            raise RuntimeError('Missing eigenvalue keys from bands data.')

        raw_eigs = np.asarray(bands[eigs_key]) - bands.get('fermi_energy', 0)
        raw_weights = np.ones_like(raw_eigs)
        if 'kpoint_weights' in bands:
            for sind, _ in enumerate(bands[eigs_key]):
                for kind, _ in enumerate(bands[eigs_key][sind]):
                    raw_weights[sind, kind, :] = bands['kpoint_weights'][kind]

        if len(raw_weights) != 1:
            if len(raw_weights) > 2:
                raise NotImplementedError('Non-collinear spin not supported')
            spin_dos = dict()
            keys = ['up', 'down']
            for sind, _ in enumerate(raw_weights):
                spin_dos[keys[sind]], energies = DensityOfStates._cheap_broaden(
                    bands[eigs_key][sind].flatten(),
                    weights=raw_weights[sind].flatten(),
                    gaussian_width=gaussian_width
                )

            return spin_dos, energies

        dos, energies = DensityOfStates._cheap_broaden(
            raw_eigs.flatten(), weights=raw_weights.flatten(), gaussian_width=gaussian_width
        )

        return dos, energies

    @staticmethod
    def _cheap_broaden(eigs, weights=None, gaussian_width=None):
        """ Quickly broaden and bin a set of eigenvalues.

        Parameters:
            eigs (numpy.ndarray): eigenvalue array.
            weights (numpy.ndarray): array of weights.

        Keyword arguments:
            gaussian_width (float): width of gaussian broadening
                to apply.

        Returns:
            Two arrays containing the DOS and energies.

        """
        if gaussian_width is None:
            gaussian_width = 0.1

        hist, energies = np.histogram(eigs, weights=weights, bins=1001)
        if gaussian_width == 0:
            return hist, energies

        # shift bin edges to bin centres
        energies -= energies[1] - energies[0]
        energies = energies[:-1]
        new_energies = np.reshape(energies, (1, len(energies)))
        new_energies = new_energies - np.reshape(energies, (1, len(energies))).T
        dos = np.sum(hist * np.exp(-(new_energies)**2 / gaussian_width), axis=1)
        dos = np.divide(dos, np.sqrt(2 * np.pi * gaussian_width**2))

        return dos, energies


class VibrationalDOS(DensityOfStates):
    """ Specific class for phonon DOS data, including free energy integration. """

    @property
    def zero_point_energy_from_dos(self):
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

    @property
    def zero_point_energy(self):
        """ Computes and return the zero-point energy from frequency data. """

        if 'eigs_q' not in self._data:
            raise RuntimeError('Unable to compute free energies without frequency data.')

        eigs = self.eigs[0] * INVERSE_CM_TO_EV
        min_energy = np.min(eigs)
        freq_cutoff = 1e-3
        if min_energy < freq_cutoff:
            warnings.warn(
                'Imaginary frequency phonons found in this structure, free energy '
                'calculation will be unreliable, using 0 K as lower limit of integration.'
            )

        zpe = 0.0
        for mode_ind in range(self.num_modes):
            for qpt_ind in range(self.num_qpoints):
                contrib = eigs[mode_ind][qpt_ind]
                if 'kpoint_weights' in self._data:
                    contrib *= self.kpoint_weights[qpt_ind]
                else:
                    contrib /= self._data['num_qpoints']
                zpe += contrib

        return 0.5 * zpe

    def vibrational_free_energy(self, temperatures=None):
        """ Computes and returns the vibrational contribution to the free
        energy, including zero-point energy, from the phonon frequencies.

        Parameters:
            temperatures (list): list or array of temperatures to compute
                G(T) at.

        Returns:
            (np.ndarray, np.ndarray): temperature and energy array.

        """

        if temperatures is None:
            temperatures = np.linspace(0, 600, num=5)

        try:
            _ = len(temperatures)
        except AttributeError:
            temperatures = [temperatures]

        if 'eigs_q' not in self._data:
            raise RuntimeError('Unable to compute free energies without frequency data.')

        eigs = self._data['eigs_q'][0] * INVERSE_CM_TO_EV

        temperatures = np.asarray(temperatures)
        free_energy = np.zeros_like(temperatures, dtype=np.float64)

        min_energy = np.min(eigs)
        freq_cutoff = 1e-3
        if min_energy < freq_cutoff:
            warnings.warn(
                'Imaginary frequency phonons found in this structure, free energy '
                'calculation will be unreliable, using 0 K as lower limit of integration.'
            )

        free_energy += self.zero_point_energy

        for ind, temperature in enumerate(temperatures):
            if temperature < 1e-6:
                continue

            kT = KELVIN_TO_EV * temperature

            for mode_ind in range(self._data['num_modes']):
                for qpt_ind in range(self._data['num_qpoints']):
                    freq = eigs[mode_ind][qpt_ind]
                    if freq > freq_cutoff and freq / kT < 32:
                        contrib = kT * np.log(1 - np.exp(-freq/kT))
                        if 'kpoint_weights' in self._data:
                            contrib *= self._data['kpoint_weights'][qpt_ind]
                        else:
                            contrib /= self._data['num_qpoints']

                        free_energy[ind] += contrib

        if len(temperatures) == 1:
            return free_energy[0]

        return temperatures, free_energy

    def vibrational_free_energy_from_dos(self, temperatures=None):
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
            min_energy = 1e-3

        for ind, temperature in enumerate(temperatures):

            # if 0 K is requested, return 0 and move on
            if temperature == 0:
                free_energy[ind] = 0.0
                errs[ind] = 0.0
                continue

            kT = KELVIN_TO_EV * temperature

            def integrand(omega):
                return self.vdos_function(omega) * np.log(1 - np.exp(-omega/kT))

            result = scipy.integrate.quad(
                integrand,
                min_energy,
                max_energy
            )

            free_energy[ind] = kT * result[0]
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
