# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements classes to store and manipulate electronic and
vibrational DOS, with or without projection data.

"""

import warnings
import copy
import functools
import numpy as np
import scipy.integrate
import scipy.interpolate

from matador.orm.orm import DataContainer
from matador.utils.chem_utils import KELVIN_TO_EV, INVERSE_CM_TO_EV

from .dispersion import Dispersion

EPS = 1e-6
MIN_PHONON_FREQ = -0.01
FREQ_CUTOFF = 1e-12


class DensityOfStates(Dispersion, DataContainer):
    """Generic class for density of states."""

    required_keys = ["dos", "energies"]

    def __init__(self, *args, **kwargs):
        """Initialise the DOS and trim the DOS data arrays.

        Parameters:
            data (dict/Dispersion): dictionary containing the phonon dos data, or
                a dispersion object to convert.

        """

        if kwargs.get("gaussian_width") is not None:
            self.gaussian_width = kwargs["gaussian_width"]

        if args and isinstance(args[0], dict):
            data = args[0]
        else:
            data = kwargs
        # as we can also construct a DOS from arbitarary kpoint/energy data,
        # check that we've been passed this first
        if isinstance(data, Dispersion) or (
            isinstance(data, dict)
            and not any(key in data for key in ["spin_dos", "dos", "pdos"])
        ):
            data = self._from_dispersion(data)
        elif isinstance(data, DensityOfStates):
            data = copy.deepcopy(DensityOfStates._data)

        # Attempted workaround for old OPTADOS bug where spin-polarized PDOS would
        # ignore set_efermi_zero. Requires any other DOS data to have been computed
        # on same grid with OptaDOS, and assumes Fermi level is the same for each spin.
        # https://github.com/optados-developers/optados/issues/24
        if (
            "pdos" in data
            and len(set(proj[2] for proj in data["pdos"]["projectors"])) == 2
        ):
            # check for energies either side of 0, if none are found, try to shift to other fermi values
            if len(data["pdos"]["energies"]) == len(data.get("energies", [])) and (
                np.max(data["pdos"]["energies"]) < 0
                or np.min(data["pdos"]["energies"]) > 0
            ):
                correction = np.max(data["pdos"]["energies"]) - np.max(data["energies"])
                data["pdos"]["energies"] -= correction
                warnings.warn(
                    f"""Corrected PDOS energies with difference between total DOS and PDOS grid.
Guessed Fermi level = {correction:4.3f} eV, assumed constant for both spins.
In the future, please use an updated version of OptaDOS.
This correction is to account for the OptaDOS bug #24:
https://github.com/optados-developers/optados/issues/24
"""
                )

        # trigger generation of dos key from spin dos
        if "dos" not in data and "spin_dos" in data:
            data["dos"] = np.asarray(data["spin_dos"]["up"]) + np.asarray(
                data["spin_dos"]["down"]
            )

        if "dos" not in data and "pdos" in data:
            data["dos"] = np.sum(
                data["pdos"]["pdos"][proj] for proj in data["pdos"]["projectors"]
            )
            data["energies"] = data["pdos"]["energies"]
            warnings.warn(
                "Total DOS created from sum of projected DOS, which may not be at all reliable."
            )

        if "spin_dos" not in data and "pdos" in data:
            spin_channels = set(proj[2] for proj in data["pdos"]["projectors"])
            if len(spin_channels) == 2:
                data["spin_dos"] = {}

                # mask negative contributions with 0
                for proj in data["pdos"]["projectors"]:
                    data["pdos"]["pdos"][proj] = np.ma.masked_where(
                        data["pdos"]["pdos"][proj] < 0,
                        data["pdos"]["pdos"][proj],
                        copy=True,
                    )
                    np.ma.set_fill_value(data["pdos"]["pdos"][proj], 0)
                    data["pdos"]["pdos"][proj] = np.ma.filled(
                        data["pdos"]["pdos"][proj]
                    )

                for channel in spin_channels:
                    data["spin_dos"][channel] = np.sum(
                        data["pdos"]["pdos"][proj]
                        for proj in data["pdos"]["projectors"]
                        if proj[2] == channel
                    )
                data["energies"] = data["pdos"]["energies"]
                warnings.warn(
                    "Total spin DOS created from sum of projected DOS, which may not be at all reliable."
                )

        super().__init__(data)
        self._trim_dos()

    def _trim_dos(self):
        """Trim the density of states/frequencies to only include the non-zero
        section of the DOS.

        """
        dos = self._data["dos"]
        first_index = np.argmax(dos > EPS)
        last_index = len(dos) - np.argmax(dos[::-1] > EPS)
        self._trimmed_dos = dos[first_index:last_index]
        self._trimmed_energies = self._data["energies"][first_index:last_index]

    def _from_dispersion(self, data, **kwargs):
        """Convert a Dispersion instance to a DOS."""
        _data = {}
        dos, energies = self.bands_as_dos(data, gaussian_width=self.gaussian_width)

        for key in data:
            _data[key] = data[key]

        _data["dos"] = dos
        _data["energies"] = energies

        if "Vibrational" not in self.__class__.__name__:
            warnings.warn("Loaded DOS from .bands file with naive Gaussian smearing.")

        return _data

    @property
    def sample_dos(self):
        """Return the calculated density of states, trimmed at each end to
        only include non-zero values.

        """
        return self._trimmed_dos

    @property
    def sample_energies(self):
        """Return the energies corresponding to the trimmed DOS."""
        return self._trimmed_energies

    def plot_dos(self, **kwargs):
        """Plot the density of states."""
        from matador.plotting.spectral_plotting import plot_spectral

        _kwargs = {
            "plot_dos": True,
            "plot_pdos": "pdos" in self,
            "plot_bandstructure": False,
            "phonons": "Vibrational" in self.__class__.__name__,
        }
        _kwargs.update(kwargs)
        plot_spectral(self, **_kwargs)

    @staticmethod
    def bands_as_dos(bands, gaussian_width=0.1):
        """Convert bands data to DOS data."""
        if "eigs_s_k" in bands:
            eigs_key = "eigs_s_k"
        elif "eigs_q" in bands:
            eigs_key = "eigs_q"
        else:
            raise RuntimeError("Missing eigenvalue keys from bands data.")

        raw_eigs = np.asarray(bands[eigs_key]) - bands.get("fermi_energy", 0)
        raw_weights = np.ones_like(raw_eigs)
        if "kpoint_weights" in bands:
            for sind, _ in enumerate(bands[eigs_key]):
                for kind, _ in enumerate(bands[eigs_key][sind][0]):
                    raw_weights[sind, :, kind] = bands["kpoint_weights"][kind]

        if len(raw_weights) != 1:
            if len(raw_weights) > 2:
                raise NotImplementedError("Non-collinear spin not supported")
            spin_dos = dict()
            keys = ["up", "down"]
            for sind, _ in enumerate(raw_weights):
                spin_dos[keys[sind]], energies = DensityOfStates._cheap_broaden(
                    bands[eigs_key][sind].flatten(),
                    weights=raw_weights[sind].flatten(),
                    gaussian_width=gaussian_width,
                )

                if "spin_fermi_energy" in bands:
                    energies -= bands["spin_fermi_energy"][0]

            return spin_dos, energies

        dos, energies = DensityOfStates._cheap_broaden(
            raw_eigs.flatten(),
            weights=raw_weights.flatten(),
            gaussian_width=gaussian_width,
        )

        if "fermi_energy" in bands:
            energies -= bands["fermi_energy"]

        return dos, energies

    @staticmethod
    def _cheap_broaden(eigs, weights=None, gaussian_width=None):
        """Quickly broaden and bin a set of eigenvalues.

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
        dos = np.sum(hist * np.exp(-((new_energies) ** 2) / gaussian_width), axis=1)
        dos = np.divide(dos, np.sqrt(2 * np.pi * gaussian_width**2))

        return dos, energies


class VibrationalDOS(DensityOfStates):
    """Specific class for phonon DOS data, including free energy integration."""

    gaussian_width = 10 * INVERSE_CM_TO_EV

    @property
    def debye_temperature(self):
        """Returns the Debye temperature in K."""
        return self.debye_freq / KELVIN_TO_EV

    @property
    def debye_freq(self):
        """Returns the Debye frequency in eV."""
        return np.max(self.eigs)

    @property
    def zpe(self):
        """The zero-point energy per atom as computed from frequency data."""
        if "zero_point_energy_per_atom" not in self._data:
            if "eigs_q" not in self._data:
                raise RuntimeError("Unable to compute ZPE without frequency data.")

            zpe = self._compute_zero_point_energy(
                self.eigs, self.num_kpoints, kpoint_weights=self.kpoint_weights
            )
            self["zero_point_energy"] = zpe
            self["zero_point_energy_per_atom"] = zpe / (self.num_modes / 3)

        return self._data["zero_point_energy_per_atom"]

    @staticmethod
    def _compute_zero_point_energy(eigs, num_kpoints, kpoint_weights=None):
        """Computes and returns the zero-point energy of the cell
        in eV from frequency data.

        Parameters:
            eigs (np.ndarray): phonon eigenvalues (in eV) array
                in any shape. If `kpoint_weights` is passed, then
                at least one axis of eigs must match the length
                of `kpoint_weights`.
            num_kpoints (int): the number of kpoints at which
                these eigenvalues were calculated.

        Keyword arguments:
            kpoint_weights (np.ndarray): array of weights to use
                for each kpoint.

        """
        min_energy = np.min(eigs)
        if min_energy < MIN_PHONON_FREQ:
            warnings.warn(
                "Imaginary frequency phonons found in this structure {:.1f}, ZPE "
                "calculation will be unreliable, using 0 eV as lower limit of integration.".format(
                    min_energy
                )
            )

        if kpoint_weights is not None:
            # if kpoint weights are all the same, discard them
            # and just normalise by number of kpoints
            if len(np.unique(kpoint_weights)) == 1:
                kpoint_weights = None
            else:
                eigs_shape = np.shape(eigs)
                if not any(len(kpoint_weights) == axis for axis in eigs_shape):
                    raise RuntimeError(
                        "Unable to match eigs with shape {} with kpoint weights of length {}".format(
                            eigs_shape, len(kpoint_weights)
                        )
                    )
        _eigs = np.copy(eigs)
        if kpoint_weights is not None:
            _eigs = _eigs * kpoint_weights
        else:
            _eigs /= num_kpoints

        return 0.5 * np.sum(np.ma.masked_where(_eigs < 0.0, _eigs, copy=False))

    def vibrational_free_energy(self, temperatures=None):
        """Computes and returns the vibrational contribution to the free
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
        except TypeError:
            temperatures = [temperatures]

        if "eigs_q" not in self._data:
            raise RuntimeError(
                "Unable to compute free energies without frequency data."
            )

        temperatures = np.asarray(temperatures)
        free_energy = np.zeros_like(temperatures, dtype=np.float64)

        min_energy = np.min(self._data["eigs_q"][0])
        if min_energy < MIN_PHONON_FREQ:
            warnings.warn(
                "Imaginary frequency phonons found in this structure {:.1f}, free energy "
                "calculation will be unreliable, using {} eV as lower limit of integration.".format(
                    min_energy, FREQ_CUTOFF
                )
            )

        for ind, temperature in enumerate(temperatures):
            free_energy[ind] = self.compute_free_energy(temperature)

        if len(temperatures) == 1:
            return free_energy[0]

        return temperatures, free_energy

    @functools.lru_cache(100)
    def compute_free_energy(self, temperature):
        """Compute the vibrational free energy at the given temperature, using
        lru_cache to avoid doing much extra work. Uses minimum temperature cutoff
        of 1e-9, below which it returns just the ZPE (unless T < 0 K).

        Raises:
            RuntimeError: if temperature is < 0 K.

        Returns:
            float: vibrational free energy per atom, including ZP correction.

        """

        free_energy = 0.0
        kT = KELVIN_TO_EV * temperature
        if temperature < 0.0:
            raise RuntimeError(
                "Not calculating free energies at T = {} K < 0 K".format(temperature)
            )

        if temperature < 1e-9:
            return self.zpe

        for mode_ind in range(self.num_modes):
            for qpt_ind in range(self.num_qpoints):
                freq = self._data["eigs_q"][0][mode_ind][qpt_ind]
                if freq > FREQ_CUTOFF and freq / kT < 32:
                    contrib = kT * np.log(1 - np.exp(-freq / kT))
                    if "kpoint_weights" in self._data:
                        contrib *= self.kpoint_weights[qpt_ind]
                    else:
                        contrib /= self.num_qpoints

                    free_energy += contrib

        # normalize by number of atoms
        free_energy /= self.num_modes / 3

        # add on zpe per atom
        free_energy += self.zpe

        return free_energy

    def vibrational_free_energy_from_dos(self, temperatures=None):
        """Computes the vibrational contribution to the free energy
        at a given set of temperatures.

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
        if min_energy < 0.01:
            warnings.warn(
                "Imaginary frequency phonons found in this structure {:.1f}, free energy "
                "calculation will be unreliable, using {} eV as lower limit of integration.".format(
                    min_energy, FREQ_CUTOFF
                ),
                Warning,
            )
            min_energy = FREQ_CUTOFF

        for ind, temperature in enumerate(temperatures):

            # if 0 K is requested, return 0 and move on
            if temperature == 0:
                free_energy[ind] = 0.0
                errs[ind] = 0.0
                continue

            kT = KELVIN_TO_EV * temperature

            def integrand(omega):
                return self.vdos_function(omega) * np.log(1 - np.exp(-omega / kT))

            result = scipy.integrate.quad(integrand, min_energy, max_energy)

            free_energy[ind] = kT * result[0]
            errs[ind] = result[1]

        if len(temperatures) == 1:
            return free_energy[0]

        return temperatures, free_energy

    @property
    def vdos_function(self):
        """From the data arrays :attr:`sample_energies` and :attr:`sample_dos`,
        return an interpolated function to integrate.

        """
        return scipy.interpolate.interp1d(
            self.sample_energies,
            self.sample_dos,
            fill_value=(0, 0),
            bounds_error=False,
            copy=False,
        )

    def plot_free_energy(self, temperatures=None, ax=None, **kwargs):
        """Plot G(T) on the array of given temperatures. Default T is [0, 800].

        Keyword arguments:
            temperatures (list/np.ndarray): list or array of temperatures to plot.
                If the array/list has length 2, use these as the start and endpoints
                with 21 plotting points.
            ax (matplotlib.pyplot.Axis): axis object to plot onto.

        """
        from matador.plotting.temperature_plotting import plot_free_energy

        return plot_free_energy(self, temperatures=temperatures, ax=ax, **kwargs)


class ElectronicDOS(DensityOfStates):
    """Specific class for electronic DOS data."""

    gaussian_width = 0.01
