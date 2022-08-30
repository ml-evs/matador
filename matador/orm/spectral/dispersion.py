# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements classes to store and manipulate electronic and
vibrational bandstructures, with or without projection data.

"""

import numpy as np
from matador.orm.spectral.spectral import Spectral
from matador.utils.chem_utils import INVERSE_CM_TO_EV, KELVIN_TO_EV

EPS = 1e-4


class Dispersion(Spectral):
    """Parent class for continuous spectra in reciprocal space, i.e.
    electronic and vibrational bandstructures.

    Note:
        This class speaks of "k-points" as general reciprocal space points
        used to display the dispersion curves; these correspond to CASTEP's
        phonon_kpoints or spectral_kpoints, and not the k-points
        used to generate the underlying wavefunction or dynamical matrix.

    """

    def find_full_kpt_branch(self):
        """Find all branch indices from branch start indices."""
        branch_inds = []
        for ind, start_ind in enumerate(self.kpoint_branch_start[:-1]):
            branch_inds.append(
                list(range(start_ind, self.kpoint_branch_start[ind + 1]))
            )
        branch_inds.append(list(range(self.kpoint_branch_start[-1], self.num_kpoints)))

        if not sum([len(branch) for branch in branch_inds]) == self.num_kpoints:
            raise RuntimeError(
                "Error parsing kpoints: number of kpoints does "
                "not match number in branches"
            )

        return branch_inds

    def set_branches_and_spacing(self):
        """Set the relevant kpoint spacing and branch attributes."""
        branch_start, spacing = self.find_kpoint_branches()
        self._data["kpoint_path_spacing"] = spacing
        self._data["kpoint_branch_start"] = branch_start

    def find_kpoint_branches(self):
        """Separate a kpoint path into discontinuous branches,

        Returns:
            list[list[int]]: list of lists containing the indices of
                the discontinuous kpoint branches.
            float: estimated k-point spacing from median of their
                separations.

        """
        kpt_diffs = np.linalg.norm(np.diff(self.kpoint_path_cartesian, axis=0), axis=1)
        spacing = np.median(kpt_diffs)
        # add 0 as its the start of the first path, then add all indices
        # have to add 1 to the where to get the start rather than end of the branch
        branch_start = [0] + (np.where(kpt_diffs > 3 * spacing)[0] + 1).tolist()
        return branch_start, spacing

    def linearise_path(self, preserve_kspace_distance=False):
        """For a given k-point path, normalise the spacing between points, mapping
        it onto [0, 1].

        Keyword arguments:
            preserve_kspace_distance (bool): if True, point separation
                will be determined by their actual separation in
                reciprocal space, otherwise they will be evenly spaced.

        Returns:
            np.ndarray: 3xN array containing k-points mapped onto [0, 1].

        """
        path = [0]
        for branch in self.kpoint_branches:
            for ind, kpt in enumerate(self.kpoint_path[branch]):
                if ind != len(branch) - 1:
                    if preserve_kspace_distance:
                        diff = np.sqrt(
                            np.sum((kpt - self.kpoint_path[branch[ind + 1]]) ** 2)
                        )
                    else:
                        diff = 1.0
                    path.append(path[-1] + diff)
        path = np.asarray(path)
        path /= np.max(path)
        if len(path) != self.num_kpoints - len(self.kpoint_branches) + 1:
            raise RuntimeError("Linearised kpoint path has wrong number of kpoints!")

        return path

    def reorder_bands(self):
        """Reorder the bands of this Dispersion object directly."""
        self._data["eigs_s_k"] = self.get_band_reordering(
            self.eigs, self.kpoint_branches
        )

    @staticmethod
    def get_band_reordering(eigs, kpoint_branches):
        """Recursively reorder eigenvalues such that bands join up correctly,
        based on local gradients.

        Parameters:
            dispersion (numpy.ndarray): array containing eigenvalues as a
                function of q/k
            branches (:obj:`list` of :obj:`int`): list containing branches of
                k/q-point path

        Returns:
            numpy.ndarray: reordered branches.

        """

        sorted_eigs = np.array(eigs, copy=True)
        num_bands = np.shape(sorted_eigs)[1]

        for channel_ind, channel in enumerate(eigs):
            eigs = channel
            for branch_ind, branch in enumerate(kpoint_branches):
                eigs_branch = channel[:, branch]
                converged = False
                counter = 0
                i_cached = 0
                while not converged and counter < len(branch):
                    counter += 1
                    for i in range(i_cached + 1, len(branch) - 1):
                        guess = 2 * eigs_branch[:, i] - eigs_branch[:, i - 1]
                        argsort_guess = np.argsort(guess)
                        if np.any(
                            np.argsort(guess) != np.argsort(eigs_branch[:, i + 1])
                        ):
                            tmp_copy = np.array(channel, copy=True)
                            for ind, mode in enumerate(
                                np.argsort(eigs_branch[:, i]).tolist()
                            ):
                                eigs_branch[mode, i + 1 :] = tmp_copy[:, branch][
                                    argsort_guess[ind], i + 1 :
                                ]
                            for other_branch in kpoint_branches[branch_ind:]:
                                eigs_other_branch = channel[:, other_branch]
                                for ind, mode in enumerate(
                                    np.argsort(channel[:, i]).tolist()
                                ):
                                    eigs_other_branch[mode] = tmp_copy[:, other_branch][
                                        argsort_guess[ind]
                                    ]
                            channel[:, other_branch] = eigs_other_branch
                            channel[:, branch] = eigs_branch
                            i_cached = i
                            break
                    else:
                        converged = True

            sorted_eigs[channel_ind] = channel.reshape(1, num_bands, len(channel[0]))

        return sorted_eigs

    def plot_dispersion(self, **kwargs):
        """Make a plot of the band structure, with projections, if found."""
        from matador.plotting.spectral_plotting import plot_spectral

        _kwargs = {
            "plot_dos": False,
            "plot_bandstructure": True,
            "plot_pdis": "projector_weights" in self,
            "phonons": "Vibrational" in self.__class__.__name__,
        }
        _kwargs.update(kwargs)
        plot_spectral(self, **_kwargs)

    def _reshaped_eigs(self, eigs, shape):
        """Attempts to reshape the eigenvalues into the desired shape.

        Parameters:
            eigs (np.ndarray): the eigs to reshape.
            shape (tuple): the desired shape.

        Returns:
            np.ndarray: the reshaped eigs.

        """
        raise NotImplementedError(
            "Wrong eigenvalue shape passed, and reshape function is not yet implemented. "
            "Eigs should have shape {}, not {}".format(shape, np.shape(eigs))
        )


class ElectronicDispersion(Dispersion):
    """Class that stores electronic dispersion data. Attributes are
    all implemented as properties based on underlying raw data.

    Attributes:
        source (list): list of source files.
        num_kpoints (int): number of kpoints.
        num_spins (int): number of spin channels.
        num_bands (int): number of bands.
        num_electrons (int): number of bands.
        eigs_s_k (numpy.ndarray):  eigenvalue array of shape
            (num_spins, num_bands, num_kpoints).
        kpoint_path (numpy.ndarray): array of shape (num_kpoints, 3)
            containing the k-point path in fractional coordinates.
        kpoint_path_cartesian (numpy.ndarray): as above, in Cartesian
            coordinates.
        fermi_energy (float): Fermi energy in eV (takes average of spin
            channels if more than one is present).
        spin_fermi_energy (list[float]): Fermi energy for each spin channel
            in eV.
        band_gap (float): smallest band gap across spin channels/momenta.
        spin_band_gap (list[float]): band gap per spin channel.
        projectors (list[tuple]): list of projector labels in format
            (`element`, `l-channel`).
        projector_weights (numpy.ndarray): array of projector_weights with shape
            (num_kpoints, num_bands, num_projectors)

    Note:
        projector_weights will only work with one spin channel.

    """

    required_keys = [
        "num_kpoints",
        "num_spins",
        "num_bands",
        "num_electrons",
        "eigs_s_k",
        "kpoint_path",
        "lattice_cart",
        "fermi_energy",
    ]

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        if kwargs.get("projection_data") is not None:
            projection_data = kwargs.get("projection_data")
            self._data["projectors"] = projection_data["projectors"]
            self._data["projector_weights"] = projection_data["projector_weights"]

        else:
            self._data["projectors"] = self._data.get("projectors")
            self._data["projector_weights"] = self._data.get("projector_weights")

        if self._data.get("projectors") is not None:
            if self.num_spins != 1:
                raise NotImplementedError(
                    "Projected dispersion not implemented" " for multiple spin channels"
                )
            # only want to take projectors and projector_weights from this data
            proj_shape = np.shape(self._data["projector_weights"])
            expected_shape = (self.num_kpoints, self.num_bands, self.num_projectors)
            if proj_shape != expected_shape:
                raise RuntimeError(
                    f"Incompatible shape of projector weights: {proj_shape}, was expecting {expected_shape}"
                )

        shape = (self.num_spins, self.num_bands, self.num_kpoints)
        if np.shape(self._data["eigs_s_k"]) != shape:
            self._data["eigs_s_k"] = self._reshaped_eigs(self._data["eigs_s_k"], shape)

    @property
    def num_spins(self):
        """Number of spin channels in spectrum."""
        return self._data["num_spins"]

    @property
    def num_electrons(self):
        """Number of electrons."""
        return self._data["num_electrons"]

    @property
    def eigs_s_k(self):
        """Array of electronic eigenvalues with shape
        (num_spins, num_bands, num_kpoints).

        """
        return self._data["eigs_s_k"]

    @property
    def eigs(self):
        """Alias for `self.eigs_s_k`."""
        return self.eigs_s_k

    @property
    def fermi_energy(self):
        """Return the Fermi energy as described in the raw data."""
        return self._data["fermi_energy"] or np.mean(
            self._data.get("spin_fermi_energy")
        )

    @property
    def spin_fermi_energy(self):
        """Return the Fermi energy as described in the raw data."""
        return self._data.get("spin_fermi_energy")

    @property
    def band_gap(self):
        """Return the band gap of the system."""
        if not self._data.get("band_gap"):
            self.set_gap_data()
        return self._data["band_gap"]

    @property
    def band_gap_path_inds(self):
        """Return the indices of the k-points that comprise the smallest
        band gap.

        """
        if not self._data.get("band_gap_path_inds"):
            self.set_gap_data()
        return self._data["band_gap_path_inds"]

    @property
    def spin_band_gap(self):
        """Return the band gap for each spin channel."""
        if not self._data.get("spin_band_gap"):
            self.set_gap_data()
        return self._data["spin_band_gap"]

    @property
    def spin_band_gap_path_inds(self):
        """Return the indices of the k-points that comprise the smallest
        band gap for each spin channel.

        """
        if not self._data.get("spin_band_gap_path_inds"):
            self.set_gap_data()
        return self._data["spin_band_gap_path_inds"]

    def set_gap_data(self):
        """Loop over bands to set the band gap, VBM, CBM, their
        positions and the smallest direct gap inside self._data,
        for each spin channel. Sets self.band_gap to be the smallest of
        the band gaps across all spin channels.

        """
        spin_keys = [
            "spin_band_gap",
            "spin_band_gap_path",
            "spin_band_gap_path_inds",
            "spin_valence_band_min",
            "spin_conduction_band_max",
            "spin_gap_momentum",
            "spin_direct_gap",
            "spin_direct_gap_path",
            "spin_direct_gap_path_inds",
            "spin_direct_valence_band_min",
            "spin_direct_conduction_band_max",
        ]

        for key in spin_keys:
            self._data[key] = self.num_spins * [None]

        for ispin in range(self.num_spins):
            vbm = -1e10
            cbm = 1e10
            cbm_pos = []
            vbm_pos = []
            # calculate indirect gap
            for _, branch in enumerate(self.kpoint_branches):
                for nb in range(self.num_bands):
                    band = (
                        self.eigs_s_k[ispin][nb][branch] - self.spin_fermi_energy[ispin]
                    )
                    argmin = np.argmin(band)
                    argmax = np.argmax(band)
                    if vbm + EPS < band[argmax] < 0:
                        vbm = band[argmax]
                        vbm_pos = [branch[argmax]]
                    elif vbm - EPS <= band[argmax] < 0:
                        vbm = band[argmax]
                        vbm_pos.extend([branch[argmax]])
                    if cbm - EPS > band[argmin] > 0:
                        cbm = band[argmin]
                        cbm_pos = [branch[argmin]]
                    elif cbm + EPS >= band[argmin] > 0:
                        cbm = band[argmin]
                        cbm_pos.extend([branch[argmin]])
                    if band[argmin] + EPS / 2 < 0 < band[argmax] - EPS / 2:
                        vbm = 0
                        cbm = 0
                        vbm_pos = [0]
                        cbm_pos = [0]
                        break

            smallest_diff = 1e10
            for _cbm_pos in cbm_pos:
                for _vbm_pos in vbm_pos:
                    if abs(_vbm_pos - _cbm_pos) < smallest_diff:
                        tmp_cbm_pos = _cbm_pos
                        tmp_vbm_pos = _vbm_pos
                        smallest_diff = abs(_vbm_pos - _cbm_pos)
            cbm_pos = tmp_cbm_pos
            vbm_pos = tmp_vbm_pos
            self._data["spin_valence_band_min"][ispin] = (
                vbm + self.spin_fermi_energy[ispin]
            )
            self._data["spin_conduction_band_max"][ispin] = (
                cbm + self.spin_fermi_energy[ispin]
            )
            self._data["spin_band_gap"][ispin] = cbm - vbm
            self._data["spin_band_gap_path"][ispin] = [
                self.kpoint_path[cbm_pos],
                self.kpoint_path[vbm_pos],
            ]
            self._data["spin_band_gap_path_inds"][ispin] = [cbm_pos, vbm_pos]
            self._data["spin_gap_momentum"][ispin] = np.linalg.norm(
                self.kpoint_path_cartesian[cbm_pos]
                - self.kpoint_path_cartesian[vbm_pos]
            )

            # calculate direct gap
            direct_gaps = np.zeros(self.num_kpoints)
            direct_cbms = np.zeros(self.num_kpoints)
            direct_vbms = np.zeros(self.num_kpoints)
            for ind, _ in enumerate(self.kpoint_path):
                direct_cbm = 1e10
                direct_vbm = -1e10
                for nb in range(self.num_bands):
                    band_eig = (
                        self.eigs_s_k[ispin][nb][ind] - self.spin_fermi_energy[ispin]
                    )
                    if direct_vbm <= band_eig < EPS:
                        direct_vbm = band_eig
                    if direct_cbm >= band_eig > EPS:
                        direct_cbm = band_eig
                direct_gaps[ind] = direct_cbm - direct_vbm
                direct_cbms[ind] = direct_cbm
                direct_vbms[ind] = direct_vbm
            self._data["spin_direct_gap"][ispin] = np.min(direct_gaps)
            self._data["spin_direct_conduction_band_max"][ispin] = (
                direct_cbms[np.argmin(direct_gaps)] + self.spin_fermi_energy[ispin]
            )
            self._data["spin_direct_valence_band_min"][ispin] = (
                direct_vbms[np.argmin(direct_gaps)] + self.spin_fermi_energy[ispin]
            )
            self._data["spin_direct_gap"][ispin] = np.min(direct_gaps)
            self._data["spin_direct_gap_path"][ispin] = 2 * [
                self.kpoint_path[np.argmin(direct_gaps)]
            ]
            self._data["spin_direct_gap_path_inds"][ispin] = 2 * [
                np.argmin(direct_gaps)
            ]

            if (
                np.abs(
                    self._data["spin_direct_gap"][ispin]
                    - self._data["spin_band_gap"][ispin]
                )
                < EPS
            ):
                self._data["spin_valence_band_min"][ispin] = (
                    direct_vbm + self.spin_fermi_energy[ispin]
                )
                self._data["spin_conduction_band_max"][ispin] = (
                    direct_cbm + self.spin_fermi_energy[ispin]
                )
                self._data["spin_band_gap_path_inds"][ispin] = self._data[
                    "spin_direct_gap_path_inds"
                ][ispin]
                cbm_pos = self._data["spin_direct_gap_path_inds"][ispin][0]
                vbm_pos = self._data["spin_direct_gap_path_inds"][ispin][1]
                self._data["spin_band_gap_path"][ispin] = self._data[
                    "spin_direct_gap_path"
                ][ispin]
                self._data["spin_gap_momentum"][ispin] = np.linalg.norm(
                    self.kpoint_path_cartesian[cbm_pos]
                    - self.kpoint_path_cartesian[vbm_pos]
                )

        spin_gap_index = np.argmin(self._data["spin_band_gap"])
        # use smallest spin channel gap data for standard non-spin data access
        for key in spin_keys:
            self._data[key.replace("spin_", "")] = self._data[key][spin_gap_index]

    def new_from_trimmed_path(self, k_start_ind=0, k_end_ind=None):
        """Returns a new ElectronicDispersion object with the kpoint
        path trimmed by the provided indices.

        """
        _new_data_dict = {}
        if k_end_ind is None:
            k_end_ind = len(self.kpoint_path)

        _new_data_dict["kpoint_path"] = self.kpoint_path[k_start_ind:k_end_ind]
        _new_data_dict["eigs_s_k"] = self.eigs[:, :, k_start_ind:k_end_ind]
        _new_data_dict["num_bands"] = self.num_bands
        _new_data_dict["num_electrons"] = self.num_electrons
        _new_data_dict["num_spins"] = self.num_spins
        _new_data_dict["num_kpoints"] = len(_new_data_dict["kpoint_path"])
        _new_data_dict["lattice_cart"] = self.lattice_cart
        _new_data_dict["fermi_energy"] = self.fermi_energy
        _new_data_dict["spin_fermi_energy"] = self.spin_fermi_energy

        return ElectronicDispersion(**_new_data_dict)


class VibrationalDispersion(Dispersion):
    """Class that stores vibrational dispersion data. Attributes are
    all implemented as properties based on underlying raw data.

    Attributes:
        source (list): list of source files.
        num_kpoints (int): number of kpoints.
        num_atoms (int): number of atoms.
        num_modes (int): number of phonon modes.
        eigs (numpy.ndarray):  eigenvalue array of shape
            (1, num_modes, num_kpoints), in frequency units below, with
            first index denoting the single "spin channel" for phonons.
        freq_unit (str): human-readable frequency unit used for eig array.
        kpoint_path (numpy.ndarray): array of shape (num_kpoints, 3)
            containing the k-point path in fractional coordinates.
        kpoint_path_cartesian (numpy.ndarray): as above, in Cartesian
            coordinates.

    """

    required_keys = ["num_kpoints", "num_modes", "eigs_q", "kpoint_path"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        shape = (1, self.num_modes, self.num_qpoints)
        if np.shape(self._data["eigs_q"]) != shape:
            self._data["eigs_q"] = self._reshaped_eigs(self._data["eigs_q"], shape)

    @property
    def num_atoms(self):
        """Number of atoms in cell."""
        return self._data["num_atoms"]

    @property
    def num_modes(self):
        """Number of phonon modes."""
        return self._data["num_modes"]

    @property
    def num_bands(self):
        """Alias for number of modes."""
        return self._data["num_modes"]

    @property
    def eigs_q(self):
        """Eigenvalues in frequency units `self.freq_unit`, with shape
        (1, num_modes, num_kpoints).

        """
        return np.asarray(self._data["eigs_q"])

    @property
    def softest_mode_freq(self):
        """The frequency of the softest mode in the calculation.
        Negative modes correspond to imaginary frequencies.

        """
        return np.min(self.eigs)

    @property
    def debye_temperature(self):
        """Returns the Debye temperature in K."""
        return self.debye_freq * INVERSE_CM_TO_EV / KELVIN_TO_EV

    @property
    def debye_freq(self):
        """Returns the Debye frequency in cm^-1."""
        return np.max(self.eigs)
