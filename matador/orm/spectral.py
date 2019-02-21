# coding: utf-8
# Distributed under the terms of the MIT license.

""" This file implements the Spectrum class and its
children, used for storing collections of electronic
and vibrational reciprocal spcae spectra, with or without
projector_weights.

"""

import copy
import numpy as np
from matador.utils.cell_utils import real2recip, frac2cart

EPS = 1e-6


class Spectrum:
    """ Container object for different types of electronic and vibrational
    spectra, e.g. DOS and dispersion.

    """
    # TODO write this class
    def __init__(self):
        raise NotImplementedError


class Dispersion:
    """ Parent class for continuous spectra in reciprocal space, i.e.
    electronic and vibrational bandstructures.

    """
    def __init__(self, data):
        """ Initalise copy of raw data. """
        self._data = copy.deepcopy(data)

    def __getitem__(self, key):
        """ Allow properties to be used with key access.

        Raises a KeyError if property or key can't be found.

        """
        try:
            return getattr(self, key)
        except AttributeError:
            pass

        return self._data[key]


class ElectronicDispersion(Dispersion):
    """ Class that stores electronic dispersion data. Attributes are
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
    def __init__(self, data, projection_data=None):
        super().__init__(data)
        if projection_data is not None:
            if self.num_spins != 1:
                raise NotImplementedError('Projected dispersion not implemented'
                                          ' for multiple spin channels')
            # only want to take projectors and projector_weights from this data
            if self.num_kpoints != len(projection_data['kpoints']):
                raise RuntimeError('k-point mismatch between bandstructure '
                                   'and projection data.')

            self._data['projectors'] = projection_data['projectors']
            self._data['projector_weights'] = projection_data['pdis']
            assert np.shape(projection_data['pdis']) == (self.num_kpoints, self.num_bands, self.num_projectors)

        else:
            self._data['projectors'] = None
            self._data['projector_weights'] = None

    @property
    def source(self):
        """ List of source files. """
        return self._data['source']

    @property
    def num_kpoints(self):
        """ Number of k-points sampled. """
        return self._data['num_kpoints']

    @property
    def num_spins(self):
        """ Number of spin channels in spectrum. """
        return self._data['num_spins']

    @property
    def num_bands(self):
        """ Number of eigenvalues per spin, per k-point. """
        return self._data['num_bands']

    @property
    def num_electrons(self):
        """ Number of electrons. """
        return self._data['num_electrons']

    @property
    def eigs_s_k(self):
        """ Array of electronic eigenvalues with shape
        (num_spins, num_bands, num_kpoints).

        """
        return self._data['eigs_s_k']

    @property
    def kpoint_path(self):
        """ The fractional sampling path in reciprocal space. """
        return self._data['kpoint_path']

    @property
    def kpoint_path_cartesian(self):
        """ The reicprocal space sampling path in Cartesian coordinates. """
        return np.asarray(frac2cart(real2recip(self.lattice_cart),
                                    self.kpoint_path))

    @property
    def lattice_cart(self):
        """ The Cartesian lattice vectors of the real space lattice. """
        return self._data['lattice_cart']

    @property
    def fermi_energy(self):
        """ Return the Fermi energy as described in the raw data. """
        return self._data['fermi_energy'] or np.mean(self._data.get('spin_fermi_energy'))

    @property
    def spin_fermi_energy(self):
        """ Return the Fermi energy as described in the raw data. """
        return self._data.get('spin_fermi_energy')

    @property
    def kpoint_branch_start(self):
        """ Return the indices of the start of branches. """
        if not self._data.get('kpoint_branch_start'):
            self.set_branches_and_spacing()
        return self._data['kpoint_branch_start']

    @property
    def kpoint_path_spacing(self):
        """ An estimated kpoint spacing. """
        if not self._data.get('kpoint_path_spacing'):
            self.set_branches_and_spacing()
        return self._data['kpoint_path_spacing']

    @property
    def kpoint_branches(self):
        """ Return the k-point branches in the older format, which
        contained a list of lists of continous indices.

        """
        if not self._data.get('kpoint_branches'):
            self._data['kpoint_branches'] = self.find_full_kpt_branch()
        return self._data['kpoint_branches']

    @property
    def band_gap(self):
        """ Return the band gap of the system. """
        if not self._data.get('band_gap'):
            self.set_gap_data()
        return self._data['band_gap']

    @property
    def band_gap_path_inds(self):
        """ Return the indices of the k-points that comprise the smallest
        band gap.

        """
        if not self._data.get('band_gap_path_inds'):
            self.set_gap_data()
        return self._data['band_gap_path_inds']

    @property
    def spin_band_gap(self):
        """ Return the band gap for each spin channel. """
        if not self._data.get('spin_band_gap'):
            self.set_gap_data()
        return self._data['spin_band_gap']

    @property
    def spin_band_gap_path_inds(self):
        """ Return the indices of the k-points that comprise the smallest
        band gap for each spin channel.

        """
        if not self._data.get('spin_band_gap_path_inds'):
            self.set_gap_data()
        return self._data['spin_band_gap_path_inds']

    @property
    def projectors(self):
        """ Return list of projector labels in the format
        (`element`, `l-channel`).

        """
        return self._data['projectors']

    @property
    def projector_weights(self):
        """ Return the array of projector weights per eigval, with shape
        (num_projectors, num_kpoints, num_bands).

        """
        return self._data['projector_weights']

    @property
    def num_projectors(self):
        """ Return the number of projectors. """
        if self.projectors is None:
            return 0
        return len(self.projectors)

    def find_full_kpt_branch(self):
        """ Find all branch indices from branch start indices. """
        branch_inds = []
        for ind, start_ind in enumerate(self.kpoint_branch_start[:-1]):
            branch_inds.append(list(range(start_ind,
                                          self.kpoint_branch_start[ind+1])))
        branch_inds.append(list(range(self.kpoint_branch_start[-1],
                                      self.num_kpoints)))

        if not sum([len(branch) for branch in branch_inds]) == self.num_kpoints:
            raise RuntimeError('Error parsing kpoints: number of kpoints does '
                               'not match number in branches')

        return branch_inds

    def set_branches_and_spacing(self):
        """ Set the relevant kpoint spacing and branch attributes. """
        branch_start, spacing = self.find_kpoint_branches()
        self._data['kpoint_path_spacing'] = spacing
        self._data['kpoint_branch_start'] = branch_start

    def find_kpoint_branches(self):
        """ Separate a kpoint path into discontinuous branches,

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
        branch_start = [0] + (np.where(kpt_diffs > 3*spacing)[0] + 1).tolist()
        return branch_start, spacing

    def set_gap_data(self):
        """ Loop over bands to set the band gap, VBM, CBM, their
        positions and the smallest direct gap inside self._data,
        for each spin channel. Sets self.band_gap to be the smallest of
        the band gaps across all spin channels.

        """
        spin_keys = ['spin_band_gap', 'spin_band_gap_path', 'spin_band_gap_path_inds',
                     'spin_valence_band_min', 'spin_conduction_band_max', 'spin_gap_momentum',
                     'spin_direct_gap', 'spin_direct_gap_path', 'spin_direct_gap_path_inds',
                     'spin_direct_valence_band_min', 'spin_direct_conduction_band_max']

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
                    band = self.eigs_s_k[ispin][nb][branch] - self.spin_fermi_energy[ispin]
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
                    if band[argmin] < 0 < band[argmax]:
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
            self._data['spin_valence_band_min'][ispin] = vbm + self.spin_fermi_energy[ispin]
            self._data['spin_conduction_band_max'][ispin] = cbm + self.spin_fermi_energy[ispin]
            self._data['spin_band_gap'][ispin] = cbm - vbm
            self._data['spin_band_gap_path'][ispin] = [self.kpoint_path[cbm_pos], self.kpoint_path[vbm_pos]]
            self._data['spin_band_gap_path_inds'][ispin] = [cbm_pos, vbm_pos]
            self._data['spin_gap_momentum'][ispin] = np.linalg.norm(self.kpoint_path_cartesian[cbm_pos] - self.kpoint_path_cartesian[vbm_pos])

            # calculate direct gap
            direct_gaps = np.zeros(self.num_kpoints)
            direct_cbms = np.zeros(self.num_kpoints)
            direct_vbms = np.zeros(self.num_kpoints)
            for ind, _ in enumerate(self.kpoint_path):
                direct_cbm = 1e10
                direct_vbm = -1e10
                for nb in range(self.num_bands):
                    band_eig = self.eigs_s_k[ispin][nb][ind] - self.spin_fermi_energy[ispin]
                    if direct_vbm <= band_eig < 0:
                        direct_vbm = band_eig
                    if direct_cbm >= band_eig > 0:
                        direct_cbm = band_eig
                direct_gaps[ind] = direct_cbm - direct_vbm
                direct_cbms[ind] = direct_cbm
                direct_vbms[ind] = direct_vbm
            self._data['spin_direct_gap'][ispin] = np.min(direct_gaps)
            self._data['spin_direct_conduction_band_max'][ispin] = direct_cbms[np.argmin(direct_gaps)] + self.spin_fermi_energy[ispin]
            self._data['spin_direct_valence_band_min'][ispin] = direct_vbms[np.argmin(direct_gaps)] + self.spin_fermi_energy[ispin]
            self._data['spin_direct_gap'][ispin] = np.min(direct_gaps)
            self._data['spin_direct_gap_path'][ispin] = 2 * [self.kpoint_path[np.argmin(direct_gaps)]]
            self._data['spin_direct_gap_path_inds'][ispin] = 2 * [np.argmin(direct_gaps)]

            if np.abs(self._data['spin_direct_gap'][ispin] - self._data['spin_band_gap'][ispin]) < EPS:
                self._data['spin_valence_band_min'][ispin] = direct_vbm + self.spin_fermi_energy[ispin]
                self._data['spin_conduction_band_max'][ispin] = direct_cbm + self.spin_fermi_energy[ispin]
                self._data['spin_band_gap_path_inds'][ispin] = self._data['spin_direct_gap_path_inds'][ispin]
                cbm_pos = self._data['spin_direct_gap_path_inds'][ispin][0]
                vbm_pos = self._data['spin_direct_gap_path_inds'][ispin][1]
                self._data['spin_band_gap_path'][ispin] = self._data['spin_direct_gap_path'][ispin]
                self._data['spin_gap_momentum'][ispin] = np.linalg.norm(self.kpoint_path_cartesian[cbm_pos] - self.kpoint_path_cartesian[vbm_pos])

        spin_gap_index = np.argmin(self._data['spin_band_gap'])
        # use smallest spin channel gap data for standard non-spin data access
        for key in spin_keys:
            self._data[key.replace('spin_', '')] = self._data[key][spin_gap_index]


# TODO finish these
# class VibrationalDispersion(Dispersion):
    # def __init__(self):
        # pass

# class DensityOfStates:
    # def __init__(self, dos_dict=None):
        # pass

# class VibrationalDOS(DensityOfStates):
    # def __init__(self):
        # pass

# class ElectronicDOS(DensityOfStates):
    # def __init__(self):
        # pass
