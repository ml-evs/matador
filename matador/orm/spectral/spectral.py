import numpy as np
from matador.utils.cell_utils import real2recip, frac2cart
from matador.orm.orm import DataContainer


class Spectral(DataContainer):
    """
    Note:
        This class speaks of "k-points" as general reciprocal space points
        used to display the dispersion curves; these correspond to CASTEP's
        phonon_kpoints or spectral_kpoints, and not the k-points
        used to generate the underlying wavefunction or dynamical matrix.

    """

    @property
    def eigs(self):
        """ Alias for the correct eigenvalue array. """
        if 'Vibrational' in self.__class__.__name__:
            return self._data['eigs_q']
        return self._data['eigs_s_k']

    @property
    def lattice_cart(self):
        """ The Cartesian lattice vectors of the real space lattice. """
        return self._data['lattice_cart']

    @property
    def num_kpoints(self):
        """ Number of dispersion k-points sampled. """
        return self._data['num_kpoints']

    @property
    def num_qpoints(self):
        """ Alias for number of kpoints. """
        return self.num_kpoints

    @property
    def projectors(self):
        """ Return list of projector labels in the format
        `(element, l-channel)`.

        """
        return self._data.get('projectors')

    @property
    def num_modes(self):
        """ Number of eigenvalues per q/k-point. """
        return self._data['num_modes']

    @property
    def num_bands(self):
        """ Number of eigenvalues per q/k-point. """
        if 'Vibrational' in self.__class__.__name__:
            return self._data['num_modes']
        return self._data['num_bands']

    @property
    def projector_weights(self):
        """ Return the array of projector weights per eigval, with shape
        (num_projectors, num_kpoints, num_bands).

        """
        return self._data.get('projector_weights')

    @property
    def num_projectors(self):
        """ Return the number of projectors. """
        if self.projectors is None:
            return 0
        return len(self.projectors)

    @property
    def kpoint_branches(self):
        """ Return the k-point branches in the older format, which
        contained a list of lists of continous indices.

        """
        if self._data.get('kpoint_branches') is None:
            self._data['kpoint_branches'] = self.find_full_kpt_branch()
        return self._data['kpoint_branches']

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
    def kpoint_path(self):
        """ The fractional sampling path in reciprocal space. """
        return np.asarray(self._data['kpoint_path'])

    @property
    def kpoint_weights(self):
        if 'kpoint_weights' in self._data:
            return np.asarray(self._data['kpoint_weights'])
        return None

    @property
    def kpoint_path_cartesian(self):
        """ The reicprocal space sampling path in Cartesian coordinates. """
        return np.asarray(frac2cart(real2recip(self.lattice_cart),
                                    self.kpoint_path))

    @property
    def num_spins(self):
        """ Dummy number of spins. """
        return 1

    @property
    def spin_fermi_energy(self):
        """ Dummy Fermi energy per spin channel. """
        return [0]
