# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the PXRD class for simulating powder XRD spectrum
of a crystal.

"""

import copy
import os
import numpy as np

from matador.similarity.fingerprint import Fingerprint, FingerprintFactory
from matador.utils.cell_utils import frac2cart, standardize_doc_cell
from matador.utils.chem_utils import get_formula_from_stoich

THETA_TOL = 1e-5


class PXRD(Fingerprint):
    """ This class for computes powder X-ray diffraction patterns  of a
    given crystal for a certain incident wavelength. The cell is
    standardised with spglib before computing PXRD.

    Attributes:
        self.peak_positions (numpy.ndarray): peak positions as values
            in 2θ
        self.spectrum (numpy.ndarray): gaussian-broadened spectrum at
            values of self.two_thetas.
        self.two_thetas (numpy.ndarray): continuous space of 2θ values
            corresponding to sample points of self.spectrum.

    """
    # TODO: group peaks and return amplitudes and positions
    # TODO: vectorise loop over q-vectors
    # TODO: check atomic scattering data
    def __init__(self, doc, lazy=False, plot=False, *args, **kwargs):
        """ Set up the PXRD, and compute it, if lazy is False.

        Parameters:
            doc (dict): matador document to compute PXRD for.

        Keyword arguments:
            lazy (bool): whether to compute PXRD or just set it up.
            plot (bool): whether to display PXRD as a plot.
            gaussian_width (float): width of gaussians for broadening.
            wavelength (float): incident X-ray wavelength
            two_theta_resolution (float): resolution of grid 2θ
                used for plotting.

        """
        prop_defaults = {'wavelength': 1.5406,
                         'gaussian_width': 0.1,
                         'two_theta_resolution': 0.001}

        options = copy.deepcopy(prop_defaults)
        options.update(kwargs)
        self.wavelength = options['wavelength']
        self.gaussian_width = options['gaussian_width']
        self.two_theta_resolution = options['two_theta_resolution']

        self.doc = standardize_doc_cell(doc)
        self.formula = get_formula_from_stoich(self.doc['stoichiometry'], tex=True)
        self.spg = self.doc['space_group']

        species = list({species for species in self.doc['atom_types']})

        self.atomic_scattering_coeffs = {}
        data_file = os.path.dirname(os.path.realpath(__file__)) + '/../data/atomic_scattering_factors.dat'
        with open(data_file, 'r') as f:
            flines = f.readlines()
        for line in flines:
            x = line.split()
            label = x[0]
            if label in species:
                a_inds = [1, 3, 5, 7]
                b_inds = [2, 4, 6, 8]
                a = np.array([float(x[ind]) for ind in a_inds])
                b = np.array([float(x[ind]) for ind in b_inds])
                c = float(x[9])
                self.atomic_scattering_coeffs[label] = [a, b, c]

        if not lazy:
            self.calculate()
            if plot:
                self.plot()

    def calc_pxrd(self):
        two_theta_bounds = [10, 140]
        lattice_abc = np.array(self.doc['lattice_abc'])
        self.lattice_cart = lattice_cart = np.array(self.doc['lattice_cart'])

        # find allowed reciprocal lattice points within limiting sphere
        Ns = np.floor((2 / self.wavelength) * lattice_abc[0, :] *
                      np.sin((two_theta_bounds[1] / 2) * np.pi / 180)).astype(int)

        R_inv = np.linalg.inv(lattice_cart)
        qs = []
        for i in range(-Ns[0], Ns[0]+1):
            for j in range(-Ns[1], Ns[1]+1):
                for k in range(-Ns[2], Ns[2]+1):
                    qs.append(np.dot(R_inv, 2*np.pi * np.array([i, j, k])))
        qs = np.asarray(qs)

        # compute Bragg condition to find peak locations
        sin_tau = np.linalg.norm(qs, axis=1) * self.wavelength / (4 * np.pi)
        sin_tau[sin_tau > 1] = 0
        taus = 2 * np.arcsin(sin_tau)
        qs = qs[np.where(taus > THETA_TOL)]
        taus = taus[taus > THETA_TOL]

        # compute structure factor S(q) as sum of atomic scattering factors
        S_q = np.zeros_like(taus)
        for ind, q_vector in enumerate(qs):
            S_q[ind] = np.abs(self.structure_factor(q_vector))**2

        # apply Lorentz correction for polarisation and finite size effects
        S_q *= 0.5 * (1 + np.cos(taus) ** 2) / (np.sin(taus) * np.sin(0.5 * taus))

        # thermal correction: TODO: what is this correction? surely depends on Debye-Waller
        B = 1
        S_q *= np.exp(-B * np.sin(taus)**2 / self.wavelength**2)**2

        # create histogram and broaden onto 2 theta space in degrees
        self.peak_positions = (180 / np.pi) * taus
        self.two_thetas = np.arange(two_theta_bounds[0],
                                    two_theta_bounds[1] + self.two_theta_resolution,
                                    self.two_theta_resolution)
        self.spectrum, bins = np.histogram(self.peak_positions, bins=self.two_thetas, weights=S_q)
        if self.gaussian_width > 0:
            self.spectrum = self._broadening_unrolled(self.spectrum, self.two_thetas, self.gaussian_width)
        self.spectrum /= np.max(self.spectrum)

    def calculate(self):
        self.calc_pxrd()

    def structure_factor(self, q_vector):
        F_s = 0
        q_mag = np.linalg.norm(q_vector)
        for ind, species in enumerate(self.doc['atom_types']):
            phase = np.dot(q_vector, frac2cart(self.lattice_cart, self.doc['positions_frac'][ind]))
            atomic_factor = self.atomic_scattering_factor(q_mag, species)
            F_s += np.exp(1j * phase) * atomic_factor

        return F_s

    def atomic_scattering_factor(self, q_mag, species):
        # return fit for particular atom at given q
        a = self.atomic_scattering_coeffs[species][0]
        b = self.atomic_scattering_coeffs[species][1]
        c = self.atomic_scattering_coeffs[species][2]
        f = c
        for i in range(4):
            f += a[i] * np.exp(-b[i] * (q_mag / (4 * np.pi))**2)

        return f

    def plot(self):
        from matador.plotting.pxrd_plotting import plot_pxrd
        plot_pxrd(self)


class PXRDFactory(FingerprintFactory):
    fingerprint = PXRD
    default_key = 'pxrd'
