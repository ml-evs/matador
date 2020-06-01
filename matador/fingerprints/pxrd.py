# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the PXRD class for simulating powder XRD pattern
of a crystal.

"""

import itertools
import os
from typing import Tuple

import numpy as np

from matador.fingerprints.fingerprint import Fingerprint, FingerprintFactory
from matador.crystal import Crystal
from matador.utils.cell_utils import standardize_doc_cell, real2recip
from matador.utils.chem_utils import get_formula_from_stoich

THETA_TOL = 1e-5


class PXRD(Fingerprint):
    """ This class for computes powder X-ray diffraction patterns  of a
    given crystal for a certain incident wavelength. The cell is
    standardised with spglib before computing PXRD.

    This calculation takes into account atomic scattering factors, Lorentz
    polarisation and thermal broadening (with Debye-Waller factors set to
    1). Note: this class does not perform any q-dependent peak broadening,
    and instead uses a simple Lorentzian broadening. The default width
    of 0.03 provides good agreement with e.g. GSAS-II's default CuKa setup.
    Only one wavelength can be used at a time, but multiple patterns could be
    combined post hoc.

    Attributes:
        self.peak_positions (numpy.ndarray): sorted peak positions as
            values in 2θ
        self.hkls (numpy.ndarray): Miller indices correspnding to peaks,
            sorted by peak angle.
        self.peak_intensities (numpy.ndarray): intensity of each peak.
        self.pattern (numpy.ndarray): Lorentzian-broadened pattern at
            values of self.two_thetas.
        self.two_thetas (numpy.ndarray): continuous space of 2θ values
            corresponding to sample points of self.pattern.

    """
    def __init__(
        self,
        doc,
        wavelength: float = 1.5406,
        lorentzian_width: float = 0.03,
        two_theta_resolution: float = 0.01,
        two_theta_bounds: Tuple[float, float] = (0, 90),
        theta_m: float = 0.0,
        scattering_factors: str = "RASPA",
        lazy=False,
        plot=False,
        progress=False,
        *args,
        **kwargs
    ):
        """ Set up the PXRD, and compute it, if lazy is False.

        Parameters:
            doc (dict/Crystal): matador document to compute PXRD for.

        Keyword arguments:
            lorentzian_width (float): width of Lorentzians for broadening (DEFAULT: 0.03)
            wavelength (float): incident X-ray wavelength
                (DEFAULT: CuKa, 1.5406).
            theta_m (float): the monochromator angle in degrees (DEFAULT: 0)
            two_theta_resolution (float): resolution of grid 2θ
                used for plotting.
            two_theta_bounds (tuple of float): values between which
                to compute the PXRD pattern.
            scattering_factors (str): either "GSAS" or "RASPA" (default),
                which set of atomic scattering factors to use.
            lazy (bool): whether to compute PXRD or just set it up.
            plot (bool): whether to display PXRD as a plot.

        """
        self.wavelength = wavelength
        self.lorentzian_width = lorentzian_width
        self.two_theta_resolution = two_theta_resolution
        self.two_theta_bounds = list(two_theta_bounds)
        self.theta_m = theta_m
        self.scattering_factors = scattering_factors
        self.progress = progress

        if self.two_theta_bounds[0] < THETA_TOL:
            self.two_theta_bounds[0] = THETA_TOL

        if np.min(doc.get('site_occupancy', [1.0])) < 1.0:
            print("System has partial occupancy, not refining with spglib.")
            self.doc = Crystal(doc)
        else:
            self.doc = Crystal(standardize_doc_cell(doc, primitive=True))

        self.formula = get_formula_from_stoich(self.doc['stoichiometry'], tex=True)
        self.spg = self.doc['space_group']

        species = list(set(self.doc['atom_types']))

        # this could be cached across PXRD objects but is much faster than the XRD calculation itself
        if self.scattering_factors == "GSAS":
            from matador.data import GSAS_ATOMIC_SCATTERING_COEFFS
            self.atomic_scattering_coeffs = {spec: GSAS_ATOMIC_SCATTERING_COEFFS[spec] for spec in species}
        elif self.scattering_factors == "RASPA":
            from matador.data import RASPA_ATOMIC_SCATTERING_COEFFS
            self.atomic_scattering_coeffs = {spec: RASPA_ATOMIC_SCATTERING_COEFFS[spec] for spec in species}
        else:
            raise RuntimeError(
                "No set of scattering factors matched: {}. Please use 'GSAS' or 'RASPA'."
                .format(self.scattering_factors)
            )

        if not lazy:
            self.calculate()
            if plot:
                self.plot()

    def calc_pxrd(self):
        """ Calculate the PXRD pattern. """

        # set crystallographic data
        lattice_abc = np.asarray(self.doc.lattice_abc)
        lattice_cart = np.asarray(self.doc.lattice_cart)
        positions_abs = np.asarray(self.doc.positions_abs)
        site_occupancies = np.asarray(self.doc.site_occupancies)

        # find allowed reciprocal lattice points within limiting sphere
        min_r, max_r = [2 / self.wavelength * np.sin(np.pi / 180 * t / 2) for t in self.two_theta_bounds]
        Ns = np.floor(max_r * lattice_abc[0, :]).astype(int)

        recip = np.asarray(real2recip(lattice_cart)).T
        qs = np.zeros((2*sum(Ns), 3), dtype=np.float64)
        hkls = np.asarray(
            list(itertools.product(
                range(-Ns[0], Ns[0] + 1),
                range(-Ns[1], Ns[1] + 1),
                range(-Ns[2], Ns[2] + 1)
            )), dtype=np.float64
        )
        hkls = hkls[np.argsort(np.linalg.norm(hkls, axis=-1))]
        qs = np.dot(recip, hkls.T).T

        # filter out by theta bounds
        q_mags = np.linalg.norm(qs, axis=-1)
        allowed = np.where(np.logical_and(q_mags <= max_r * 2 * np.pi, q_mags >= min_r * 2 * np.pi))
        qs = qs[allowed]
        hkls = hkls[allowed]

        # compute Bragg condition to find peak locations
        sin_tau = np.linalg.norm(qs, axis=1) * self.wavelength / (4 * np.pi)
        sin_tau[sin_tau > 1] = 0
        taus = 2 * np.arcsin(sin_tau)

        # compute structure factor S(q) as sum of atomic scattering factors
        S_q = np.zeros_like(taus)

        if self.progress:
            import tqdm
            bar = tqdm.tqdm
        else:
            def bar(x):
                return x

        for ind, q_vector in bar(enumerate(qs)):
            # accumulate atomic scattering factors
            atomic_factor = {}
            for species in set(self.doc.atom_types):
                atomic_factor[species] = self.atomic_scattering_factor(q_mags[ind], species)
            factors = np.array([atomic_factor[species] for species in self.doc.atom_types])
            F_s = np.sum(np.exp(1j * positions_abs @ q_vector) * factors * site_occupancies)
            S_q[ind] = np.abs(F_s)**2

        # apply Lorentz correction for polarisation and finite size effects
        S_q *= 2 * (1 + np.cos(taus) ** 2 * np.cos(2 * self.theta_m) ** 2) / (np.sin(taus) * np.sin(0.5 * taus))

        # thermal correction assuming no Debye-Waller factor
        S_q *= np.exp(-np.sin(taus)**2 / self.wavelength**2)**2
        S_q /= np.max(S_q)

        # create histogram and broaden onto 2 theta space in degrees
        self.peak_positions = (180 / np.pi) * taus
        self.hkls = hkls
        self.two_thetas = np.arange(self.two_theta_bounds[0],
                                    self.two_theta_bounds[1] + self.two_theta_resolution,
                                    self.two_theta_resolution)

        self.pattern, bins = np.histogram(self.peak_positions, bins=self.two_thetas, weights=S_q)

        if self.lorentzian_width > 0:
            self.pattern = self._broadening_unrolled(
                self.pattern, self.two_thetas, self.lorentzian_width, broadening_type='lorentzian'
            )
        else:
            # shift and clip the last two theta value if we didnt do broadening
            self.two_thetas = self.two_thetas[:-1] + self.two_theta_resolution / 2

        self.pattern /= np.max(self.pattern)

        order = np.argsort(self.peak_positions)
        self.hkls = hkls[order]
        self.peak_intensities = S_q
        self.peak_positions = self.peak_positions[order]

        # alias old name for compatibility...
        self.spectrum = self.pattern

    def calculate(self):
        """ Alias for calculating the PXRD pattern. """
        self.calc_pxrd()

    def atomic_scattering_factor(self, q_mag, species):
        """ Return fit for particular atom at given q-vector.

        Parameters:
            q_mag (float): magnitude of the q_vector.
            species (str): the element label.

        Returns:
            float: the atomic scattering factor.

        """
        a = self.atomic_scattering_coeffs[species][0]
        b = self.atomic_scattering_coeffs[species][1]
        c = self.atomic_scattering_coeffs[species][2]
        return c + np.sum(a * np.exp(-b * (q_mag / (4 * np.pi))**2))

    def plot(self, **kwargs):
        """ Wrapper function to plot the PXRD pattern. """
        from matador.plotting.pxrd_plotting import plot_pxrd
        plot_pxrd(self, **kwargs)

    def save_pattern(pxrd, fname):
        """ Write a file to `fname` that contains the xy coordinates of
        the PXRD pattern.

        """

        if os.path.isfile(fname):
            raise RuntimeError(f"Requested filename {fname} already exists!")

        from matador import __version__
        header = f""" PXRD pattern computed with matador {__version__}.
Input file:
{pxrd.doc.source[0]}

Structure:
{pxrd.doc}

Settings:
wavelength = {pxrd.wavelength} Å
theta_m = {pxrd.theta_m} degrees
lorentzian_width = {pxrd.lorentzian_width} degrees

2θ (degrees),\t\t\tRelative intensity"""

        np.savetxt(fname, np.vstack([pxrd.two_thetas, pxrd.pattern]).T, header=header, fmt='%.14e', delimiter='\t')

    def save_peaks(pxrd, fname):
        """ Write a file to `fname` that contains the peak list. """

        if os.path.isfile(fname):
            raise RuntimeError(f"Requested filename {fname} already exists!")

        from matador import __version__
        header = f""" PXRD peaks computed with matador {__version__}.
Input file:
{pxrd.doc.source[0]}

Structure:
{pxrd.doc}

Settings:
wavelength = {pxrd.wavelength} Å
theta_m = {pxrd.theta_m} degrees
lorentzian_width = {pxrd.lorentzian_width} degrees

<hkl>,\t\tPeak position (degrees)"""

        np.savetxt(
            fname,
            np.vstack(
                [pxrd.hkls[:, 0], pxrd.hkls[:, 1], pxrd.hkls[:, 2], pxrd.peak_positions]
            ).T,
            header=header,
            fmt=['% d', '% d', '% d', '%.14e'],
            delimiter='\t'
        )


class PXRDFactory(FingerprintFactory):
    fingerprint = PXRD
    default_key = 'pxrd'
