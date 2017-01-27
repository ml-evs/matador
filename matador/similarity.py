# coding: utf-8
""" This file defines various measures and ways of calculating
the similarity between two structures.

TO-DO:
    * standardize overlap integrals
    * otf calculation of required num_images
    * non-diagonal supercells
    * generic wrapper for a cursor to be screened

"""

# matador modules
from matador.utils.cell_utils import frac2cart
# external libraries
import numpy as np
# standard library
from itertools import product, combinations_with_replacement


class PDF(object):
    """ This class implements the calculation and comparison of pair
    distribution functions.
    """
    def __init__(self, doc, **kwargs):
        """ Initialise parameters.

        dr             : bin width for PDF (Angstrom)
        gaussian_width : width of Gaussian smearing (Angstrom)
        num_images     : number of unit cell images include in PDF calculation
        rmax           : maximum distance cutoff for PDF (Angstrom)

        """
        if kwargs.get('dr') is None:
            self.dr = 0.001
        else:
            self.dr = kwargs['dr']
        if kwargs.get('gaussian_width') is None:
            self.gaussian_width = 0.01
        else:
            self.gaussian_width = kwargs['gaussian_width']
        if kwargs.get('num_images') is None:
            self.num_images = 1
        else:
            self.num_images = kwargs['num_images']
        if kwargs.get('rmax') is None:
            self.rmax = 10
        else:
            self.rmax = kwargs['rmax']
        self.r_space = np.arange(0, self.rmax, self.dr)
        self.Gr = np.zeros((int(self.rmax / self.dr)))
        self.elem_Gr = dict()
        for comb in combinations_with_replacement(set(doc['atom_types']), 2):
            self.elem_Gr[tuple(set(comb))] = np.zeros((int(self.rmax / self.dr)))
        self.lattice = np.asarray(doc['lattice_cart'])
        self.atoms = frac2cart(doc['lattice_cart'], doc['positions_frac'])
        self.types = doc['atom_types']
        self.num_atoms = len(self.atoms)
        self.image_atoms = np.copy(self.atoms)
        if not kwargs.get('lazy'):
            self._calc_pdf()

    def _calc_pdf(self):
        """ Calculate PDF of a matador document.

        TO-DO: vectorise and element-projected.

        """
        for i in range(self.num_atoms):
            for j in range(i+1, self.num_atoms):
                d_ij = np.sqrt(np.sum((self.atoms[i] - self.atoms[j])**2))
                if d_ij <= self.rmax:
                    self.Gr += 2*np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / self.num_atoms
                    self.elem_Gr[tuple(set((self.types[i], self.types[j])))] += 2*np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / (self.num_atoms * (self.num_images+1)**3)
        # iterate over image cells
        trans = np.zeros((3))
        for prod in product(range(-self.num_images, self.num_images+1), repeat=3):
            if prod == (0, 0, 0):
                continue
            trans = 0
            for ind, multi in enumerate(prod):
                trans += self.lattice[ind] * multi
            for i in range(self.num_atoms):
                for j in range(self.num_atoms):
                    d_ij = np.sqrt(np.sum((self.atoms[i] - self.atoms[j] - trans)**2))
                    if d_ij < self.rmax:
                        self.Gr += np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / self.num_atoms
                        self.elem_Gr[tuple(set((self.types[i], self.types[j])))] += np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / (self.num_atoms * (self.num_images+1)**3)
        return

    def plot_projected(self):
        """ Plot projected PDFs. """
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(111)
        ax1.plot(self.r_space, self.Gr, lw=0.5, c='k', ls='--', label='total')
        for key in self.elem_Gr:
            ax1.plot(self.r_space, self.elem_Gr[key], label=key)
        ax1.legend()
        return


class PDFOverlap(object):
    def __init__(self, pdf_A, pdf_B):
        self.pdf_A = pdf_A
        self.pdf_B = pdf_B
        self.pdf_overlap()

    def pdf_overlap(self):
        """ Calculate the overlap of two PDFs via
        a simple meshed sum of their difference.
        """
        self.fine_space = np.arange(0, self.pdf_A.rmax, self.pdf_A.dr/2.0)
        self.fine_Gr_A = np.interp(self.fine_space, self.pdf_A.r_space, self.pdf_A.Gr)
        self.fine_Gr_B = np.interp(self.fine_space, self.pdf_B.r_space, self.pdf_B.Gr)
        self.overlap_fn = self.fine_Gr_A - self.fine_Gr_B
        self.overlap_int = np.trapz(np.abs(self.overlap_fn), dx=self.pdf_A.dr/2.0)

    def pdf_convolve(self):
        """ Calculate the convolution of two PDFs.
        """
        self.convolution = np.convolve(self.pdf_A.Gr, self.pdf_B.Gr)

    def pdf_diff_plot(self):
        """ Simple plot for comparing two PDF's. """
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(self.pdf_A.r_space, self.pdf_A.Gr)
        ax1.plot(self.pdf_B.r_space, self.pdf_B.Gr)
        ax1.set_xlabel('$r$ (Angstrom)')
        ax1.set_ylabel('$G(r)$')
        ax2.axhline(0, ls='--', c='k', lw=0.5)
        ax2.plot(self.fine_space, self.overlap_fn, ls='--')
        ax2.set_ylim(-0.5*ax1.get_ylim()[1], 0.5*ax1.get_ylim()[1])
        ax2.set_xlabel('$r$ (Angstrom)')
        ax2.set_ylabel('$G(r)$')
        return
