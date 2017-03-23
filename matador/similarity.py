# coding: utf-8
""" This file defines various measures and ways of calculating
the similarity between two structures.

TO-DO:
    * standardize overlap integrals
    * otf calculation of required num_images
    * non-diagonal supercells
    * generic wrapper for a cursor to be screened
    * comparing PDFs at different parameters

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

        dr             : bin width for PDF (Angstrom) (DEFAULT: 0.001)
        gaussian_width : width of Gaussian smearing (Angstrom) (DEFAULT: 0.05)
        num_images     : number of unit cell images include in PDF calculation (DEFAULT: 3)
        rmax           : maximum distance cutoff for PDF (Angstrom) (DEFAULT: 15)

        """
        if kwargs.get('dr') is None:
            self.dr = 0.001
        else:
            self.dr = kwargs['dr']
        if kwargs.get('gaussian_width') is None:
            self.gaussian_width = 0.05
        else:
            self.gaussian_width = kwargs['gaussian_width']
        if kwargs.get('num_images') is None:
            self.num_images = 3
        else:
            self.num_images = kwargs['num_images']
        if kwargs.get('rmax') is None:
            self.rmax = 15
        else:
            self.rmax = kwargs['rmax']
        self.r_space = np.arange(0, self.rmax, self.dr)
        self.Gr = np.zeros((int(self.rmax / self.dr)))
        self.elem_Gr = dict()
        for comb in combinations_with_replacement(set(doc['atom_types']), 2):
            self.elem_Gr[tuple(set(comb))] = np.zeros((int(self.rmax / self.dr)))
        self.lattice = np.asarray(doc['lattice_cart'])
        self.atoms = np.asarray(frac2cart(doc['lattice_cart'], doc['positions_frac']))
        self.types = doc['atom_types']
        self.label = ' '.join(doc['text_id'])
        self.num_atoms = len(self.atoms)
        self.volume = doc['cell_volume']
        self.image_atoms = np.copy(self.atoms)
        if not kwargs.get('lazy'):
            self._calc_pdf()

    def _calc_pdf(self):
        """ Calculate PDF of a matador document.

        TO-DO: vectorise.

        """
        for i in range(self.num_atoms):
            for j in range(i+1, self.num_atoms):
                d_ij = np.sqrt(np.sum((self.atoms[i] - self.atoms[j])**2))
                if d_ij <= self.rmax:
                    self.Gr += 2*np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / (self.num_atoms * (self.num_images+1)**3)
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
                        self.Gr += np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / (self.num_atoms * (self.num_images+1)**3)
                        self.elem_Gr[tuple(set((self.types[i], self.types[j])))] += np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / (self.num_atoms * (self.num_images+1)**3)
        return

    def get_sim_distance(self, pdf_B):
        """ Return the similarity between two PDFs. """
        return PDFOverlap(self, pdf_B, rescale='bond').similarity_distance

    def plot_projected(self, keys=None):
        """ Plot projected PDFs. """
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(111)
        ax1.plot(self.r_space, self.Gr, lw=1, ls='--', label='total')
        if keys is None:
            keys = [key for key in self.elem_Gr]
        for key in keys:
            ax1.plot(self.r_space, self.elem_Gr[key], label='-'.join(key))
        ax1.legend(loc=1)
        ax1.set_ylabel('$g(r)$')
        ax1.set_xlabel('$r$ (Angstrom)')
        return

    def plot_pdf(self):
        """ Plot projected PDFs. """
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(111)
        ax1.plot(self.r_space, self.Gr, lw=2)
        ax1.set_ylabel('$g(r)$')
        ax1.set_xlabel('$r$ (Angstrom)')
        return


class PDFOverlap(object):
    def __init__(self, pdf_A, pdf_B, rescale=None):
        self.pdf_A = pdf_A
        self.pdf_B = pdf_B
        self.fine_dr = self.pdf_A.dr/2.0
        self.rescale = rescale
        # initialise with large number
        self.similarity_distance = 1e10
        self.pdf_overlap()

    def pdf_overlap(self):
        """ Calculate the overlap of two PDFs via
        a simple meshed sum of their difference.
        """
        self.fine_space = np.arange(0, self.pdf_A.rmax, self.fine_dr)
        self.fine_Gr_A = np.interp(self.fine_space, self.pdf_A.r_space, self.pdf_A.Gr)
        self.fine_Gr_B = np.interp(self.fine_space, self.pdf_B.r_space, self.pdf_B.Gr)
        # discard last quarter before rmax to remove noise
        if self.rescale is not None:
            val = 0.2
            # scaling factor here is essentially equalising the shortest bond length
            bond_rescaling_factor = np.argmax(self.fine_Gr_B > val) / np.argmax(self.fine_Gr_A > val)
            # scaling factor here is normalising to number density
            density_rescaling_factor = pow((self.pdf_B.volume / self.pdf_B.num_atoms) / (self.pdf_A.volume / self.pdf_A.num_atoms), 1/3)
            if self.rescale == 'density':
                rescale_factor = density_rescaling_factor
            else:
                if self.rescale != 'bond':
                    print('Rescaling mode not specified, performing default (bond rescaling).')
                rescale_factor = bond_rescaling_factor
            self.fine_Gr_A = np.interp(self.fine_space, rescale_factor*self.fine_space, self.fine_Gr_A)
        self.fine_Gr_A = self.fine_Gr_A[:int(len(self.fine_space)*0.75)]
        self.fine_Gr_B = self.fine_Gr_B[:int(len(self.fine_space)*0.75)]
        self.fine_space = self.fine_space[:int(len(self.fine_space)*0.75)]
        self.overlap_fn = self.fine_Gr_A - self.fine_Gr_B
        self.worst_case_overlap_int = np.trapz(np.abs(self.fine_Gr_A), dx=self.pdf_A.dr/2.0) + \
            np.trapz(np.abs(self.fine_Gr_B), dx=self.pdf_B.dr/2.0)
        self.overlap_int = np.trapz(np.abs(self.overlap_fn), dx=self.pdf_A.dr/2.0)
        self.similarity_distance = self.overlap_int / self.worst_case_overlap_int

    def pdf_convolve(self, mode='same'):
        """ Calculate the convolution of two PDFs.
        """
        self.convolution = np.convolve(self.fine_Gr_A, self.fine_Gr_B, mode=mode)

    def plot_diff(self):
        """ Simple plot for comparing two PDF's. """
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(self.fine_space, self.fine_Gr_A, label=self.pdf_A.label)
        ax1.plot(self.fine_space, self.fine_Gr_B, label=self.pdf_B.label)
        ax1.legend(loc=1)
        ax1.set_xlabel('$r$ (Angstrom)')
        ax1.set_ylabel('$g(r)$')
        ax2.axhline(0, ls='--', c='k', lw=0.5)
        ax2.plot(self.fine_space, self.overlap_fn, ls='-')
        ax2.set_ylim(-0.5*ax1.get_ylim()[1], 0.5*ax1.get_ylim()[1])
        ax2.set_xlabel('$r$ (Angstrom)')
        ax2.set_ylabel('$g(r)$')
        return

    def plot_convolution(self):
        """ Plot the convolution of two PDFs. """
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(211)
        ax1.plot(np.arange(len(self.convolution), 0, step=-1) * self.fine_dr / 2.0,
                 self.convolution)
        ax1.set_ylabel('$g_A(r) \\ast g_B(r)$')
        ax1.set_xlabel('$\\Delta$ (Angstrom)')
        return
