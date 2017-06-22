# coding: utf-8
""" This file defines various measures and ways of calculating
the similarity between two structures.

TO-DO:
    * implement more similarity distances (a la Oganov)
    * otf calculation of required num_images
    * non-diagonal supercells
    * element-projected Fortran PDF calculator

"""

# matador modules
from matador.utils.cell_utils import frac2cart, cart2abc
# external libraries
import numpy as np
from scipy.spatial.distance import cdist
# standard library
from itertools import product, combinations_with_replacement
from math import ceil
from copy import deepcopy


class PDF(object):
    """ This class implements the calculation and comparison of pair
    distribution functions.

    Args:
        doc            : matador document to calculate PDF of.
        dr             : bin width for PDF (Angstrom) (DEFAULT: 0.1)
        gaussian_width : width of Gaussian smearing (Angstrom) (DEFAULT: 0.01)
        num_images     : number of unit cell images include in PDF calculation (DEFAULT: 'auto')
        rmax           : maximum distance cutoff for PDF (Angstrom) (DEFAULT: 15)
        calculator     : F or None, for Fortran or Python calculator (DEFAULT: None)
        projected      : True/False, optionally calculate the element-projected PDF.

    """
    def __init__(self, doc, **kwargs):
        """ Initialise parameters. """
        if 'sim_calc_args' in kwargs:
            kwargs = kwargs['sim_calc_args']
        if kwargs.get('dr') is None:
            self.dr = 0.1
        else:
            self.dr = kwargs['dr']
        if kwargs.get('gaussian_width') is None:
            self.gaussian_width = 0.01
        else:
            self.gaussian_width = kwargs['gaussian_width']
        if kwargs.get('rmax') is None:
            self.rmax = 15
        else:
            self.rmax = kwargs['rmax']
        if kwargs.get('num_images') is None:
            self.num_images = 'auto'
        else:
            self.num_images = kwargs['num_images']
        if kwargs.get('calculator') is None:
            self._calc_pdf = self._calc_py_pdf
        else:
            self._calc_pdf = self._calc_fortran_pdf
        if kwargs.get('style') is not None:
            self.style = kwargs.get('style')
        else:
            self.style = 'smear'
        if kwargs.get('debug'):
            self.debug = True
        else:
            self.debug = False
        if kwargs.get('low_mem'):
            self.low_mem = True
        else:
            self.low_mem = False
        if kwargs.get('max_num_images'):
            self.max_num_images = kwargs.get('max_num_images')
        else:
            self.max_num_images = 1000
        self.doc = doc
        self.lattice = np.asarray(doc['lattice_cart'])
        self.poscart = np.asarray(frac2cart(doc['lattice_cart'], doc['positions_frac']))
        self.types = doc['atom_types']
        if 'text_id' in doc:
            self.label = ' '.join(doc['text_id'])
        else:
            self.label = 'null'
        self.num_atoms = len(self.poscart)
        self.volume = doc['cell_volume']
        self.number_density = self.num_atoms / self.volume
        self._set_image_trans_vectors()
        if not kwargs.get('lazy'):
            self._calc_pdf()
            if kwargs.get('projected'):
                self._calc_py_projected_pdf()

    def _calc_fortran_pdf(self):
        """ Calculate PDF of a matador document with Fortran calculator. """
        # from similarity.pdf.pdf_calculator import pdf_calc
        raise NotImplementedError

    def _calc_distances(self, poscart, poscart_B=None, debug=False):
        """ Calculate PBC distances with cdist.

        Input:

            poscart    : np.array of Cartesian atomic coordinates.
            poscart_B  : (OPTIONAL) positions of a second type of atoms,
                         where only A-B distances will be calculated.

        Returns:

            Sets and returns self.distances to d_ij matrix,
            with values > rmax < 1e-12 removed.


        """
        distances = np.array([])
        if poscart_B is None:
            poscart_B = deepcopy(poscart)
        for prod in self.image_vec:
            trans = np.zeros((3))
            for ind, multi in enumerate(prod):
                trans += self.lattice[ind] * multi
            distances = np.append(distances, cdist(poscart+trans, poscart_B))
        # mask by rmax/0 and remove masked values
        distances = np.ma.masked_where(distances > self.rmax, distances)
        distances = np.ma.masked_where(distances < 1e-12, distances)
        if debug:
            print('Calculated: {}, Used: {}, Ignored: {}'.format(len(distances),
                                                                 np.ma.count(distances),
                                                                 np.ma.count_masked(distances)))
        self.distances = distances.compressed()
        return self.distances

    def _calc_py_pdf(self, debug=False):
        """ Wrapper function to calculate distances and output
        a broadened and normalised PDF.

        Returns:

            Sets self.Gr and self.r_space to G(r) and r respectively.

        """
        if self.debug:
            import time
            start = time.time()
        distances = self._calc_distances(self.poscart, debug=self.debug)
        if self.debug:
            end = time.time()
            print('Calculated distances in {} s'.format(end-start))
        if self.debug:
            start = time.time()
        self.r_space = np.arange(0, self.rmax+self.dr, self.dr)
        self.Gr = self._set_broadened_normalised_pdf(distances,
                                                     style=self.style,
                                                     gaussian_width=self.gaussian_width)
        if self.debug:
            end = time.time()
            print('Calculated broadening and normalised in {} s'.format(end-start))

    def _set_broadened_normalised_pdf(self, distances, style='smear', gaussian_width=0.01):
        """ Broaden the values provided as distances and return
        G(r) and r_space of the normalised PDF.

        Input:

            distances      : used to calculate PDF
            style          : either 'smear' or 'histogram'
            gaussian_width : smearing width in Angstrom^1/2

        Requires:

            self.r_space to be set.

        Returns and sets:

            Gr             : G(r), the PDF of supplied distances

        """
        hist = np.zeros_like(self.r_space, dtype=int)
        Gr = np.zeros_like(self.r_space)
        for d_ij in self.distances:
            hist[ceil(d_ij/self.dr)] += 1
        if style == 'histogram' or gaussian_width == 0:
            # if hist, normalise G(r) by ideal gas then be done
            Gr = np.divide(hist,
                           4*np.pi * (self.r_space + self.dr)**2 * self.dr * self.num_atoms * self.number_density)
        elif style == 'smear':
            if not self.low_mem:
                try:
                    # otherwise, stack some Gaussians on that PDF
                    new_space = np.reshape(self.r_space, (1, len(self.r_space))) - np.reshape(self.r_space, (1, len(self.r_space))).T
                    Gr = np.sum(hist*np.exp(-(new_space)**2 / gaussian_width), axis=1)
                except MemoryError:
                    print('Memory usage too high; consider decreasing dr or turning on low_mem mode.')
            else:
                if self.debug:
                    print('Using low memory mode...')
                # if run out of memory, use low memory mode
                new_space = np.reshape(self.r_space, (1, len(self.r_space))) - np.reshape(self.distances, (1, len(self.distances))).T
                Gr = np.sum(np.exp(-(new_space)**2 / gaussian_width), axis=0)
            # normalise G(r) by Gaussian integral and then ideal gas
            Gr = np.divide(Gr,
                           np.sqrt(np.pi * gaussian_width) *
                           4*np.pi * (self.r_space + self.dr)**2 * self.num_atoms * self.number_density)
        return Gr

    def _calc_py_projected_pdf(self):
        """ Calculate broadened and normalised element-projected PDF of a matador document.

        Sets self.elem_Gr of e.g. Li2Zn3 to

            {
                ('Li', 'Li'): G_{Li-Li}(r),
                ('Li', 'Zn'): G_{Li-Zn}(r),
                ('Zn', 'Zn'): G_{Zn-Zn}(r)
            }


        """
        # initalise dict of element pairs with correct keys
        self.r_space = np.arange(0, self.rmax+self.dr, self.dr)
        self.elem_Gr = dict()
        for comb in combinations_with_replacement(set(self.doc['atom_types']), 2):
            self.elem_Gr[tuple(set(comb))] = np.zeros_like(self.r_space)

        distances = dict()
        for elem_type in self.elem_Gr:
            poscart = [self.poscart[i] for i in range(len(self.poscart)) if self.doc['atom_types'][i] == elem_type[0]]
            poscart_B = ([self.poscart[i] for i in range(len(self.poscart)) if self.doc['atom_types'][i] == elem_type[1]]
                         if len(elem_type) == 2 else None)
            distances[elem_type] = self._calc_distances(poscart, poscart_B=poscart_B)

            self.elem_Gr[elem_type] = self._set_broadened_normalised_pdf(distances[elem_type],
                                                                         gaussian_width=self.gaussian_width)

    def _calc_py_deprecated_projected_pdf(self):
        """ DEPCRECATED: calculate element-projected PDF of a matador document with Python calculator. """
        raise DeprecationWarning
        self.elem_Gr = dict()
        for comb in combinations_with_replacement(set(self.doc['atom_types']), 2):
            self.elem_Gr[tuple(set(comb))] = np.zeros((int(self.rmax / self.dr)))
        self._deprecated_Gr = np.zeros((int(self.rmax / self.dr)))
        for i in range(self.num_atoms):
            for j in range(i+1, self.num_atoms):
                d_ij = np.sqrt(np.sum((self.poscart[i] - self.poscart[j])**2))
                if d_ij <= self.rmax:
                    self._deprecated_Gr += 2*np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / (self.num_atoms * (self.num_images+1)**3)
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
                    d_ij = np.sqrt(np.sum((self.poscart[i] - self.poscart[j] - trans)**2))
                    if d_ij <= self.rmax:
                        self._deprecated_Gr += np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / (self.num_atoms * (self.num_images+1)**3)
                        self.elem_Gr[tuple(set((self.types[i], self.types[j])))] += np.exp(-(self.r_space - d_ij)**2 / self.gaussian_width) / (self.num_atoms * (self.num_images+1)**3)
        return

    def _set_image_trans_vectors(self):
        """ Sets self.image_vec to a list/generator of image translation vectors,
        based on self.num_images.

        If self.num_images is an integer, create all 3-member integer combinations
        up to the value.

        If self.num_images is 'auto', create all translation vectors up to length self.rmax.

        e.g. self.image_vec = [[1, 0, 1], [0, 1, 1], [1, 1, 1]].

        """
        if self.num_images is 'auto':
            self.image_vec = set()
            any_in_sphere = True
            # find longest combination of single LV's
            max_trans = 0
            for prod in product(range(-1, 2), repeat=3):
                trans = np.zeros((3))
                for ind, multi in enumerate(prod):
                    trans += self.lattice[ind] * multi
                if np.sqrt(np.sum(trans**2)) > max_trans:
                    max_trans = np.sqrt(np.sum(trans**2))
            test_num_images = 1
            while any_in_sphere:
                any_in_sphere = False
                for prod in product(range(-test_num_images, test_num_images+1), repeat=3):
                    if prod in self.image_vec:
                        continue
                    trans = np.zeros((3))
                    for ind, multi in enumerate(prod):
                        trans += self.lattice[ind] * multi
                    if np.sqrt(np.sum(trans**2)) <= self.rmax+self.dr+0.5*max_trans:
                        self.image_vec.add(prod)
                        any_in_sphere = True
                test_num_images += 1
                if test_num_images > self.max_num_images:
                    print(self.image_vec)
                    print('Something has probably gone wrong; required images reached {}.'.format(self.max_num_images))
                    print(self.doc['text_id'])
                    if 'lattice_abc' in self.doc:
                        print(self.doc['lattice_abc'])
                    else:
                        print(cart2abc(self.doc['lattice_cart']))
                    raise RuntimeError
        else:
            self.image_vec = product(range(-self.num_images, self.num_images+1), repeat=3)

    def get_sim_distance(self, pdf_B, projected=False):
        """ Return the similarity between two PDFs. """
        return PDFOverlap(self, pdf_B, projected=projected).similarity_distance

    def pdf(self):
        """ Return G(r) and the r_space for easy plotting. """
        try:
            return (self.r_space, self.Gr)
        except:
            return (None, None)

    def plot_projected_pdf(self, keys=None):
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

    def plot_pdf(self, other_pdfs=None):
        """ Plot projected PDFs, with optional list of
        tuples [(r_space, Gr), ...] of other PDFs.
        """
        import matplotlib.pyplot as plt
        import seaborn
        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(111)
        ax1.plot(self.r_space, self.Gr, lw=2, label=self.label)
        if other_pdfs is not None:
            for pdf in other_pdfs:
                if isinstance(pdf, PDF):
                    ax1.plot(pdf.r_space, pdf.Gr, lw=2, label=pdf.label)
                elif isinstance(pdf, tuple):
                    ax1.plot(pdf[0], pdf[1], lw=2)
                else:
                    raise RuntimeError
        ax1.set_ylabel('$g(r)$')
        ax1.set_xlabel('$r$ (Angstrom)')
        plt.legend()
        plt.show()
        return


class PDFOverlap(object):
    """ Calculate the PDFOverlap between two PDF objects,
    pdf_A and pdf_B, with number density rescaling.

    Args:
        pdf_A/B : two PDF objects to compare.
        projected : if True, attempt to use projected PDFs.


    """
    def __init__(self, pdf_A, pdf_B, projected=False):
        """ Perform the overlap and similarity distance calculations. """
        self.pdf_A = pdf_A
        self.pdf_B = pdf_B
        self.fine_dr = self.pdf_A.dr/2.0
        # initialise with large number
        self.similarity_distance = 1e10
        self.overlap_int = 0
        if projected:
            if isinstance(pdf_A.elem_Gr, dict) and isinstance(pdf_B.elem_Gr, dict):
                self.projected_pdf_overlap()
            else:
                print('Projected PDFs missing, continuing with total.')
                self.pdf_overlap()
        else:
            self.pdf_overlap()

    def pdf_overlap(self):
        """ Calculate the overlap of two PDFs via
        a simple meshed sum of their difference.
        """
        self.overlap_int = 0
        self.similarity_distance = 1e10
        self.fine_space = np.arange(0, self.pdf_A.rmax, self.fine_dr)
        self.fine_Gr_A = np.interp(self.fine_space, self.pdf_A.r_space, self.pdf_A.Gr)
        self.fine_Gr_B = np.interp(self.fine_space, self.pdf_B.r_space, self.pdf_B.Gr)
        # scaling factor here is normalising to number density
        density_rescaling_factor = pow((self.pdf_B.volume / self.pdf_B.num_atoms) / (self.pdf_A.volume / self.pdf_A.num_atoms), 1/3)
        rescale_factor = density_rescaling_factor
        self.fine_Gr_A = np.interp(self.fine_space, rescale_factor*self.fine_space, self.fine_Gr_A)
        self.fine_Gr_A = self.fine_Gr_A[:int(len(self.fine_space)*0.75)]
        self.fine_Gr_B = self.fine_Gr_B[:int(len(self.fine_space)*0.75)]
        self.fine_space = self.fine_space[:int(len(self.fine_space)*0.75)]
        self.overlap_fn = self.fine_Gr_A - self.fine_Gr_B
        self.worst_case_overlap_int = np.trapz(np.abs(self.fine_Gr_A), dx=self.pdf_A.dr/2.0) + \
            np.trapz(np.abs(self.fine_Gr_B), dx=self.pdf_B.dr/2.0)
        self.overlap_int = np.trapz(np.abs(self.overlap_fn), dx=self.pdf_A.dr/2.0)
        self.similarity_distance = self.overlap_int / self.worst_case_overlap_int

    def projected_pdf_overlap(self):
        """ Calculate the overlap of two projected PDFs via
        a simple meshed sum of their difference.
        """
        self.fine_space = np.arange(0, self.pdf_A.rmax, self.fine_dr)
        self.overlap_int = 0
        self.similarity_distance = 1e10
        elems = set(key for key in self.pdf_A.elem_Gr)
        if elems != set(key for key in self.pdf_B.elem_Gr):
            for key in self.pdf_B.elem_Gr:
                elems.add(key)
        # pad out missing elements with zero PDFs
        for key in elems:
            if key not in self.pdf_A.elem_Gr:
                self.pdf_A.elem_Gr[key] = np.zeros_like(self.pdf_A.r_space)
            if key not in self.pdf_B.elem_Gr:
                self.pdf_B.elem_Gr[key] = np.zeros_like(self.pdf_B.r_space)
        self.fine_elem_Gr_A, self.fine_elem_Gr_B = dict(), dict()
        for key in elems:
            self.fine_elem_Gr_A[key] = np.interp(self.fine_space, self.pdf_A.r_space, self.pdf_A.elem_Gr[key])
            self.fine_elem_Gr_B[key] = np.interp(self.fine_space, self.pdf_B.r_space, self.pdf_B.elem_Gr[key])
        # scaling factor here is normalising to number density
        density_rescaling_factor = pow((self.pdf_B.volume / self.pdf_B.num_atoms) / (self.pdf_A.volume / self.pdf_A.num_atoms), 1/3)
        rescale_factor = density_rescaling_factor
        for key in elems:
            self.fine_elem_Gr_A[key] = np.interp(self.fine_space, rescale_factor*self.fine_space, self.fine_elem_Gr_A[key])
        for key in elems:
            self.fine_elem_Gr_A[key] = self.fine_elem_Gr_A[key][:int(len(self.fine_space)*0.75)]
            self.fine_elem_Gr_B[key] = self.fine_elem_Gr_B[key][:int(len(self.fine_space)*0.75)]
        self.fine_space = self.fine_space[:int(len(self.fine_space)*0.75)]
        self.overlap_fn = dict()
        for key in elems:
            self.overlap_fn[key] = self.fine_elem_Gr_A[key] - self.fine_elem_Gr_B[key]
        self.worst_case_overlap_int = dict()
        for key in elems:
            self.worst_case_overlap_int[key] = np.trapz(np.abs(self.fine_elem_Gr_A[key]), dx=self.pdf_A.dr/2.0) + \
                np.trapz(np.abs(self.fine_elem_Gr_B[key]), dx=self.pdf_B.dr/2.0)
        for key in elems:
            self.overlap_int += np.trapz(np.abs(self.overlap_fn[key]), dx=self.pdf_A.dr/2.0) / self.worst_case_overlap_int[key]
        self.similarity_distance = self.overlap_int / len(elems)

    def pdf_convolve(self, mode='same'):
        """ Calculate the convolution of two PDFs.
        """
        self.convolution = np.convolve(self.fine_Gr_A, self.fine_Gr_B, mode=mode)

    def plot_diff(self):
        """ Simple plot for comparing two PDF's. """
        import matplotlib.pyplot as plt
        import seaborn
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
        plt.show()
        return

    def plot_convolution(self):
        """ Plot the convolution of two PDFs. """
        import matplotlib.pyplot as plt
        import seaborn
        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(211)
        ax1.plot(np.arange(len(self.convolution), 0, step=-1) * self.fine_dr / 2.0,
                 self.convolution)
        ax1.set_ylabel('$g_A(r) \\ast g_B(r)$')
        ax1.set_xlabel('$\\Delta$ (Angstrom)')
        plt.show()
        return
