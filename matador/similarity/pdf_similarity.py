# coding: utf-8
""" This file defines PDF and PDFOverlap as ways of calculating
the similarity between two structures.
"""

from itertools import product, combinations_with_replacement
from math import ceil
from copy import deepcopy
import time

import numpy as np
from scipy.spatial.distance import cdist

from matador.utils.cell_utils import frac2cart, cart2abc, cart2volume
from matador.utils.cell_utils import standardize_doc_cell
from matador.utils.print_utils import print_notify
from matador.similarity.fingerprint import Fingerprint


class PDF(Fingerprint):
    """ This class implements the calculation and comparison of pair
    distribution functions.

    Attributes:
        r_space (ndarray) : 1-D array containing real space grid
        gr (ndarray): 1-D array containing total PDF
        dr (float): real-space grid spacing in Å
        rmax (float): extent of real-space grid in Å
        label (str): structure label
        elem_gr (dict): dict with pairs of element symbol keys, containing 1-D arrays of projected PDFs (if calculated)
        number_density (float): number density for renormalisation and comparison with other PDFs
        kwargs (dict): arguments used to create PDF

    """

    def __init__(self, doc, **kwargs):
        """ Initialise parameters and run PDF (unless lazy=True).

        Parameters:

            doc (dict) : matador document to calculate PDF of

        Keyword Arguments:

            dr (float) : bin width for PDF (Angstrom) (DEFAULT: 0.01)
            gaussian_width (float) : width of Gaussian smearing (Angstrom) (DEFAULT: 0.01)
            num_images (int/str) : number of unit cell images include in PDF calculation (DEFAULT: 'auto')
            max_num_images (int) : cutoff number of unit cells before crashing (DEFAULT: 50)
            rmax (float) : maximum distance cutoff for PDF (Angstrom) (DEFAULT: 15)
            projected (bool) : optionally calculate the element-projected PDF
            standardize (bool) : optionally standardize cell before calculating PDF
            lazy (bool) : if True, calculator is not called when initializing PDF object
            timing (bool) : if True, print the total time taken to calculate the PDF

        """

        prop_defaults = {'dr': 0.01, 'gaussian_width': 0.01, 'rmax': 15, 'num_images': 'auto',
                         'style': 'smear', 'debug': False, 'timing': False, 'low_mem': False, 'projected': True,
                         'max_num_images': 50, 'standardize': False}

        # read and store kwargs
        self.kwargs = prop_defaults
        if 'sim_calc_args' in kwargs:
            self.kwargs.update(kwargs['sim_calc_args'])
        self.kwargs.update(kwargs)

        # standardize cell
        structure = deepcopy(doc)
        if self.kwargs.get('standardize'):
            structure = standardize_doc_cell(structure)

        # private variables
        self._num_images = self.kwargs.get('num_images')
        self._lattice = np.asarray(structure['lattice_cart'])
        self._poscart = np.asarray(frac2cart(structure['lattice_cart'], structure['positions_frac'])).reshape(-1, 3)
        self._types = structure['atom_types']
        self._num_atoms = len(self._poscart)
        self._volume = cart2volume(self._lattice)
        self._image_vec = None

        # public variables
        self.rmax = self.kwargs.get('rmax')
        self.number_density = self._num_atoms / self._volume
        self.dr = self.kwargs.get('dr')
        self.r_space = None
        self.gr = None
        self.elem_gr = None

        if 'text_id' in structure:
            self.label = ' '.join(structure['text_id'])
        else:
            self.label = 'null'

        if not kwargs.get('lazy'):
            if self.kwargs.get('timing'):
                start = time.time()
            self.calc_pdf()
            if self.kwargs.get('timing'):
                end = time.time()
                print('PDF calculated in {:.3f} s'.format(end - start))

    def calc_pdf(self):
        """ Wrapper to calculate PDF with current settings. """
        if self._image_vec is None:
            self._set_image_trans_vectors()
        if self.kwargs.get('projected'):
            if self.elem_gr is None:
                self._calc_projected_pdf()
        if self.gr is None:
            self._calc_unprojected_pdf()

    def _calc_distances(self, poscart, poscart_b=None):
        """ Calculate PBC distances with cdist.

        Parameters:
            poscart (ndarray): array of absolute atomic coordinates.

        Keyword Arguments:
            poscart_b (ndarray): absolute positions of a second type of atoms,
                                  where only A-B distances will be calculated.

        Returns:
            distances (ndarray): pair d_ij matrix with values > rmax < 1e-12 removed.


        """
        if self.kwargs.get('debug'):
            start = time.time()
        distances = np.array([])
        if poscart_b is None:
            poscart_b = deepcopy(poscart)
        for prod in self._image_vec:
            trans = np.zeros((3))
            for ind, multi in enumerate(prod):
                trans += self._lattice[ind] * multi
            distances = np.append(distances, cdist(poscart + trans, poscart_b))
        # mask by rmax/0 and remove masked values
        distances = np.ma.masked_where(distances > self.rmax, distances)
        distances = np.ma.masked_where(distances < 1e-12, distances)
        if self.kwargs.get('debug'):
            print('Calculated: {}, Used: {}, Ignored: {}'.format(len(distances),
                                                                 np.ma.count(distances),
                                                                 np.ma.count_masked(distances)))
        distances = distances.compressed()

        if self.kwargs.get('debug'):
            end = time.time()
            print('Calculated distances in {} s'.format(end - start))

        return distances

    def _calc_unprojected_pdf(self):
        """ Wrapper function to calculate distances and output
        a broadened and normalised PDF. Sets self.gr and self.r_space
        to G(r) and r respectively.

        """
        if self.elem_gr is not None:
            self._calc_unprojected_pdf_from_projected()
        else:
            distances = self._calc_distances(self._poscart)
            self.r_space = np.arange(0, self.rmax + self.dr, self.dr)
            self.gr = self._set_broadened_normalised_pdf(distances,
                                                         style=self.kwargs.get('style'),
                                                         gaussian_width=self.kwargs.get('gaussian_width'))

    def _calc_unprojected_pdf_from_projected(self):
        """" Reconstruct full PDF from projected. """
        self.gr = np.zeros_like(self.r_space)
        for key in self.elem_gr:
            self.gr += self.elem_gr[key]

    def _set_broadened_normalised_pdf(self, distances, style='smear', gaussian_width=0.01):
        """ Broaden the values provided as distances and return
        G(r) and r_space of the normalised PDF.

        Parameters:
            distances (numpy.ndarray): distances used to calculate PDF

        Keyword arguments:
            style (str): either 'smear' or 'histogram'
            gaussian_width (float): smearing width in Angstrom^1/2

        Returns:
            gr (np.ndarray): G(r), the PDF of supplied distances

        """

        if self.kwargs.get('debug'):
            start = time.time()

        hist = np.zeros_like(self.r_space, dtype=int)
        gr = np.zeros_like(self.r_space)
        for d_ij in distances:
            hist[ceil(d_ij / self.dr)] += 1
        if style == 'histogram' or gaussian_width == 0:
            # if hist, normalise G(r) by ideal gas then be done
            norm = 4 * np.pi * (self.r_space + self.dr)**2 * self.dr * self._num_atoms * self.number_density
            gr = np.divide(hist, norm)
        # otherwise do normal smearing
        else:
            if self.kwargs.get('low_mem'):
                if self.kwargs.get('debug'):
                    print('Using low memory mode...')
                new_space = np.reshape(self.r_space,
                                       (1, len(self.r_space))) - np.reshape(distances, (1, len(distances))).T
                gr = np.sum(np.exp(-(new_space)**2 / gaussian_width), axis=0)
            else:
                try:
                    new_space = np.reshape(self.r_space,
                                           (1, len(self.r_space))) - np.reshape(self.r_space, (1, len(self.r_space))).T
                    gr = np.sum(hist * np.exp(-(new_space)**2 / gaussian_width), axis=1)
                except MemoryError:
                    # if run out of memory, use low memory mode
                    if self.kwargs.get('debug'):
                        print('Ran out of memory, using low memory mode...')
                    self.kwargs['low_mem'] = True
                    new_space = np.reshape(self.r_space,
                                           (1, len(self.r_space))) - np.reshape(distances, (1, len(distances))).T
                    gr = np.sum(np.exp(-(new_space)**2 / gaussian_width), axis=0)

            # normalise G(r) by Gaussian integral and then ideal gas
            gr = np.divide(gr,
                           np.sqrt(np.pi * gaussian_width) *
                           4*np.pi * (self.r_space + self.dr)**2 * self._num_atoms * self.number_density)

        if self.kwargs.get('debug'):
            end = time.time()
            print('Calculated broadening and normalised in {} s'.format(end - start))

        return gr

    def _calc_projected_pdf(self):
        """ Calculate broadened and normalised element-projected PDF of a matador document.
        Sets self.elem_gr of e.g. Li2Zn3 to

            {
                ('Li', 'Li'): G_{Li-Li}(r),
                ('Li', 'Zn'): G_{Li-Zn}(r),
                ('Zn', 'Zn'): G_{Zn-Zn}(r)
            }


        """
        # initalise dict of element pairs with correct keys
        self.r_space = np.arange(0, self.rmax + self.dr, self.dr)
        elem_gr = dict()
        for comb in combinations_with_replacement(set(self._types), 2):
            elem_gr[tuple(set(comb))] = np.zeros_like(self.r_space)

        for elem_type in elem_gr:
            poscart = [self._poscart[i] for i in range(len(self._poscart)) if self._types[i] == elem_type[0]]
            poscart_b = ([self._poscart[i] for i in range(len(self._poscart)) if self._types[i] == elem_type[1]]
                         if len(elem_type) == 2 else None)
            distances = self._calc_distances(poscart, poscart_b=poscart_b)
            style = self.kwargs.get('style')
            gw = self.kwargs.get('gaussian_width')
            elem_gr[elem_type] = (
                len(elem_type) * self._set_broadened_normalised_pdf(distances, style=style, gaussian_width=gw)
            )

        self.elem_gr = elem_gr

    def _set_image_trans_vectors(self):
        """ Sets self._image_vec to a list/generator of image translation vectors,
        based on self._num_images.

        If self._num_images is an integer, create all 3-member integer combinations
        up to the value.

        If self._num_images is 'auto', create all translation vectors up to length self.rmax.

        e.g. self._image_vec = [[1, 0, 1], [0, 1, 1], [1, 1, 1]].

        """
        if self.kwargs.get('debug'):
            start = time.time()

        if self._num_images == 'auto':
            self._image_vec = set()
            any_in_sphere = True
            # find longest combination of single LV's
            max_trans = 0
            trans = np.zeros((3))
            for prod in product(range(-1, 2), repeat=3):
                trans = 0
                for ind, multi in enumerate(prod):
                    trans += self._lattice[ind] * multi
                if np.sqrt(np.sum(trans**2)) > max_trans:
                    max_trans = np.sqrt(np.sum(trans**2))
            first_attempt = 3
            test_num_images = deepcopy(first_attempt)
            while any_in_sphere:
                any_in_sphere = False
                for prod in product(range(-test_num_images, test_num_images + 1), repeat=3):
                    if prod in self._image_vec:
                        continue
                    trans = 0
                    for ind, multi in enumerate(prod):
                        trans += self._lattice[ind] * multi
                    if np.sqrt(np.sum(trans**2)) <= self.rmax + self.dr + max_trans:
                        self._image_vec.add(prod)
                        any_in_sphere = True
                test_num_images += 1
                if test_num_images > self.kwargs.get('max_num_images'):
                    print('Something has probably gone wrong; required images reached {}.'
                          .format(self.kwargs.get('max_num_images')))
                    print('lattice_abc:')
                    print('Continuing with num_images = 1')
                    self._num_images = 1
                    self._image_vec = list(product(range(-self._num_images, self._num_images + 1), repeat=3))
                    print(cart2abc(self._lattice))
                    break
        else:
            self._image_vec = list(product(range(-self._num_images, self._num_images + 1), repeat=3))
        if self.kwargs.get('debug'):
            end = time.time()
            print('Image vectors length = {}'.format(len(self._image_vec)))
            print('Set image trans vectors in {} s'.format(end - start))

    def get_sim_distance(self, pdf_b, projected=False):
        """ Return the similarity between two PDFs. """
        return PDFOverlap(self, pdf_b, projected=projected).similarity_distance

    def pdf(self):
        """ Return G(r) and the r_space for easy plotting. """
        try:
            return (self.r_space, self.gr)
        except AttributeError:
            return (None, None)

    def plot_projected_pdf(self, keys=None, other_pdfs=None, cmap='Dark2'):
        """ Plot projected PDFs.

        Keyword arguments:
            keys (list): plot only a subset of projections, e.g. [('K', )].
            other_pdfs (list of PDF): other PDFs to plot.
            cmap (str): name of matplotlib colourmap.

        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set(style='whitegrid', font_scale=1.2)
        fig = plt.figure(figsize=(8, 5))
        ax1 = fig.add_subplot(111)
        ax1.plot(self.r_space, self.gr, lw=1, zorder=100000, ls='-', label='total {}'.format(self.label), c='k')
        if keys is None:
            keys = [key for key in self.elem_gr]
        for key in keys:
            ax1.plot(self.r_space, self.elem_gr[key], label='-'.join(key) + ' {}'.format(self.label))
        if other_pdfs is not None:
            if isinstance(other_pdfs, PDF):
                other_pdfs = [other_pdfs]
            for pdf in other_pdfs:
                if isinstance(pdf, PDF):
                    ax1.plot(pdf.r_space, pdf.gr, lw=1, zorder=99999, ls='--',
                             label='total {}'.format(pdf.label), c='k')
                    for key in keys:
                        ax1.plot(self.r_space, pdf.elem_gr[key], ls='--',
                                 label='-'.join(key) + ' {}'.format(pdf.label))
                elif isinstance(pdf, tuple):
                    ax1.plot(pdf[0], pdf[1], lw=2, alpha=1, ls='--')
                else:
                    raise RuntimeError
        ax1.legend(loc=1)
        ax1.set_ylabel('$g(r)$')
        ax1.set_xlabel('$r$ (Angstrom)')
        plt.show()

    def plot_pdf(self, other_pdfs=None, cmap='Dark2'):
        """ Plot PDFs.

        Keyword arguments:
            other_pdfs (list of PDF): other PDFs to add to the plot.
            cmap (str): name of matplotlib colourmap.


        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set(style='whitegrid', font_scale=1.2)
        fig = plt.figure(figsize=(8, 5))
        ax1 = fig.add_subplot(111)
        ax1.plot(self.r_space, self.gr, lw=2, label=self.label, c='k')
        ax1.set_ylabel('Pair distribution function, $g(r)$')
        ax1.set_xlim(0, self.rmax)
        if other_pdfs is not None:
            if isinstance(other_pdfs, PDF):
                other_pdfs = [other_pdfs]
            for pdf in other_pdfs:
                if isinstance(pdf, PDF):
                    ax1.plot(pdf.r_space, pdf.gr, lw=2, label=pdf.label, alpha=1, ls='--')
                elif isinstance(pdf, tuple):
                    ax1.plot(pdf[0], pdf[1], lw=2, alpha=1, ls='--')
                else:
                    raise RuntimeError
        ax1.set_xlabel('$r$ (Angstrom)')
        ax1.legend()
        plt.show()


class PDFFactory:
    """ This class computes PDF objects from a list of structures,
    as concurrently as possible. The PDFs are stored under the `pdf`
    key inside each structure dict.

    Attributes:
        nprocs (int): number of concurrent processes.

    """

    def __init__(self, cursor, debug=False, concurrency='pool', **pdf_args):
        """ Compute PDFs over n processes, where n is set by either
        SLURM_NTASKS, OMP_NUM_THREADS or physical core count.

        Parameters:
            cursor (list of dict): list of matador structures

        Keyword arguments:
            concurrency (str): either 'pool' or 'queue'
            pdf_args (dict): arguments to pass to the PDF calculator

        """
        import multiprocessing as mp
        import os
        # create list of empty (lazy) PDF objects
        if 'lazy' in pdf_args:
            del pdf_args['lazy']
        for doc in cursor:
            doc['pdf'] = PDF(doc, lazy=True, **pdf_args)

        # how many processes to use? either SLURM_NTASKS, OMP_NUM_THREADS or total num CPUs
        if os.environ.get('SLURM_NTASKS') is not None:
            self.nprocs = int(os.environ.get('SLURM_NTASKS'))
            env = '$SLURM_NTASKS'
        elif os.environ.get('OMP_NUM_THREADS') is not None:
            self.nprocs = int(os.environ.get('OMP_NUM_THREADS'))
            env = '$OMP_NUM_THREADS'
        else:
            self.nprocs = mp.cpu_count()
            env = 'core count'
        print_notify('Running {} jobs on {} processes in {} mode, set by {}.'.format(len(cursor),
                                                                                     self.nprocs,
                                                                                     concurrency,
                                                                                     env))

        if debug:
            print('Initialising worker {}...'.format(concurrency))

        start = time.time()
        if concurrency == 'queue':
            queue = mp.Queue()
            # split cursor into subcursors
            pdfs_per_proc = int(round((len(cursor) / self.nprocs)))
            subcursors = []
            for i in range(self.nprocs):
                if i == self.nprocs - 1:
                    subcursors.append(deepcopy(cursor[i * pdfs_per_proc:]))
                else:
                    subcursors.append(deepcopy(cursor[i * pdfs_per_proc:(i + 1) * pdfs_per_proc]))

            processes = [mp.Process(target=calc_pdf_queue_wrapper,
                                    args=(subcursors[i], i, queue))
                         for i in range(self.nprocs)]
            for proc in processes:
                proc.start()

            results_cursor = dict()
            for i in range(self.nprocs):
                results_cursor.update(queue.get())

            assert len(results_cursor) == self.nprocs
            pdf_cursor = []
            for i in range(self.nprocs):
                pdf_cursor.append(results_cursor[i])
            pdf_cursor = [doc for subcursor in pdf_cursor for doc in subcursor]

        elif concurrency == 'pool':
            pool = mp.Pool(processes=self.nprocs)
            pdf_cursor = []
            pool.map_async(calc_pdf_pool_wrapper, cursor, callback=pdf_cursor.extend, error_callback=print)
            pool.close()
            pool.join()

        if len(pdf_cursor) != len(cursor):
            raise RuntimeError('There was an error calculating the desired PDFs')

        for ind, doc in enumerate(cursor):
            cursor[ind]['pdf'] = pdf_cursor[ind]['pdf']

        elapsed = time.time() - start
        if debug:
            print('Compute time: {:.4f} s'.format(elapsed))
            print('Work complete!')


def calc_pdf_queue_wrapper(cursor, i, queue):
    """ Evaluate PDFs of a cursor where a lazy init of each doc's PDF object has
    already been made. The result is parcelled into a dictionary with key i
    and pushed to the queue.

    Parameters:
        cursor (list of dict): list of matador structures with empty PDF
            objects stored under `pdf`.
        i (int): position of cursor in overall subcursor array
        queue (multiprocessing.Queue): processing queue.

    """
    for ind, _ in enumerate(cursor):
        cursor[ind]['pdf'].calc_pdf()
    queue.put({i: cursor})


def calc_pdf_pool_wrapper(doc):
    """ Evaluate PDF of a structure where a lazy init of the doc's PDF object has
    already been made.

    Parameters:
        doc (dict): matador structures with empty PDF

    """
    doc['pdf'].calc_pdf()
    return {'pdf': doc['pdf']}


class PDFOverlap:
    """ Calculate the PDFOverlap between two PDF objects,
    pdf_a and pdf_b, with number density rescaling.

    Attributes:
        pdf_a (PDF): first PDF to compare.
        pdf_b (PDF): second PDF to compare.
        fine_dr (float): fine grid scale on which to compare.
        similarity_distance (float): "distance" between PDFs.
        overlap_int (float): the value of the overlap integral.

    """

    def __init__(self, pdf_a, pdf_b, projected=False):
        """ Perform the overlap and similarity distance calculations.

        Parameters:
            pdf_a (PDF): first PDF to compare.
            pdf_b (PDF): second PDF to compare.

        Keyword arguments:
            projected : if True, attempt to use projected PDFs.

        """
        self.pdf_a = pdf_a
        self.pdf_b = pdf_b
        self.fine_dr = self.pdf_a.dr / 2.0
        # initialise with large number
        self.similarity_distance = 1e10
        self.overlap_int = 0
        if projected:
            if isinstance(pdf_a.elem_gr, dict) and isinstance(pdf_b.elem_gr, dict):
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
        self.fine_space = np.arange(0, self.pdf_a.rmax, self.fine_dr)
        self.fine_gr_a = np.interp(self.fine_space, self.pdf_a.r_space, self.pdf_a.gr)
        self.fine_gr_b = np.interp(self.fine_space, self.pdf_b.r_space, self.pdf_b.gr)
        # scaling factor here is normalising to number density
        density_rescaling_factor = pow(self.pdf_b.number_density / (self.pdf_a.number_density), 1 / 3)
        rescale_factor = density_rescaling_factor
        self.fine_gr_a = np.interp(self.fine_space, rescale_factor * self.fine_space, self.fine_gr_a)
        self.fine_gr_a = self.fine_gr_a[:int(len(self.fine_space) * 0.75)]
        self.fine_gr_b = self.fine_gr_b[:int(len(self.fine_space) * 0.75)]
        self.fine_space = self.fine_space[:int(len(self.fine_space) * 0.75)]
        overlap_fn = self.fine_gr_a - self.fine_gr_b
        worst_case_overlap_int = np.trapz(np.abs(self.fine_gr_a), dx=self.pdf_a.dr/2.0) + \
            np.trapz(np.abs(self.fine_gr_b), dx=self.pdf_b.dr/2.0)
        self.overlap_int = np.trapz(np.abs(overlap_fn), dx=self.pdf_a.dr / 2.0)
        self.similarity_distance = self.overlap_int / worst_case_overlap_int

    def projected_pdf_overlap(self):
        """ Calculate the overlap of two projected PDFs via
        a simple meshed sum of their difference.

        """
        self.fine_space = np.arange(0, self.pdf_a.rmax, self.fine_dr)
        self.overlap_int = 0
        self.similarity_distance = 1e10
        elems = set(key for key in self.pdf_a.elem_gr)
        if elems != set(key for key in self.pdf_b.elem_gr):
            for key in self.pdf_b.elem_gr:
                elems.add(key)
        # pad out missing elements with zero PDFs
        for key in elems:
            if key not in self.pdf_a.elem_gr:
                self.pdf_a.elem_gr[key] = np.zeros_like(self.pdf_a.r_space)
            if key not in self.pdf_b.elem_gr:
                self.pdf_b.elem_gr[key] = np.zeros_like(self.pdf_b.r_space)
        self.fine_elem_gr_a, self.fine_elem_gr_b = dict(), dict()
        for key in elems:
            self.fine_elem_gr_a[key] = np.interp(self.fine_space, self.pdf_a.r_space, self.pdf_a.elem_gr[key])
            self.fine_elem_gr_b[key] = np.interp(self.fine_space, self.pdf_b.r_space, self.pdf_b.elem_gr[key])
        # scaling factor here is normalising to number density
        density_rescaling_factor = pow((self.pdf_b.number_density) / (self.pdf_a.number_density), 1 / 3)
        rescale_factor = density_rescaling_factor
        for key in elems:
            self.fine_elem_gr_a[key] = (
                np.interp(self.fine_space, rescale_factor * self.fine_space, self.fine_elem_gr_a[key])
            )
        for key in elems:
            self.fine_elem_gr_a[key] = self.fine_elem_gr_a[key][:int(len(self.fine_space) * 0.75)]
            self.fine_elem_gr_b[key] = self.fine_elem_gr_b[key][:int(len(self.fine_space) * 0.75)]
        self.fine_space = self.fine_space[:int(len(self.fine_space) * 0.75)]

        for key in elems:
            overlap_fn = self.fine_elem_gr_a[key] - self.fine_elem_gr_b[key]
            worst_case_a = np.trapz(np.abs(self.fine_elem_gr_a[key]), dx=self.pdf_a.dr / 2.0)
            worst_case_b = np.trapz(np.abs(self.fine_elem_gr_b[key]), dx=self.pdf_b.dr / 2.0)
            worst_case_overlap = worst_case_a + worst_case_b
            overlap = np.trapz(np.abs(overlap_fn), dx=self.pdf_a.dr / 2.0)
            self.overlap_int += overlap / worst_case_overlap

        self.similarity_distance = self.overlap_int / len(elems)

    def plot_diff(self, cmap='Dark2'):
        """ Simple plot for comparing two PDF's. """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import seaborn as sns
        sns.set(style='whitegrid', font_scale=1.2)

        plt.figure(figsize=(8, 6))
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        gs.update(hspace=0)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1)

        ax2.set_xlabel('$r$ (\\AA)')
        ax1.set_ylabel('$g(r)$')
        ax2.set_ylabel('$g_a(r) - g_b(r)$')
        ax2.axhline(0, ls='--', c='k', lw=0.5)
        ax1.set_xlim(0, np.max(self.fine_space))

        ax1.plot(self.fine_space, self.fine_gr_a, label=self.pdf_a.label, c=colours[0])
        ax1.plot(self.fine_space, self.fine_gr_b, label=self.pdf_b.label, c=colours[1])

        plt.setp(ax1.get_xticklabels(), visible=False)
        ax2.set_ylim(-0.5 * ax1.get_ylim()[1], 0.5 * ax1.get_ylim()[1])

        ax1.legend(loc=0)
        ax2.plot(self.fine_space, self.overlap_fn, ls='-', c=colours[2])
        ax2.set_ylim(-0.5 * ax1.get_ylim()[1], 0.5 * ax1.get_ylim()[1])
        plt.tight_layout()
        plt.show()

    def plot_projected_diff(self, cmap='Dark2'):
        """ Simple plot for comparing two PDF's. """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import seaborn as sns
        sns.set(style='whitegrid', font_scale=1.2)

        plt.figure(figsize=(8, 6))
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        gs.update(hspace=0)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1], sharex=ax1)
        ax2.set_xlabel('$r$ (\\AA)')
        ax1.set_ylabel('$g(r)$')
        ax2.set_ylabel('$g_a(r) - g_b(r)$')
        ax2.axhline(0, ls='--', c='k', lw=0.5)
        ax1.set_xlim(0, np.max(self.fine_space))
        for ind, key in enumerate(self.fine_elem_gr_a):
            ax1.plot(self.fine_space, self.fine_elem_gr_a[key],
                     label='-'.join(key) + ' {}'.format(self.pdf_a.label),
                     c=colours[ind])
            ax1.plot(self.fine_space, self.fine_elem_gr_b[key],
                     label='-'.join(key) + ' {}'.format(self.pdf_b.label),
                     ls='--', c=colours[ind])
            ax2.plot(self.fine_space, self.fine_elem_gr_a[key] - self.fine_elem_gr_b[key],
                     label='-'.join(key) + ' diff',
                     c=colours[ind], ls='-')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax2.set_ylim(-0.5 * ax1.get_ylim()[1], 0.5 * ax1.get_ylim()[1])
        ax1.legend(loc=0)
        ax2.legend(loc=2)
        plt.tight_layout()
        plt.show()
