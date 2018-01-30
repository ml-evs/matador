# coding: utf-8
""" This file defines various measures and ways of calculating
the similarity between two structures.
"""

# matador modules
from matador.utils.cell_utils import frac2cart, cart2abc, cart2volume
from matador.utils.cell_utils import standardize_doc_cell
from matador.utils.print_utils import print_notify

# external libraries
import numpy as np
from scipy.spatial.distance import cdist
# standard library
from itertools import product, combinations_with_replacement
from math import ceil
from copy import deepcopy
import time


class PDF:
    """ This class implements the calculation and comparison of pair
    distribution functions.

    Args:
        | doc            : dict, matador document to calculate PDF of
        | dr             : float, bin width for PDF (Angstrom) (DEFAULT: 0.1)
        | gaussian_width : float, width of Gaussian smearing (Angstrom) (DEFAULT: 0.01)
        | num_images     : int/str, number of unit cell images include in PDF calculation (DEFAULT: 'auto')
        | max_num_images : int, cutoff number of unit cells before crashing (DEFAULT: 50)
        | rmax           : float, maximum distance cutoff for PDF (Angstrom) (DEFAULT: 15)
        | projected      : bool, optionally calculate the element-projected PDF
        | standardize    : bool, optionally standardize cell before calculating PDF
        | lazy           : bool, if True, calculator is not called when initializing PDF object
        | timing         : bool, if True, print the total time taken to calculate the PDF

    """
    def __init__(self, doc, **kwargs):
        """ Initialise parameters and run PDF. """
        prop_defaults = {'dr': 0.1, 'gaussian_width': 0.01, 'rmax': 15, 'num_images': 'auto',
                         'style': 'smear', 'debug': False, 'timing': False, 'low_mem': False, 'projected': True,
                         'max_num_images': 50, 'standardize': False}

        self.__dict__.update(prop_defaults)
        if 'sim_calc_args' in kwargs:
            kwargs = kwargs['sim_calc_args']
        self.__dict__.update(kwargs)

        self.doc = doc
        if kwargs.get('standardize'):
            self.doc = standardize_doc_cell(self.doc)

        self.lattice = np.asarray(doc['lattice_cart'])
        self.poscart = np.asarray(frac2cart(doc['lattice_cart'], doc['positions_frac']))
        self.types = doc['atom_types']
        if 'text_id' in doc:
            self.label = ' '.join(doc['text_id'])
        else:
            self.label = 'null'
        self.num_atoms = len(self.poscart)
        if 'cell_volume' not in doc:
            self.volume = cart2volume(doc['lattice_cart'])
        else:
            self.volume = doc['cell_volume']
        self.number_density = self.num_atoms / self.volume
        if not kwargs.get('lazy'):
            if self.timing:
                start = time.time()
            self.calc_pdf()
            if self.timing:
                end = time.time()
                print('PDF calculated in {:.3f} s'.format(end - start))
            if kwargs.get('projected'):
                self._calc_projected_pdf()

    def _calc_distances(self, poscart, poscart_B=None, debug=False):
        """ Calculate PBC distances with cdist.

        Input:

            | poscart    : np.ndarray, of Cartesian atomic coordinates.

        Args:

            | poscart_B  : np.ndarray, positions of a second type of atoms,
                           where only A-B distances will be calculated.

        Returns:

            | Sets and returns self.distances to d_ij matrix,
              with values > rmax < 1e-12 removed.


        """
        if self.debug:
            start = time.time()
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

        if self.debug:
            end = time.time()
            print('Calculated distances in {} s'.format(end-start))

        return self.distances

    def calc_pdf(self, debug=False):
        """ Wrapper function to calculate distances and output
        a broadened and normalised PDF.

        Returns:

            | Sets self.Gr and self.r_space to G(r) and r respectively.

        """
        self._set_image_trans_vectors()
        distances = self._calc_distances(self.poscart, debug=self.debug)
        self.r_space = np.arange(0, self.rmax+self.dr, self.dr)
        self.Gr = self._set_broadened_normalised_pdf(distances,
                                                     style=self.style,
                                                     gaussian_width=self.gaussian_width)

    def _set_broadened_normalised_pdf(self, distances, style='smear', gaussian_width=0.01):
        """ Broaden the values provided as distances and return
        G(r) and r_space of the normalised PDF.

        Input:

            | distances      : np.ndarray, used to calculate PDF
            | style          : str, either 'smear' or 'histogram'
            | gaussian_width : float, smearing width in Angstrom^1/2

        Requires:

            self.r_space to be set.

        Returns and sets:

            Gr             : np.ndarray, G(r), the PDF of supplied distances

        """

        if self.debug:
            start = time.time()

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

        if self.debug:
            end = time.time()
            print('Calculated broadening and normalised in {} s'.format(end-start))

        return Gr

    def _calc_projected_pdf(self):
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

    def _set_image_trans_vectors(self):
        """ Sets self.image_vec to a list/generator of image translation vectors,
        based on self.num_images.

        If self.num_images is an integer, create all 3-member integer combinations
        up to the value.

        If self.num_images is 'auto', create all translation vectors up to length self.rmax.

        e.g. self.image_vec = [[1, 0, 1], [0, 1, 1], [1, 1, 1]].

        """
        if self.debug:
            start = time.time()

        if self.num_images == 'auto':
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
            first_attempt = 3
            test_num_images = deepcopy(first_attempt)
            while any_in_sphere:
                any_in_sphere = False
                for prod in product(range(-test_num_images, test_num_images+1), repeat=3):
                    if prod in self.image_vec:
                        continue
                    trans = np.zeros((3))
                    for ind, multi in enumerate(prod):
                        trans += self.lattice[ind] * multi
                    if np.sqrt(np.sum(trans**2)) <= self.rmax+self.dr+max_trans:
                        self.image_vec.add(prod)
                        any_in_sphere = True
                test_num_images += 1
                if test_num_images > self.max_num_images:
                    print('Something has probably gone wrong; required images reached {}.'.format(self.max_num_images))
                    print('text_id:')
                    print(self.doc['text_id'])
                    print('lattice_abc:')
                    if 'lattice_abc' in self.doc:
                        print(self.doc['lattice_abc'])
                    else:
                        print(cart2abc(self.doc['lattice_cart']))
                    break
        else:
            self.image_vec = product(range(-self.num_images, self.num_images+1), repeat=3)
        if self.debug:
            end = time.time()
            print('Set image trans vectors in {} s'.format(end-start))

    def get_sim_distance(self, pdf_B, projected=False):
        """ Return the similarity between two PDFs. """
        return PDFOverlap(self, pdf_B, projected=projected).similarity_distance

    def pdf(self):
        """ Return G(r) and the r_space for easy plotting. """
        try:
            return (self.r_space, self.Gr)
        except:
            return (None, None)

    def plot_projected_pdf(self, keys=None, other_pdfs=None):
        """ Plot projected PDFs. """
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8, 5))
        ax1 = fig.add_subplot(111)
        ax1.plot(self.r_space, self.Gr, lw=1, ls='--', label='total')
        if keys is None:
            keys = [key for key in self.elem_Gr]
        for key in keys:
            ax1.plot(self.r_space, self.elem_Gr[key], label='-'.join(key))
        ax1.legend(loc=1)
        ax1.set_ylabel('$g(r)$')
        ax1.set_xlabel('$r$ (Angstrom)')
        plt.show()
        return

    def plot_pdf(self, other_pdfs=None):
        """ Plot projected PDFs, with optional list of
        tuples [(r_space, Gr), ...] of other PDFs.
        """
        import matplotlib.pyplot as plt
        try:
            import seaborn as sns
        except:
            pass
        fig = plt.figure(figsize=(8, 5))
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
        ax1.set_ylabel('Pair distribution function, $g(r)$')
        ax1.set_xlim(0, self.rmax)
        ax1.set_xlabel('$r$ (Angstrom)')
        try:
            sns.despine()
        except:
            pass
        plt.legend()
        plt.show()
        return


class PDFFactory:
    """ This class computes PDF objects from a list of structures,
    as concurrently as possible. The PDFs are stored under the `pdf`
    key inside each structure dict.
    """
    def __init__(self, cursor, debug=False, concurrency='pool', **pdf_args):
        """ Compute PDFs over n processes, where n is set by either
        SLURM_NTASKS, OMP_NUM_THREADS or physical core count.

        Input:

            | cursor: list(dict), list of matador structures

        Args:

            | concurrency: str, either 'pool' or 'queue'
            | pdf_args: dict, arguments to pass to the PDF calculator

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

        import time
        from copy import deepcopy
        start = time.time()
        if concurrency is 'queue':
            queue = mp.Queue()
            # split cursor into subcursors
            pdfs_per_proc = int(round((len(cursor)/self.nprocs)))
            subcursors = []
            for i in range(self.nprocs):
                if i == self.nprocs - 1:
                    subcursors.append(deepcopy(cursor[i*pdfs_per_proc:]))
                else:
                    subcursors.append(deepcopy(cursor[i*pdfs_per_proc:(i+1)*pdfs_per_proc]))

            processes = [mp.Process(target=calc_pdf_queue_wrapper,
                                    args=(subcursors[i], i, queue))
                         for i in range(self.nprocs)]
            for proc in processes:
                proc.start()

            results_cursor = dict()
            [results_cursor.update(queue.get()) for i in range(self.nprocs)]
            assert len(results_cursor) == self.nprocs
            pdf_cursor = []
            for i in range(self.nprocs):
                pdf_cursor.append(results_cursor[i])
            pdf_cursor = [doc for subcursor in pdf_cursor for doc in subcursor]

        elif concurrency is 'pool':
            pool = mp.Pool(processes=self.nprocs)
            pdf_cursor = []
            pool.map_async(calc_pdf_pool_wrapper, cursor, callback=pdf_cursor.extend)
            pool.close()
            pool.join()

        assert len(pdf_cursor) == len(cursor)

        for ind, doc in enumerate(cursor):
            doc['pdf'] = pdf_cursor[ind]['pdf']

        cursor = pdf_cursor
        elapsed = time.time() - start

        if debug:
            print('Compute time: {:.4f} s'.format(elapsed))
            print('Work complete!')


def calc_pdf_queue_wrapper(cursor, i, queue):
    """ Evaluate PDFs of a cursor where a lazy init of each doc's PDF object has
    already been made. The result is parcelled into a dictionary with key i
    and pushed to the queue.

    Input:

        | cursor: list(dict), list of matador structures with empty PDF
                  objects stored under `pdf`.
        | i     : int, position of cursor in overall subcursor array
        | queue : mp.Queue, processing queue.

    """
    for ind, doc in enumerate(cursor):
        cursor[ind]['pdf'].calc_pdf()
    queue.put({i: cursor})


def calc_pdf_pool_wrapper(doc):
    """ Evaluate PDF of a structure where a lazy init of the doc's PDF object has
    already been made.

    Input:

        | doc: dict, matador structures with empty PDF

    """
    doc['pdf'].calc_pdf()
    return {'pdf': doc['pdf']}


class PDFOverlap:
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
        try:
            import seaborn
        except:
            pass
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
        try:
            import seaborn
        except:
            pass
        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(211)
        ax1.plot(np.arange(len(self.convolution), 0, step=-1) * self.fine_dr / 2.0,
                 self.convolution)
        ax1.set_ylabel('$g_A(r) \\ast g_B(r)$')
        ax1.set_xlabel('$\\Delta$ (Angstrom)')
        plt.show()
        return
