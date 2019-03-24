# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the base class for all "fingerprints", which
here refers to any object derived from purely structural features of a
crystal, e.g. pair distribution functions (PDF) or simulated powder
X-ray diffraction (PXRD) spectra.

"""

import abc
import multiprocessing as mp
import os
import time
import sys
import numba
import numpy as np

from matador.utils.print_utils import print_notify


class Fingerprint(abc.ABC):

    fingerprint = None
    default_key = None
    @abc.abstractmethod
    def __init__(self, doc, lazy=True, *args, **kwargs):
        pass

    @abc.abstractmethod
    def calculate(self):
        pass

    # TODO: wrap these broadening methods with heuristics to decide which to use
    # TODO: Lorentzian broadening as an option

    @staticmethod
    @numba.njit
    def _broadening_space_dominated(distances, r_space, gaussian_width):
        """ Add Gaussian broadening to the PDF by convolving distances with
        the radial space and summing. More memory-efficient if len(r_space)
        is less than len(distances).

        Parameters:
            distances (numpy.ndarray): array of pair-wise distances.
            r_space (numpy.ndarray): radial grid
            gaussian_width (float): amount of gaussian broadening.

        Returns:
            gr (numpy.ndarray): the unnormalised PDF.

        """
        new_space = (np.reshape(r_space, (1, len(r_space))) -
                     np.reshape(distances, (1, len(distances))).T)
        gr = np.sum(np.exp(-(new_space / gaussian_width)**2), axis=0)
        return gr

    @staticmethod
    @numba.njit
    def _broadening_distance_dominated(hist, r_space, gaussian_width):
        """ Add Gaussian broadening to the PDF by convolving the distance histogram with
        the radial space and summing. Potentially more memory-efficient than the alternative
        implementation if len(distances) > len(r_space).

        Parameters:
            hist (numpy.ndarray): histogram of pairwise frequencies.
            r_space (numpy.ndarray): radial grid
            gaussian_width (float): amount of gaussian broadening.

        Returns:
            gr (numpy.ndarray): the unnormalised PDF.

        """
        new_space = (np.reshape(r_space, (1, len(r_space))) -
                     np.reshape(r_space, (1, len(r_space))).T)
        gr = np.sum(hist * np.exp(-(new_space / gaussian_width)**2), axis=1)
        return gr

    @staticmethod
    @numba.njit
    def _broadening_unrolled(hist, r_space, gaussian_width):
        """ Add Gaussian broadening to the PDF by convolving the distance histogram with
        the radial space and summing. Unrolled loop to save memory.


        Parameters:
            hist (numpy.ndarray): histogram of pairwise frequencies.
            r_space (numpy.ndarray): radial grid
            gaussian_width (float): amount of gaussian broadening.

        Returns:
            gr (numpy.ndarray): the unnormalised PDF.

        """
        gr = np.zeros_like(r_space)
        for ind, _ in enumerate(hist):
            if hist[ind] != 0:
                gr += hist[ind] * np.exp(-((r_space-r_space[ind]) / gaussian_width)**2)
        return gr


class FingerprintFactory:
    """ This class computes Fingerprint objects from a list of structures,
    using multiprocessing to perform calculations concurrently. The computed
    fingerprints are stored in each structure's dictionary under the
    default key defined by the Fingerprint objects.

    Note:
        The number of processes used to concurrency is set by the following
        hierarchy:
        SLURM_NTASKS -> OMP_NUM_THREADS -> multiprocessing.cpu_count().

    Attributes:
        nprocs (int): number of concurrent processes to be used.

    """
    # TODO: how do you prevent this init from being called when fingerprint is None?

    def __init__(self, cursor, debug=False, **fprint_args):
        """ Compute PDFs over n processes, where n is set by either
        SLURM_NTASKS, OMP_NUM_THREADS or physical core count.

        Parameters:
            cursor (list of dict): list of matador structures

        Keyword arguments:
            pdf_args (dict): arguments to pass to the PDF calculator

        """
        # create list of empty (lazy) PDF objects
        if 'lazy' in fprint_args:
            del fprint_args['lazy']

        for doc in cursor:
            doc[self.default_key] = self.fingerprint(doc, lazy=True, **fprint_args)

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
        print_notify('Running {} jobs on {} processes, set by {}.'.format(len(cursor),
                                                                          self.nprocs,
                                                                          env))

        start = time.time()
        if self.nprocs == 1:
            import tqdm
            for ind, doc in tqdm.tqdm(enumerate(cursor)):
                cursor[ind][self.default_key].calculate()
        else:
            pool = mp.Pool(processes=self.nprocs)
            import functools
            fprint_cursor = []
            results = pool.map_async(functools.partial(calc_fprint_pool_wrapper,
                                                       self.default_key),
                                     cursor,
                                     callback=fprint_cursor.extend,
                                     error_callback=print,
                                     # set chunksize to 1 as Fingerprints will often be O(N^2)
                                     # in the number of atoms: alternatively sort by N at
                                     # the start
                                     chunksize=1)
            pool.close()
            width = len(str(len(cursor)))
            total = len(cursor)
            while not results.ready():
                sys.stdout.write('{done:{width}d} / {total:{width}d}  {percentage:3d}%\r'
                                 .format(width=width, done=total-results._number_left,
                                         total=total, percentage=int(100*(total-results._number_left)/total)))
                sys.stdout.flush()
                time.sleep(1)

            if len(fprint_cursor) != len(cursor):
                raise RuntimeError('There was an error calculating the desired Fingerprint')

            for ind, doc in enumerate(cursor):
                cursor[ind][self.default_key] = fprint_cursor[ind][self.default_key]

        elapsed = time.time() - start
        if debug:
            print('Compute time: {:.4f} s'.format(elapsed))
            print('Work complete!')


def calc_fprint_pool_wrapper(doc, key):
    """ Evaluate Fingerprint of a structure where a lazy init of the
    doc's PDF object has already been made.

    Parameters:
        doc (dict): matador structures with empty PDF

    """
    doc[key].calculate()
    return {key: doc[key]}
