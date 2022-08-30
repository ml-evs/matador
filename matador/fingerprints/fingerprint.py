# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the base class for all "fingerprints", which
here refers to any object derived from purely structural features of a
crystal, e.g. pair distribution functions (PDF) or simulated powder
X-ray diffraction (PXRD) spectra.

"""

# TODO: wrap these broadening methods with heuristics to decide which to use

import abc
import multiprocessing as mp
import os
import time
import sys
import functools

import psutil
import numba
import numpy as np

from matador.utils.print_utils import print_notify
from matador.crystal import Crystal


class Fingerprint(abc.ABC):

    fingerprint = None
    default_key = None

    @abc.abstractmethod
    def __init__(self, doc, lazy=True, *args, **kwargs):
        pass

    @abc.abstractmethod
    def calculate(self):
        pass

    @staticmethod
    # @numba.njit
    def _broadening_space_dominated(
        distances, r_space, width, broadening_type="gaussian"
    ):
        """Add broadening to the PDF by convolving distances with
        the radial space and summing. More memory-efficient if len(r_space)
        is less than len(distances).

        Parameters:
            distances (numpy.ndarray): array of pair-wise distances.
            r_space (numpy.ndarray): radial grid
            width (float): amount of broadening.
            broadening_type (str): 'gaussian' or 'lorentzian'.

        Returns:
            gr (numpy.ndarray): the unnormalised PDF.

        """
        new_space = (
            np.reshape(r_space, (1, len(r_space)))
            - np.reshape(distances, (1, len(distances))).T
        )
        if broadening_type == "lorentzian":
            width /= 2
            return np.sum(1 / (1 + (new_space / width) ** 2), axis=0)

        return np.sum(np.exp(-((new_space / width) ** 2)), axis=0)

    @staticmethod
    @numba.njit
    def _broadening_distance_dominated(
        hist, r_space, width, broadening_type="gaussian"
    ):
        """Add broadening to the PDF by convolving the distance histogram with
        the radial space and summing. Potentially more memory-efficient than the alternative
        implementation if len(distances) > len(r_space).

        Parameters:
            hist (numpy.ndarray): histogram of pairwise frequencies.
            r_space (numpy.ndarray): radial grid
            width (float): amount of gaussian broadening.
            broadening_type (str): 'gaussian' or 'lorentzian'.

        Returns:
            gr (numpy.ndarray): the unnormalised PDF.

        """
        new_space = (
            np.reshape(r_space, (1, len(r_space)))
            - np.reshape(r_space, (1, len(r_space))).T
        )

        if broadening_type == "lorentzian":
            width /= 2
            return np.sum(hist / (1 + (new_space / width) ** 2), axis=1)

        return np.sum(hist * np.exp(-((new_space / width) ** 2)), axis=1)

    @staticmethod
    @numba.njit
    def _broadening_unrolled(hist, r_space, width, broadening_type="gaussian"):
        """Add broadening to the PDF by convolving the distance histogram with
        the radial space and summing. Unrolled loop to save memory.


        Parameters:
            hist (numpy.ndarray): histogram of pairwise frequencies.
            r_space (numpy.ndarray): radial grid
            width (float): amount of gaussian broadening.
            broadening_type (str): 'gaussian' or 'lorentzian'.

        Returns:
            gr (numpy.ndarray): the unnormalised PDF.

        """
        gr = np.zeros_like(r_space)

        if broadening_type == "lorentzian":
            width /= 2
            for ind, _ in enumerate(hist):
                if hist[ind] != 0:
                    gr += hist[ind] / (1 + ((r_space - r_space[ind]) / width) ** 2)

        else:
            for ind, _ in enumerate(hist):
                if hist[ind] != 0:
                    gr += hist[ind] * np.exp(-(((r_space - r_space[ind]) / width) ** 2))

        return gr


class FingerprintFactory(abc.ABC):
    """This class computes Fingerprint objects from a list of structures,
    using multiprocessing to perform calculations concurrently. The computed
    fingerprints are stored in each structure's dictionary under the
    default key defined by the Fingerprint objects.

    Note:
        The number of processes used to concurrency is set by the following
        hierarchy:
        ``$SLURM_NTASKS -> $OMP_NUM_THREADS -> psutil.cpu_count(logical=False)``.

    Attributes:
        nprocs (int): number of concurrent processes to be used.

    """

    fingerprint = None
    default_key = None

    def __init__(self, cursor, required_inds=None, debug=False, **fprint_args):
        """Compute PDFs over n processes, where n is set by either
        ``$SLURM_NTASKS``, ``$OMP_NUM_THREADS`` or physical core count.

        Parameters:
            cursor (list of dict): list of matador structures
            fingerprint (Fingerprint): class to compute for each structure

        Keyword arguments:
            pdf_args (dict): arguments to pass to the fingerprint calculator
            required_inds (list(int)): indices in cursor to skip.

        """
        if required_inds is None:
            required_inds = list(range(len(cursor)))
        elif len(required_inds) == 0:
            return
        else:
            print(
                "Skipping {} structures out of {} as no comparisons are required".format(
                    len(cursor) - len(required_inds), len(cursor)
                )
            )

        if self.fingerprint is None or self.default_key is None:
            raise NotImplementedError(
                "Do not create FingerprintFactory directly, "
                "use the appropriate sub-class!"
            )

        # create list of empty (lazy) PDF objects
        if "lazy" in fprint_args:
            del fprint_args["lazy"]

        for ind, doc in enumerate(cursor):
            if isinstance(doc, Crystal):
                doc._data.pop(self.default_key, None)
            if ind in required_inds:
                doc[self.default_key] = self.fingerprint(doc, lazy=True, **fprint_args)
            else:
                doc[self.default_key] = None

        compute_list = [doc for ind, doc in enumerate(cursor) if ind in required_inds]

        # how many processes to use? either SLURM_NTASKS, OMP_NUM_THREADS or total num CPUs
        if os.environ.get("SLURM_NTASKS") is not None:
            self.nprocs = int(os.environ.get("SLURM_NTASKS"))
            env = "$SLURM_NTASKS"
        elif os.environ.get("OMP_NUM_THREADS") is not None:
            self.nprocs = int(os.environ.get("OMP_NUM_THREADS"))
            env = "$OMP_NUM_THREADS"
        else:
            self.nprocs = psutil.cpu_count(logical=False)
            env = "core count"
        print_notify(
            "Running {} jobs on at most {} processes, set by {}.".format(
                len(required_inds), self.nprocs, env
            )
        )
        self.nprocs = min(len(compute_list), self.nprocs)

        start = time.time()
        if self.nprocs == 1:
            import tqdm

            for ind, doc in tqdm.tqdm(enumerate(cursor)):
                if cursor[ind][self.default_key] is not None:
                    cursor[ind][self.default_key].calculate()
        else:
            pool = mp.Pool(processes=self.nprocs)
            fprint_cursor = []
            # for large cursors, set chunk to at most 16
            # for smaller cursors, tend to use chunksize 1 for improved load balancing
            chunksize = min(max(1, int(0.25 * len(compute_list) / self.nprocs)), 16)
            results = pool.map_async(
                functools.partial(_calc_fprint_pool_wrapper, key=self.default_key),
                compute_list,
                callback=fprint_cursor.extend,
                error_callback=print,
                chunksize=chunksize,
            )
            pool.close()
            width = len(str(len(required_inds)))
            total = len(required_inds)
            while not results.ready():
                sys.stdout.write(
                    "{done:{width}d} / {total:{width}d}  {percentage:3d}%\r".format(
                        width=width,
                        done=total - results._number_left,
                        total=total,
                        percentage=int(100 * (total - results._number_left) / total),
                    )
                )
                sys.stdout.flush()
                time.sleep(1)

            if len(fprint_cursor) != len(required_inds):
                raise RuntimeError(
                    "There was an error calculating the desired Fingerprint"
                )

            fprint_ind = 0
            for ind, doc in enumerate(cursor):
                if ind in required_inds:
                    if isinstance(cursor[ind], Crystal):
                        cursor[ind]._data.pop(self.default_key, None)
                    cursor[ind][self.default_key] = fprint_cursor[fprint_ind][
                        self.default_key
                    ]
                    fprint_ind += 1

        elapsed = time.time() - start
        if debug:
            pool.close()
            print("Compute time: {:.4f} s".format(elapsed))
            print("Work complete!")


def _calc_fprint_pool_wrapper(doc, key=None):
    """Evaluate Fingerprint of a structure where a lazy init of the
    doc's Fingerprint object has already been made.

    Parameters:
        doc (dict): matador structures with empty PDF

    """
    doc[key].calculate()
    return {key: doc[key]}
