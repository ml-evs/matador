# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements a simple agnostic interface to various
queueing systems.

"""

import os
import abc


class QueueManager(abc.ABC):
    """Abstract base class for queue managers."""

    token = None

    def __repr__(self):
        return "{} with env: {}".format(type(self), self.env)

    def __init__(self):
        self.env = get_queue_env(self.token)

    @property
    def walltime(self):
        """Return the allotted walltime in seconds,
        returning None if not available.

        """
        try:
            return self.get_walltime()
        except Exception:
            None

    @property
    def ntasks(self):
        return self.get_ntasks()

    @property
    def max_memory(self):
        """Return the allotted memory in MB,
        returning None if not available.

        """
        try:
            return self.get_max_memory()
        except Exception:
            return None

    @property
    def array_id(self):
        return self.get_array_id()

    @abc.abstractmethod
    def get_walltime(self):
        pass

    @abc.abstractmethod
    def get_max_memory(self):
        pass

    @abc.abstractmethod
    def get_array_id(self):
        pass


def get_queue_env(token):
    """Read os.environment variables for either PBS or SLURM
    prefixes, and return a dictionary of those vars only.

    Parameter:
        token (str): choose one of either SLURM or PBS explicitly.

    Returns:
        (dict, str): dictionary of keys from the detected/specified queue, and
            a string containing either "slurm" or "pbs".

    """
    return {key: os.environ[key] for key in os.environ if key.lower().startswith(token)}


def get_queue_manager():
    """Detects whether PBS, SLURM or neither is being used
    by probing the environment variables SLURM_NTASKS and
    PBS_TASKNUM.

    Returns:
        str or None: either "slurm", "pbs" or None.

    Raises:
        SystemExit: if both SLURM and PBS were found.

    """
    queue_mgr = []
    if os.environ.get("SLURM_NTASKS") is not None:
        queue_mgr.append("slurm")
    if os.environ.get("PBS_TASKNUM") is not None:
        queue_mgr.append("pbs")

    if len(queue_mgr) > 1:
        raise RuntimeError("Both SLURM and PBS were found... aah!")
    if not queue_mgr:
        return None

    if queue_mgr[0] == "slurm":
        from matador.compute.slurm import SlurmQueueManager

        return SlurmQueueManager()

    if queue_mgr[0] == "pbs":
        from matador.compute.pbs import PBSQueueManager

        return PBSQueueManager()
