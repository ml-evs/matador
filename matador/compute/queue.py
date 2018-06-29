# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements a simple agnostic interface to various
queueing systems.

"""

import os


def get_queue_env(token=None):
    """ Read os.environment variables for either PBS or SLURM
    prefixes, and return a dictionary of those vars only.

    Keyword arguments:
        token (str): choose one of either SLURM or PBS explicitly.
Returns:
        (dict, str): dictionary of keys from the detected/specified queue, and
            a string containing either "slurm" or "pbs".

    """

    if token is not None:
        queue_mgr = token.lower()
    else:
        queue_mgr = get_queue_manager()

    return {key: os.environ[key] for key in os.environ if key.lower().startswith(queue_mgr)}


def get_queue_walltime(queue_env, queue_mgr):
    """ Query the current job in the detect queuing system in
    order to scrape its walltime.

    Parameters:
        queue_env (dict): dictionary of queue parameters.
        queue_mgr (str): either 'slurm' or 'pbs'.

    Returns:
        int: number of seconds allocated to job.

    """
    if queue_mgr == 'slurm':
        from matador.compute.slurm import get_slurm_walltime
        return get_slurm_walltime(queue_env)
    elif queue_mgr == 'pbs':
        from matador.compute.pbs import get_pbs_walltime
        return get_pbs_walltime(queue_env)
    else:
        raise SystemExit('Unable to detect queue.')


def get_queue_manager():
    """ Detects whether PBS, SLURM or neither is being used
    by probing the environment variables SLURM_NTASKS and
    PBS_TASKNUM.

    Returns:
        str or None: either "slurm", "pbs" or None.

    Raises:
        SystemExit: if both SLURM and PBS were found.

    """
    queue_mgr = []
    if os.environ.get('SLURM_NTASKS') is not None:
        queue_mgr.append('slurm')
    if os.environ.get('PBS_TASKNUM') is not None:
        queue_mgr.append('pbs')

    if len(queue_mgr) > 1:
        raise SystemExit('Both SLURM and PBS were found... aah!')
    elif not queue_mgr:
        return None
    else:
        return queue_mgr[0]
