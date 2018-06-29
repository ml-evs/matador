# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements a simple interface to basic PBS
functionality, mostly for monitoring walltime of jobs submitted
via PBS.

"""

from matador.compute.queue import get_queue_env


def get_pbs_env():
    """ Scrape PBS environment variables from current env.
    This function can be used when called inside an active PBS job.

    Returns:
        dict: dictionary containing all the currently set SLURM environment variables.

    """
    return get_queue_env(token='pbs')


def get_pbs_walltime(pbs_dict):
    """ Query available walltime with qstat on the current job.

    Parameters:
        pbs_dict (dict): pbs env parameters to query.

    Raises:
        RuntimeError: if PBS_JOBID not present in slurm env.
        subprocess.CalledProcessError: if unable to use qstat.

    Returns:
        int: maximum allowed walltime time in seconds.

    """
    import subprocess as sp
    job_id = pbs_dict.get('PBS_JOBID')
    if job_id is not None:
        output = sp.check_output('qstat -f {}'.format(job_id), shell=True).decode('utf-8').split('\n')

    output_dict = {
        line.strip().split('=')[0].strip(): ' '.join(line.strip().split('=')[1:]).strip()
        for line in output if '=' in line
    }

    walltime = output_dict.get('Resource_List.walltime')

    hrs = 0
    if '-' in walltime:
        days = int(walltime.split('-')[0])
        walltime = walltime.split('-')[1]
        hrs += days * 24

    hrs += int(walltime.split(':')[0])
    mins = int(walltime.split(':')[1])
    secs = int(walltime.split(':')[2])
    walltime_in_seconds = (60 * hrs + mins) * 60 + secs

    return walltime_in_seconds
