# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements a simple interface to basic SLURM
functionality, including creating and submitting slurm scripts and
cancelling jobs.

"""


def get_slurm_env(fail_loudly=True):
    """ Scrape SLURM environment variables from current env.
    This function can be used when called inside an active slurm job.

    Keyword arguments:
        fail_loudly (bool): option to raise SystemExit if SLURM not detected.

    Returns:
        dict: dictionary containing all the currently set SLURM environment variables.

    """
    from os import environ
    slurm_dict = {key: environ[key] for key in environ if 'slurm' in key.lower()}
    if not slurm_dict and fail_loudly:
        exit('Requested SLURM array mode, yet no SLURM settings were found. Was this process submitted as a job?')
    return slurm_dict


def get_slurm_walltime(slurm_dict):
    """ Query available walltime with scontrol on the current job.

    Parameters:
        slurm_dict (dict): slurm env parameters to query.

    Raises:
        RuntimeError: if SLURM_JOB_ID not present in slurm env.
        subprocess.CalledProcessError: if unable to use scontrol.

    Returns:
        int: maximum allowed walltime time in seconds.

    """
    import subprocess as sp
    job_id = slurm_dict.get('SLURM_JOB_ID')
    if job_id is not None:
        output = sp.check_output('scontrol show job={}'.format(job_id), shell=True).decode('utf-8')

    output_dict = {line.split('=')[0].lower(): line.split('=')[-1] for line in output.split()}

    walltime = output_dict.get('timelimit')
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


def scancel_all_matching_jobs(name=None):
    """ Cancel all of the user's jobs.

    Keyword arguments:
        name (str): optional name to pass to scancel

    Returns:
        str: output from scancel.

    """
    from os import getlogin
    import subprocess as sp
    user = getlogin()
    if name is None:
        return sp.check_output('scancel -u {}'.format(user), shell=True).decode('utf-8')

    return sp.check_output('scancel -u {} -n {}'.format(user, name), shell=True).decode('utf-8')


def submit_slurm_script(slurm_fname, depend_on_job=None, num_array_tasks=None):
    """ Submit a SLURM job.

    Parameters:
        slurm_fname (str): SLURM job file to submit.

    Keyword arguments:
        depend_on_job (int): job ID to make current job depend on.
        num_array_tasks (int): number of array tasks to submit.

    Raises:
        subprocess.CalledProcessError: if jobfile doesn't exist or has failed.

    Return:
        int: submitted SLURM job ID.

    """
    import subprocess as sp
    command = 'sbatch '
    if depend_on_job is not None:
        command += '--dependency=afterany:{} '.format(depend_on_job)
    if num_array_tasks is not None:
        assert num_array_tasks > 0
        if num_array_tasks != 1:
            command += '--array=0-{} '.format(num_array_tasks-1)
    command += '{}'.format(slurm_fname)
    slurm_output = sp.check_output(command, shell=True).decode('utf-8')
    slurm_job_id = int(slurm_output.strip().split()[-1])
    return slurm_job_id


def get_slurm_header(slurm_dict, walltime_hrs, num_nodes=None):
    """ Write a SLURM script header from a set of slurm parameters.

    Parameters:
        slurm_dict (dict): dictionary of SLURM environment variables.
        walltime_hrs (int): allowed walltime in hours

    Keyword arguments:
        num_nodes (int): overrides $SLURM_JOB_NUM_NODES with a custom value.

    Returns:
        header (str): the SLURM file header.

    """

    header = "#!/bin/bash\n"
    header += "#! SLURM file written by matador (Matthew Evans 2016-2018).\n\n"
    header += "#! Name of job:\n"
    header += "#SBATCH --job-name {}\n".format(slurm_dict['SLURM_JOB_NAME'])
    header += "#! Name of project:\n"
    header += "#SBATCH --account {}\n".format(slurm_dict['SLURM_JOB_ACCOUNT'])
    if num_nodes is None:
        num_nodes = slurm_dict['SLURM_JOB_NUM_NODES']
    header += "#! Number of nodes to allocate:\n"
    header += "#SBATCH --nodes {}\n".format(num_nodes)
    header += "#! Number of tasks to allocate:\n"
    header += "#SBATCH --ntasks {}\n".format(slurm_dict['SLURM_NTASKS'])
    header += "#! Partition:\n"
    header += "#SBATCH --partition {}\n".format(slurm_dict['SLURM_JOB_PARTITION'])
    header += "#! Walltime to allocate:\n"
    header += "#SBATCH --time {}:00:00\n".format(walltime_hrs)

    return header


def write_slurm_submission_script(slurm_fname, slurm_dict, compute_string, walltime_hrs, template=None):
    """ Write a full slurm submission script based on the
    input settings.

    Parameters:
        slurm_fname (str): the desired filename for the submission script
        slurm_dict (dict): dictionary of SLURM environment variables
        compute_string (str): the compute commands to run
        walltime_hrs (int): maximum walltime in hours

    Keyword arguments:
        template (str): filename containing job preamble, e.g. module loads

    """
    header = get_slurm_header(slurm_dict, walltime_hrs)
    if template is not None:
        with open(template, 'r') as f:
            preamble = f.readlines()
    else:
        preamble = []

    with open(slurm_fname, 'w') as f:
        f.write(header)
        f.write('\n\n')
        for line in preamble:
            f.write(line)
        f.write('\n\n')
        f.write(compute_string)
