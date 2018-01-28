""" This file implements a simple interface
to creating and submitting slurm scripts.
"""


def get_slurm_env(fail_loudly=True):
    """ Scrape SLURM environment variables from current env.
    This function can be used when called inside an active slurm job.

    Args:

        | fail_loudly: bool, raise SystemExit if SLURM not detected.

    Returns:

        | slurm_dict: dict, dictionary containing all the currently set SLURM environment variables.

    """
    from os import environ
    slurm_dict = {key: environ[key] for key in environ if 'slurm' in key.lower()}
    if len(slurm_dict) == 0 and fail_loudly:
        exit('Requested SLURM array mode, yet no SLURM settings were found. Was this process submitted as a job?')
    return slurm_dict


def scancel_all_matching_jobs(name=None):
    """ Cancel all of the user's jobs.

    Args:

        | name: str, optional name to pass to scancel

    Returns:

        | slurm_output: str, output from scancel.

    """
    from os import getlogin
    import subprocess as sp
    user = getlogin()
    if name is None:
        return sp.check_output('scancel -u {}'.format(user), shell=True).decode('utf-8')
    else:
        return sp.check_output('scancel -u {} -n {}'.format(user, name), shell=True).decode('utf-8')


def submit_slurm_script(slurm_fname, depend_on_job=None, num_array_tasks=None):
    """ Submit a SLURM job.

    Input:

        | slurm_fname: str, SLURM job file to submit.

    Args:

        | depend_on_job: int, job ID to make current job depend on.
        | num_array_tasks  : int, number of array tasks to submit.

    Raises:

        | subprocess.CalledProcessError: if jobfile doesn't exist or has failed.

    Return:

        | slurm_job_id: int, submitted SLURM job ID.

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

    Input:

        | slurm_dict: dict, dictionary of SLURM environment variables.
        | walltime_hrs: int, allowed walltime in hours

    Args:

        | num_nodes: int, overrides $SLURM_JOB_NUM_NODES with a custom value.

    Returns:

        | header: str, the SLURM file header.

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


def write_slurm_submission_script(slurm_fname, slurm_dict, compute_string, walltime_hrs,
                                  template=None, num_nodes=None):
    """ Write a full slurm submission script based on the
    input settings.

    Input:

        | slurm_fname    : str, the desired filename for the submission script
        | slurm_dict     : dict, dictionary of SLURM environment variables
        | compute_string : str, the compute commands to run
        | walltime_hrs   : int, maximum walltime in hours

    Args:

        | template: str, filename containing job preamble, e.g. module loads
        | num_nodes: int, override SLURM settings for maximum number of nodes per job

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
