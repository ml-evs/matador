# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements a simple interface to basic PBS
functionality, mostly for monitoring walltime of jobs submitted
via PBS.

"""

from matador.compute.queueing import QueueManager


class PBSQueueManager(QueueManager):
    """Wrapper for the PBS queueing system."""

    token = "pbs"

    def get_ntasks(self):
        return int(self.env["PBS_TASKNUM"])

    def get_max_memory(self):
        return None

    def get_array_id(self):
        if self.env.get("PBS_ARRAYID") is not None:
            return int(self.env["PBS_ARRAYID"])
        return None

    def get_walltime(self):
        """Query available walltime with qstat on the current job.

        Parameters:
            pbs_dict (dict): pbs env parameters to query.

        Raises:
            RuntimeError: if PBS_JOBID not present in slurm env.
            subprocess.CalledProcessError: if unable to use qstat.

        Returns:
            int: maximum allowed walltime time in seconds.

        """
        import subprocess as sp

        pbs_dict = self.env
        job_id = pbs_dict.get("PBS_JOBID")
        if job_id is not None:
            output = (
                sp.check_output("qstat -f {}".format(job_id), shell=True)
                .decode("utf-8")
                .split("\n")
            )

        output_dict = {
            line.strip()
            .split("=")[0]
            .strip(): " ".join(line.strip().split("=")[1:])
            .strip()
            for line in output
            if "=" in line
        }

        walltime = output_dict.get("Resource_List.walltime")

        hrs = 0
        if "-" in walltime:
            days = int(walltime.split("-")[0])
            walltime = walltime.split("-")[1]
            hrs += days * 24

        hrs += int(walltime.split(":")[0])
        mins = int(walltime.split(":")[1])
        secs = int(walltime.split(":")[2])
        walltime_in_seconds = (60 * hrs + mins) * 60 + secs

        return walltime_in_seconds
