# coding: utf-8
# Distributed under the terms of the MIT License.

""" The compute module contains three submodules, compute, batch and
slurm.

The compute submodule contains the :class:`ComputeTask` class for performing
continually restarted geometry optimisation and SCF calculations in
CASTEP, as well as the execution of arbitrary programs with mpirun.

The batch submodule contains the :class:`BatchRun` class for running several
independent :class:`ComputeTask` instances on a folder of structures, without
clashes.

The slurm submodule provides a wrapper to useful slurm commands, and to
writing slurm job submission files.

"""


__all__ = ['ComputeTask', 'BatchRun', 'reset_job_folder']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.compute.compute import ComputeTask
from matador.compute.batch import BatchRun, reset_job_folder
