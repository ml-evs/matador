# coding: utf-8
# Distributed under the terms of the MIT License.

""" The compute module contains two classes:

* the FullRelaxer class for performing continuously restarted
geometry optimisation and SCF calculations in CASTEP, as well
as the execution of arbitrary programs with mpirun.

* the BatchRun class for running several independent FullRelaxer instances
on a folder of structures, without clashes.

"""


__all__ = ['FullRelaxer', 'BatchRun']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.compute.compute import FullRelaxer, BatchRun
