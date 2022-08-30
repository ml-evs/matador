# coding: utf-8
# Distributed under the terms of the MIT License.

""" The workflows module contains ways of creating custom workflows for
chaining together calculations. Each custom workflow inherits from
:obj:`workflows.Workflow` and consists of several
:obj:`workflows.WorkflowStep` objects.
If the creation of the :obj:`Workflow` is wrapped in a function with
the correct signature, a whole :obj:`Workflow` can itself be used as a
:obj:`workflows.WorkflowStep`.


"""


__all__ = ["Workflow", "WorkflowStep"]

__author__ = "Matthew Evans"
__maintainer__ = "Matthew Evans"


from .workflows import Workflow, WorkflowStep
