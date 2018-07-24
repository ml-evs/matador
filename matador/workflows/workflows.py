# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements various workflows, ways
of chaining up different calculations at high-throughput.

"""

import logging
import abc


class Workflow:
    """ Workflow objects are bundles of calculations defined as
    WorkflowStep objects. Each WorkflowStep takes three arguments:
    the FullRelaxer object used to run the calculations, the calculation
    parameters (which can be modified by each step), the seed name.
    Any subclass of Workflow must implement `preprocess` and `postprocess`
    methods (even if they just return True).

    Attributes:
        relaxer (FullRelaxer): the object that calls CASTEP.
        calc_doc (dict): the interim dictionary of structural and
            calculation parameters.
        seed (str): the root seed for the calculation.
        label (str): the name of the type of the Workflow object.
        success (bool): the status of the Workflow: only set to True after
            final Workflow step completes, BEFORE post-processing.
        steps (:obj:`list` of :obj:`WorkflowStep`): list of steps to be
            completed.

    """
    def __init__(self, relaxer, calc_doc, seed):
        """ Initialise the Workflow object from a FullRelaxer, calculation
        parameters and the seed name.

        Parameters:
            relaxer (FullRelaxer): the object that will be calling CASTEP.
            calc_doc (dict): dictionary of structure and calculation
                parameters.
            seed (str): root seed for the calculation.

        Raises:
            RuntimeError: if any part of the calculation fails.

        """
        self.relaxer = relaxer
        self.calc_doc = calc_doc
        self.seed = seed
        self.label = self.__class__.__name__
        self.success = None
        self.steps = []

        logging.info('Performing Workflow of type {} on {}'.format(self.label, self.seed))

        self.preprocess()
        self.run_steps()
        self.postprocess()

    @abc.abstractmethod
    def preprocess(self):
        """ This function is run at the start of the workflow, and is
        responsible for adding WorkflowStep objects to the Workflow.

        """
        raise NotImplementedError('Please implement a preprocess method.')

    def postprocess(self):
        """ This function is run upon successful completion of all steps
        of the workflow and is responsible for cleaning up files and any
        other post-processing.

        """
        if self.success:
            logging.info('Writing results of Workflow {} run to res file and tidying up.'.format(self.label))
            self.relaxer.mv_to_completed(self.seed, keep=True)
        else:
            logging.info('Writing results of failed Workflow {} run to res file and tidying up.'.format(self.label))
            self.relaxer.mv_to_bad(self.seed)

    def add_step(self, function, name):
        """ Add a step to the workflow.

        Parameters:
            function (Function): the function to run in the step; must
                accept arguments of (self.relaxer, self.calc_doc, self.seed).
            name (str): the desired name for the step (human-readable).

        """
        self.steps.append(WorkflowStep(function, name))

    def run_steps(self):
        """ Loop over steps and run them. """
        try:
            if not self.steps:
                msg = 'No steps added to Workflow!'
                logging.error(msg)
                raise RuntimeError(msg)

            for step in self.steps:
                step.run_step(self.relaxer, self.calc_doc, self.seed)
            self.success = True

        except RuntimeError:
            self.success = False
            msg = '{} workflow exiting...'.format(self.label)
            logging.error(msg)
            raise RuntimeError(msg)


class WorkflowStep:
    """ An individual step in a Workflow, defined by a Python function
    and a name. The function will be called with arguments
    (relaxer, calc_doc, seed) with the run_step method.

    """
    def __init__(self, function, name):
        """ Construct a WorkflowStep from a function. """
        logging.debug('Constructing WorkflowStep: {}'.format(name))
        self.function = function
        self.name = name

    def run_step(self, relaxer, calc_doc, seed):
        """ Run the workflow step.

        Parameters:
            relaxer (FullRelaxer): the object that will be calling CASTEP.
            calc_doc (dict): dictionary of structure and calculation
                parameters.
            seed (str): root seed for the calculation.

        Raises:
            RuntimeError: if any step fails.

        """
        try:
            logging.info('{} starting...'.format(self.name))
            self.function(relaxer, calc_doc, seed)
            logging.info('{} completed successfully.'.format(self.name))

        except RuntimeError as exc:
            msg = '{} step failed with error {}.'.format(self.name, exc)
            logging.error(msg)
            raise RuntimeError(msg)
