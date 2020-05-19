# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements various workflows, ways
of chaining up different calculations at high-throughput.

"""

import abc
import os
import logging
from matador.utils.print_utils import dumps

LOG = logging.getLogger('run3')


class Workflow:
    """ Workflow objects are bundles of calculations defined as
    :obj:`WorkflowStep` objects. Each :obj:`WorkflowStep` takes three arguments:
    the :obj:`matador.compute.ComputeTask` object used to run the calculations, the calculation
    parameters (which can be modified by each step), the seed name.
    Any subclass of Workflow must implement `preprocess` and `postprocess`
    methods (even if they just return True).

    Attributes:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be running the computation.
        calc_doc (dict): the interim dictionary of structural and
            calculation parameters.
        seed (str): the root seed for the calculation.
        label (str): the name of the type of the Workflow object.
        success (bool): the status of the workflow. This is only set to True after
            final step completes, but BEFORE post-processing.
        steps (:obj:`list` of :obj:`WorkflowStep`): list of steps to be
            completed.

    """
    def __init__(self, relaxer, calc_doc, seed, **workflow_kwargs):
        """ Initialise the Workflow object from a :obj:`matador.compute.ComputeTask`, calculation
        parameters and the seed name.

        Parameters:
            relaxer (:obj:`matador.compute.ComputeTask`): the object that will be running the computation.
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
        self.compute_dir = relaxer.compute_dir
        self.success = None
        self.steps = []
        self.clean_after_step = []
        self.workflow_params = workflow_kwargs

        LOG.info('Performing Workflow of type {} on {}'.format(self.label, self.seed))

        self.preprocess()
        try:
            self.run_steps()
        except RuntimeError as exc:
            LOG.critical('Workflow failed: calling postprocess()')
            self._clean_up()
            raise exc

        if self.success:
            self.postprocess()

        self._clean_up()

    @abc.abstractmethod
    def preprocess(self):
        """ This function is run at the start of the workflow, and is
        responsible for adding WorkflowStep objects to the Workflow.

        """
        raise NotImplementedError('Please implement a preprocess method.')

    def postprocess(self):
        """ This OPTIONAL function is run upon successful completion of all steps
        of the workflow and can be overloaded by the subclass to perform
        any postprocessing steps. This occurs *before* cleaning up the
        directory (i.e. moving to completed/bad_castep).

        """

    def _clean_up(self, success=None):
        """ This method moves files to `completed/` or `bad_castep/` depending on the status
        of the workflow. It will use the current seed of the relaxer, so this function can be
        called at intermediate steps if this seed changes.

        """
        cwd = os.getcwd()
        if success is None:
            success = self.success

        if self.compute_dir:
            os.chdir(self.compute_dir)
            if success:
                LOG.info('Writing results from compute dir of Workflow {} run to completed folder and tidying up.'
                         .format(self.label))
                self.relaxer.mv_to_completed(self.relaxer.seed, keep=True, skip_existing=False)
            else:
                LOG.info('Writing results from compute dir of failed Workflow {} run to bad_castep folder and tidying up.'
                         .format(self.label))
                self.relaxer.mv_to_bad(self.relaxer.seed)
            os.chdir(cwd)

        if success:
            LOG.info('Writing results of Workflow {} run to completed folder and tidying up.'.format(self.label))
            self.relaxer.mv_to_completed(self.relaxer.seed, keep=True, skip_existing=True)
        else:
            LOG.info('Writing results of failed Workflow {} run to bad_castep folder and tidying up.'.format(self.label))
            self.relaxer.mv_to_bad(self.relaxer.seed)

    def add_step(self, function, name, input_exts=None, output_exts=None, clean_after=False, **func_kwargs):
        """ Add a step to the workflow.

        Parameters:
            function (Function): the function to run in the step; must
                accept arguments of (self.relaxer, self.calc_doc, self.seed).
            name (str): the desired name for the step (human-readable).

        Keyword arguments:
            clean_after (bool): whether or not to clean up after this step is called
            func_kwargs (dict): any arguments to pass to function when called.

        """
        self.steps.append(WorkflowStep(function, name,
                                       self.compute_dir, input_exts, output_exts,
                                       **func_kwargs))
        self.clean_after_step.append(clean_after)

    def run_steps(self):
        """ Loop over steps and run them. """
        try:
            if not self.steps:
                msg = 'No steps added to Workflow!'
                LOG.error(msg)
                raise RuntimeError(msg)

            for ind, step in enumerate(self.steps):
                LOG.info("Running step {step.name}: {step.function}".format(step=step))
                LOG.debug("Current state: " + dumps(self.calc_doc, indent=None))
                success = step.run_step(self.relaxer, self.calc_doc, self.seed)
                if self.clean_after_step[ind]:
                    self._clean_up(success=success)

            self.success = True

        except RuntimeError:
            self.success = False
            msg = '{} workflow exiting...'.format(self.label)
            LOG.error(msg)
            raise RuntimeError(msg)


class WorkflowStep:
    """ An individual step in a Workflow, defined by a Python function
    and a name. The function will be called with arguments
    (relaxer, calc_doc, seed) with the run_step method.

    Attributes:
        function (function): the function to call.
        name (str): the human-readable name of the step.
        compute_dir (str): the folder that relaxer will perform the calculation in.
        func_kwargs (dict): any extra kwargs to pass to the function.
        input_exts (list): list of input file extensions to cache after running.
        output_exts (list): list of output file extensions to cache after running.

    """

    success = False

    def __init__(self, function, name, compute_dir=None, input_exts=None, output_exts=None, **func_kwargs):
        """ Construct a WorkflowStep from a function. """
        LOG.debug('Constructing WorkflowStep: {}'.format(name))
        self.function = function
        self.name = name
        self.compute_dir = compute_dir
        self.func_kwargs = func_kwargs
        self.input_exts = input_exts
        self.output_exts = output_exts

    def _cache_files(self, seed, exts, mode, directory=None):
        """ Copy any files <seed>.<ext> for ext in exts to
        <seed>.<ext>_<label>.

        Parameters:
            seed (str): seed for the workflow step.
            exts (:obj:`list` of :obj:`str`): list of file extensions, including '.'.
            mode (str): either 'in' (warning printed if file missing) or 'out' (no warning).

        """
        import shutil
        import glob

        for ext in exts:
            if '*' in ext:
                srcs = glob.glob('{}{}'.format(seed, ext))
            else:
                srcs = ['{}{}'.format(seed, ext)]
            for src in srcs:
                dst = src + '_{}'.format(self.name)
                if os.path.isfile(src):
                    shutil.copy2(src, dst, follow_symlinks=True)
                    LOG.info('Backed up {} file {} to {}.'.format(mode, src, dst))
                else:
                    if mode == 'in':
                        error = 'Failed to cache input file {} for step {}.'.format(src, self.name)
                        LOG.warning(error)

    def _cache_inputs(self, seed):
        """ Save any input files for the WorkflowStep with appropriate suffix
        as determined by the WorkflowStep label. All files with <seed>.<ext>
        will be moved to <seed>.<ext>_<name>, for any <ext> inside the
        `input_exts` attribute. This is called after the WorkflowStep has
        finished, even if it does not succeed...

        Parameters:
            seed (str): seed for the workflow step.

        """
        if self.input_exts is not None:
            self._cache_files(seed, self.input_exts, 'in')

    def _cache_outputs(self, seed):
        """ Save any output files for the WorkflowStep with appropriate suffix
        as determined by the WorkflowStep label. All files with <seed>.<ext>
        will be moved to <seed>.<ext>_<name>, for any <ext> inside the
        `output_exts` attribute.

        Parameters:
            seed (str): seed for the workflow step.

        """
        if self.output_exts is not None:
            self._cache_files(seed, self.output_exts, 'out')

    def cache_files(self, seed):
        """ Wrapper for calling both _cache_inputs and _cache_outputs, without
        throwing any errors.
        """
        cwd = os.getcwd()
        if self.compute_dir is not None:
            os.chdir(self.compute_dir)

        self._cache_inputs(seed)
        self._cache_outputs(seed)

        if self.compute_dir is not None:
            os.chdir(cwd)

    def run_step(self, relaxer, calc_doc, seed):
        """ Run the workflow step.

        Parameters:
            relaxer (:obj:`matador.compute.ComputeTask`): the object that will be running the computation.
            calc_doc (dict): dictionary of structure and calculation
                parameters.
            seed (str): root seed for the calculation.

        Raises:
            RuntimeError: if any step fails.

        """
        try:
            LOG.info('WorkflowStep {} starting...'.format(self.name))
            self.success = self.function(relaxer, calc_doc, seed, **self.func_kwargs)
        except RuntimeError as exc:
            msg = 'WorkflowStep {} failed with error {}.'.format(self.name, exc)
            LOG.error(msg)
            self.success = False
            self.cache_files(seed)
            raise exc

        if self.success is None:
            LOG.info('WorkflowStep {} skipped, did you provide all the input files?'.format(self.name))
            return self.success

        if self.success:
            LOG.info('WorkflowStep {} completed successfully.'.format(self.name))
        else:
            LOG.warning('WorkflowStep {} was unsuccessful.'.format(self.name))

        self.cache_files(seed)

        return self.success
