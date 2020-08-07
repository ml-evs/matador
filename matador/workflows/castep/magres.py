# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the :class:`CastepMagresWorkflow`
class, which performs magres calculations with CASTEP in
multiple steps (only when necessary):

    1. Try to pre-relax structure (skipped if check file
       is already present).
    2. Perform an SCF with lower electronic tolerances.
    3. Calculate NMR properties, e.g. shielding and EFG.

"""


import os
import copy
import logging
from matador.workflows.workflows import Workflow
from matador.workflows.castep.common import castep_prerelax, castep_scf

LOG = logging.getLogger('run3')


def castep_full_phonon(computer, calc_doc, seed, **kwargs):
    """ Perform a "full" magres calculation on a system, i.e.
    first perform a relaxation, then do a high quality SCF
    and compute NMR properties in the same step.

    This function is a wrapper for the :class:`CastepMagresWorkflow` class.

    Parameters:
        computer (:obj:`matador.compute.ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): dictionary of structure and calculation
            parameters.
        seed (str): root seed for the calculation.

    Raises:
        RuntimeError: if any part of the calculation fails.

    Returns:
        bool: True if Workflow completed successfully, or False otherwise.

    """
    workflow = CastepMagresWorkflow(computer, calc_doc, seed, **kwargs)
    return workflow.success


class CastepMagresWorkflow(Workflow):
    """ Perform a "full" magres calculation on a system, i.e.
    first perform a relaxation in a standardised unit cell,
    then do a high quality SCF, then compute NMR properties.

    Attributes:
        computer (:obj:`ComputeTask`): the object that calls CASTEP.
        calc_doc (dict): the interim dictionary of structural and
            calculation parameters.
        seed (str): the root seed for the calculation.
        success (bool): the status of the Workflow: only set to True after
            post-processing method completes.

    """
    def preprocess(self):
        """ Decide which parts of the Workflow need to be performed,
        and set the appropriate CASTEP parameters.

        """
        # default todo
        todo = {'relax': True, 'scf': True, 'magres': True}
        # definition of steps and names
        steps = {'relax': castep_prerelax, 'magres': castep_magres_scf}

        exts = {
            'relax':
                {'input': ['.cell', '.param'], 'output': ['.castep', '-out.cell', '.*err']},
            'scf':
                {'input': ['.cell', '.param'], 'output': ['.castep']},
            'magres':
                {'input': ['.cell', '.param'], 'output': ['.castep', '.magres']}
        }

        # prepare to do pre-relax if there's no check file
        if os.path.isfile(self.seed + '.check'):
            todo['relax'] = False
            LOG.info('Restarting from {}.check, so not performing re-relaxation'.format(self.seed))

        for key in todo:
            if todo[key]:
                self.add_step(steps[key], key,
                              input_exts=exts[key].get('input'),
                              output_exts=exts[key].get('output'))


def castep_magres_scf(computer, calc_doc, seed):
    """ Run a singleshot SCF calculation with a high elec_energy_tol.

    Parameters:
        computer (:obj:`matador.compute.ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    Returns:
        bool: whether or not the SCF was successful.

    """

    return castep_scf(computer, calc_doc, seed, elec_energy_tol=1e-12)


def castep_magres(computer, calc_doc, seed):
    """ Runs a NMR properties calculation on top of a completed
    SCF calculation.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info('Performing CASTEP phonon dispersion calculation...')
    magres_doc = copy.deepcopy(calc_doc)
    magres_doc['task'] = 'magres'
    magres_doc['magres_task'] = 'nmr'
    magres_doc['continuation'] = 'default'

    required = []
    forbidden = []

    computer.validate_calc_doc(magres_doc, required, forbidden)

    return computer.run_castep_singleshot(magres_doc, seed, keep=True, intermediate=True)
