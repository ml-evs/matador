# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the CastepElasticWorkflow, which performs
elastic tensor calculations with CASTEP. So far, only bulk modulus
calculations (i.e. varying volume from equilibrium) have been implemented.

"""


import logging
import copy
import numpy as np
from matador.workflows import Workflow


def castep_elastic(relaxer, calc_doc, seed):
    """ Perform a calculation of the elastic tensor on a system.
    Currently only calculation of the bulk modulus is implemented.

    Parameters:
        relaxer (:obj:`FullRelaxer`): the object that will be calling CASTEP.
        calc_doc (dict): dictionary of structure and calculation
            parameters.
        seed (str): root seed for the calculation.

    Raises:
        RuntimeError: if any part of the calculation fails.

    Returns:
        bool: True if Workflow completed successfully, or False otherwise.

    """
    workflow = CastepElasticWorkflow(relaxer, calc_doc, seed, bulk_modulus_only=True)
    return workflow.success


class CastepElasticWorkflow(Workflow):
    """ Perform a "full" spectral calculation on a system, i.e. first
    perform an SCF then interpolate to different kpoint paths/grids to
    form DOS and dispersions. Optionally use OptaDOS for post-processing
    of DOS.

    Attributes:
        relaxer (:obj:`FullRelaxer`): the object that calls CASTEP.
        calc_doc (dict): the interim dictionary of structural and
            calculation parameters.
        seed (str): the root seed for the calculation.
        volume_rescale (list): a list of floats containing the amounts
            by which to rescale volumes when performing bulk modulus calculation.
        success (bool): the status of the Workflow: only set to True after
            post-processing method completes.

    Keyword arguments:
        bulk_modulus_only (bool): only calculate the bulk modulus.

    """
    def preprocess(self):
        """ Decide which parts of the Workflow need to be performed,
        and set the appropriate CASTEP parameters.

        """
        num_volumes = 9
        self.volume_rescale = np.cbrt(np.geomspace(0.7, 1.2, num=num_volumes, endpoint=True)).tolist()
        self.volume_rescale.append(1.0)
        self._completed_volumes = []
        logging.info('Preprocessing completed: run3 bulk modulus calculation options {}'
                     .format(self.volume_rescale))

        for volume in self.volume_rescale:
            self.add_step(castep_rescaled_volume_scf, 'rescaled_volume_scf', rescale=volume)


def castep_rescaled_volume_scf(relaxer, calc_doc, seed, rescale=1):
    """ Run a singleshot SCF calculation.

    Parameters:
        relaxer (:obj:`FullRelaxer`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    Keyword arguments:
        rescale (float): amount to rescale volume by.

    """
    assert rescale > 0
    logging.info('Performing CASTEP SCF on volume rescaled by {}.'.format(rescale**3))
    scf_doc = copy.deepcopy(calc_doc)
    for i in range(3):
        for k in range(3):
            scf_doc['lattice_cart'][i][k] *= rescale
    scf_doc['task'] = 'singlepoint'

    return relaxer.scf(scf_doc, seed, keep=True, intermediate=True)
