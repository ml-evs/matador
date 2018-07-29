# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the CastepSpectralWorkflow, which performs
spectral calculations with CASTEP in multiple steps:

1. Performs a singlepoint calculation (if check file is not found).
2. If `spectral_fine_kpoints_mp_spacing` keyword is found, interpolate
   wavefunction to DOS grid.
   - If an OptaDOS input file (.odi) with the root seedname
     is found, run OptaDOS on the resulting density of states.
3. If `spectral_fine_kpoints_path_spacing` keyword is present, create
   a bandstructure on the seekpath-generated path.

"""


import logging
import os
import copy
import shutil
import glob
from matador.workflows import Workflow


def castep_full_spectral(relaxer, calc_doc, seed):
    """ Perform a "full" spectral calculation on a system, i.e. first
    perform an SCF then interpolate to different kpoint paths/grids to
    form DOS and dispersions. Optionally use OptaDOS for post-processing
    of DOS.

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
    workflow = CastepSpectralWorkflow(relaxer, calc_doc, seed)
    return workflow.success


class CastepSpectralWorkflow(Workflow):
    """ Perform a "full" spectral calculation on a system, i.e. first
    perform an SCF then interpolate to different kpoint paths/grids to
    form DOS and dispersions. Optionally use OptaDOS for post-processing
    of DOS.

    Attributes:
        relaxer (:obj:`FullRelaxer`): the object that calls CASTEP.
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
        todo = {'scf': True, 'dos': True, 'dispersion': False}
        # definition of steps and names
        steps = {'scf': castep_spectral_scf,
                 'dos': castep_spectral_dos,
                 'dispersion': castep_spectral_dispersion}

        if (('spectral_kpoints_path' not in self.calc_doc and 'spectral_kpoints_list' not in self.calc_doc) and
                'spectral_kpoints_path_spacing' in self.calc_doc):
            todo['dispersion'] = True
        if 'spectral_kpoints_mp_spacing' in self.calc_doc:
            todo['dos'] = True

        # prepare to do scf if there's no check file
        if os.path.isfile(self.seed + '.check'):
            todo['scf'] = False
            logging.info('Restarting from {}.check, so not performing initial SCF'.format(self.seed))

        for key in todo:
            if todo[key]:
                self.add_step(steps[key], key)

        # always reduce cell to primitivve and standardise the cell so that any
        # post-processing performed after the fact will be consistent
        from matador.utils.cell_utils import cart2abc
        prim_doc, kpt_path = self.relaxer.get_seekpath_compliant_input(
            self.calc_doc, self.calc_doc.get('spectral_kpoints_path_spacing', 0.02))
        self.calc_doc.update(prim_doc)
        self.calc_doc['lattice_abc'] = cart2abc(self.calc_doc['lattice_cart'])

        if todo['dispersion']:
            self.calc_doc['spectral_kpoints_list'] = kpt_path

        logging.info('Preprocessing completed: run3 spectral options {}'.format(todo))


def castep_spectral_scf(relaxer, calc_doc, seed):
    """ Run a singleshot SCF calculation.

    Parameters:
        relaxer (:obj:`FullRelaxer`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    logging.info('Performing CASTEP spectral SCF...')
    scf_doc = copy.deepcopy(calc_doc)
    scf_doc['write_checkpoint'] = 'ALL'
    scf_doc['task'] = 'singlepoint'

    required = []
    forbidden = ['spectral_kpoints_list',
                 'spectral_kpoints_path',
                 'spectral_kpoints_mp_spacing',
                 'spectral_kpoints_path_spacing']

    relaxer.validate_calc_doc(scf_doc, required, forbidden)

    return relaxer.scf(scf_doc, seed, keep=True, intermediate=True)


def castep_spectral_dos(relaxer, calc_doc, seed):
    """ Runs a DOS interpolation on top of a completed SCF. If a single
    .odi file is found, run OptaDOS on the resulting DOS.

    Parameters:
        relaxer (:obj:`FullRelaxer`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    logging.info('Performing CASTEP spectral DOS calculation...')

    if len(glob.glob('*.odi')) == 1:
        odi_fname = glob.glob('*.odi')[0]
        optados = True
        logging.info('Detected single .odi file for OptaDOS post-processing.')
    else:
        optados = False
        logging.info('Detected {} .odi files; not performing OptaDOS post-processing.'.format(len(glob.glob('*.odi'))))

    dos_doc = copy.deepcopy(calc_doc)
    dos_doc['task'] = 'spectral'
    dos_doc['spectral_task'] = 'dos'
    if optados:
        dos_doc['pdos_calculate_weights'] = True

    required = ['spectral_kpoints_mp_spacing']
    forbidden = ['spectral_kpoints_list',
                 'spectral_kpoints_path',
                 'spectral_kpoints_path_spacing']

    relaxer.validate_calc_doc(dos_doc, required, forbidden)

    success = relaxer.scf(dos_doc, seed, keep=True, intermediate=True)

    if optados:
        shutil.copy(odi_fname, seed + '.odi')
    logging.info('Performing OptaDOS smearing calculation with parameters from {}'.format(odi_fname))
    _cache_executable = copy.deepcopy(relaxer.executable)
    relaxer.executable = 'optados'
    relaxer.parse_executable(seed)
    try:
        success = relaxer.run_generic(seed, intermediate=True)
    except Exception as exc:
        relaxer.executable = _cache_executable
        raise exc
    relaxer.executable = _cache_executable
    return success


def castep_spectral_dispersion(relaxer, calc_doc, seed):
    """ Runs a dispersion interpolation on top of a completed SCF calculation.

    Parameters:
        relaxer (:obj:`FullRelaxer`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    logging.info('Performing CASTEP phonon dispersion calculation...')
    disp_doc = copy.deepcopy(calc_doc)
    disp_doc['task'] = 'spectral'
    disp_doc['phonon_calculate_dos'] = False

    required = ['spectral_kpoints_list']
    forbidden = ['spectral_kpoints_mp_spacing',
                 'spectral_kpoints_path',
                 'spectral_kpoints_path_spacing']

    relaxer.validate_calc_doc(disp_doc, required, forbidden)

    return relaxer.scf(disp_doc, seed, keep=True, intermediate=True)
