# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the CastepSpectralWorkflow, which performs
spectral calculations with CASTEP in multiple steps:

1. Performs a singlepoint calculation (if check file is not found).
2. If `spectral_kpoints_mp_spacing` keyword is found, interpolate
   wavefunction to DOS grid.
   - If an OptaDOS input file (.odi) with the root seedname
     is found, run OptaDOS on the resulting density of states.
3. If `spectral_kpoints_path_spacing` keyword is present, create
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
        todo = {'scf': True, 'dos': False, 'dispersion': False}
        # definition of steps and names
        steps = {'scf': castep_spectral_scf,
                 'dos': castep_spectral_dos,
                 'dispersion': castep_spectral_dispersion}

        if os.path.isfile(self.seed + '.check'):
            logging.info('Found {}.check, so skipping initial SCF.'.format(self.seed))
            todo['scf'] = False

        if (
                (
                    'spectral_kpoints_path' not in self.calc_doc and
                    'spectral_kpoints_list' not in self.calc_doc
                ) and
                (
                    'spectral_kpoints_path_spacing' in self.calc_doc or
                    self.calc_doc.get('spectral_task').lower() == 'bandstructure'
                )
        ):
            todo['dispersion'] = True
        if (
                'spectral_kpoints_mp_spacing' in self.calc_doc or
                self.calc_doc.get('spectral_task').lower() == 'dos'
        ):
            todo['dos'] = True

        for key in todo:
            if todo[key]:
                self.add_step(steps[key], key)

        # always reduce cell to primitive and standardise the cell so that any
        # post-processing performed after the fact will be consistent
        from matador.utils.cell_utils import cart2abc
        self.calc_doc['spectral_kpoints_mp_spacing'] = self.calc_doc.get('spectral_kpoints_mp_spacing', 0.05)
        self.calc_doc['spectral_kpoints_path_spacing'] = self.calc_doc.get('spectral_kpoints_path_spacing', 0.05)
        prim_doc, kpt_path = self.relaxer.get_seekpath_compliant_input(
            self.calc_doc, self.calc_doc['spectral_kpoints_path_spacing'])
        self.calc_doc.update(prim_doc)
        self.calc_doc['lattice_abc'] = cart2abc(self.calc_doc['lattice_cart'])
        self.calc_doc['continuation'] = 'default'

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
    if 'spectral_task' in scf_doc:
        del scf_doc['spectral_task']

    required = []
    forbidden = ['spectral_task',
                 'spectral_kpoints_list',
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

    dos_doc = copy.deepcopy(calc_doc)
    dos_doc['task'] = 'spectral'
    dos_doc['spectral_task'] = 'dos'
    # disable checkpointing for BS/DOS by default, leaving just SCF
    dos_doc['write_checkpoint'] = 'none'
    dos_doc['pdos_calculate_weights'] = True
    if os.path.isfile('{}.bands'.format(seed)):
        shutil.copy2('{}.bands'.format(seed), '{}.bands_bak'.format(seed))

    required = ['spectral_kpoints_mp_spacing']
    forbidden = ['spectral_kpoints_list',
                 'spectral_kpoints_path',
                 'spectral_kpoints_path_spacing']

    relaxer.validate_calc_doc(dos_doc, required, forbidden)

    success = relaxer.scf(dos_doc, seed, keep=True, intermediate=True)

    shutil.move('{}.bands'.format(seed), '{}.bands_dos'.format(seed))
    if os.path.isfile('{}.bands_bak'.format(seed)):
        shutil.move('{}.bands_bak'.format(seed), '{}.bands'.format(seed))
    shutil.copy2('{}.cell'.format(seed), '{}.cell_dos'.format(seed))
    shutil.copy2('{}.param'.format(seed), '{}.param_dos'.format(seed))

    from matador.scrapers import arbitrary2dict
    from matador.export import doc2arbitrary

    odi_fname = None
    if os.path.isfile(seed + '.odi'):
        os.remove(seed + '.odi')

    if glob.glob('*.odi'):
        # dodginess: choose the odi file with the shortest name...
        shortest_fname = None
        for fname in glob.glob('*.odi'):
            if shortest_fname is None or len(fname) < len(shortest_fname):
                shortest_fname = fname
        odi_fname = shortest_fname

    if odi_fname is not None:
        odi_dict, _ = arbitrary2dict(odi_fname)
        _cache_executable = copy.deepcopy(relaxer.executable)
        _cache_core = copy.deepcopy(relaxer.ncores)
        relaxer.ncores = 1
        relaxer.executable = 'optados'
        relaxer.parse_executable(seed)

        # if pdos keyword is present, try to run a pDOS
        if 'pdos' in odi_dict:
            odi_dict['task'] = 'pdos'
            if 'pdispersion' in odi_dict:
                del odi_dict['pdispersion']

            logging.info('Performing OptaDOS pDOS calculation with parameters from {}'.format(odi_fname))
            doc2arbitrary(odi_dict, seed + '.odi', overwrite=True)
            shutil.copy2('{}.odi'.format(seed), '{}.odi_pdos'.format(seed))
            if os.path.isfile('{}.pdos_bin'.format(seed)):
                shutil.copy2('{}.pdos_bin'.format(seed), '{}.pdos_bin_pdos_bak'.format(seed))

            try:
                success = relaxer.run_generic(seed, intermediate=True)
            except Exception as exc:
                relaxer.executable = _cache_executable
                logging.warning('Failed to call optados, with error: {}'.format(exc))

        # if broadening keyword is present, try to run a normal DOS
        if 'broadening' in odi_dict:
            odi_dict['task'] = 'dos'
            if 'pdos' in odi_dict:
                del odi_dict['pdos']
            if 'pdispersion' in odi_dict:
                del odi_dict['pdispersion']

            logging.info('Performing OptaDOS DOS calculation with parameters from {}'.format(odi_fname))
            doc2arbitrary(odi_dict, seed + '.odi', overwrite=True)
            shutil.copy2('{}.odi'.format(seed), '{}.odi_dos'.format(seed))

            try:
                success = relaxer.run_generic(seed, intermediate=True)
            except Exception as exc:
                relaxer.executable = _cache_executable
                relaxer.ncores = _cache_core
                logging.warning('Failed to call optados, with error: {}'.format(exc))

            relaxer.executable = _cache_executable
            relaxer.ncores = _cache_core

    return success


def castep_spectral_dispersion(relaxer, calc_doc, seed):
    """ Runs a dispersion interpolation on top of a completed SCF calculation,
    optionally running orbitals2bands and OptaDOS projected dispersion.

    Parameters:
        relaxer (:obj:`FullRelaxer`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    logging.info('Performing CASTEP spectral dispersion calculation...')
    disp_doc = copy.deepcopy(calc_doc)
    disp_doc['task'] = 'spectral'
    disp_doc['spectral_task'] = 'bandstructure'

    # disable checkpointing for BS/DOS by default, leaving just SCF
    disp_doc['write_checkpoint'] = 'none'
    disp_doc['pdos_calculate_weights'] = True
    disp_doc['continuation'] = 'default'

    required = ['spectral_kpoints_list']
    forbidden = ['spectral_kpoints_mp_spacing',
                 'spectral_kpoints_path',
                 'spectral_kpoints_path_spacing']

    relaxer.validate_calc_doc(disp_doc, required, forbidden)
    success = relaxer.scf(disp_doc, seed, keep=True, intermediate=True)

    if disp_doc.get('write_orbitals'):
        logging.info('Planning to call orbitals2bands...')
        if os.path.isfile('{}.bands'.format(seed)):
            shutil.copy2('{}.bands'.format(seed), '{}.bands_orig'.format(seed))
        shutil.copy2('{}.cell'.format(seed), '{}.cell_bs'.format(seed))
        shutil.copy2('{}.param'.format(seed), '{}.param_bs'.format(seed))

        _cache_executable = copy.deepcopy(relaxer.executable)
        _cache_core = copy.deepcopy(relaxer.ncores)
        relaxer.ncores = 1
        relaxer.executable = 'orbitals2bands'
        try:
            success = relaxer.run_generic(seed, intermediate=True)
        except Exception as exc:
            relaxer.executable = _cache_executable
            relaxer.ncores = _cache_core
            logging.warning('Failed to call orbitals2bands, with error: {}'.format(exc))

        relaxer.ncores = _cache_core
        relaxer.executable = _cache_executable

    from matador.scrapers import arbitrary2dict
    from matador.export import doc2arbitrary

    odi_fname = None
    if os.path.isfile(seed + '.odi'):
        os.remove(seed + '.odi')

    if glob.glob('*.odi'):
        # dodginess: choose the odi file with the shortest name...
        shortest_fname = None
        for fname in glob.glob('*.odi'):
            if shortest_fname is None or len(fname) < len(shortest_fname):
                shortest_fname = fname
        odi_fname = shortest_fname

    if odi_fname is not None:
        odi_dict, _ = arbitrary2dict(odi_fname)

        # if pdispersion keyword is present, try to run a pdis
        if 'pdispersion' in odi_dict:
            odi_dict['task'] = 'pdispersion'
            if 'pdos' in odi_dict:
                del odi_dict['pdos']
            logging.info('Performing OptaDOS pDIS calculation with parameters from {}'.format(odi_fname))

            doc2arbitrary(odi_dict, seed + '.odi', overwrite=True)
            if os.path.isfile('{}.pdos_bin'.format(seed)):
                shutil.copy2('{}.pdos_bin'.format(seed), '{}.pdos_bin_pdis_bak'.format(seed))

            _cache_executable = copy.deepcopy(relaxer.executable)
            _cache_core = copy.deepcopy(relaxer.ncores)

            relaxer.ncores = 1
            relaxer.executable = 'optados'
            shutil.copy2('{}.odi'.format(seed), '{}.odi_bs'.format(seed))
            try:
                success = relaxer.run_generic(seed, intermediate=True)
            except Exception as exc:
                relaxer.executable = _cache_executable
                relaxer.ncores = _cache_core
                logging.warning('Failed to call optados, with error: {}'.format(exc))


            relaxer.ncores = _cache_core
            relaxer.executable = _cache_executable

    return success
