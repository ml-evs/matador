# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the CastepSpectralWorkflow, which performs
spectral calculations with CASTEP in multiple steps:

1. Performs a singlepoint calculation (if check file is not found).
2. If the `spectral_kpoints_mp_spacing` keyword is found, interpolate
   wavefunction to DOS grid.
   - If an OptaDOS input file (.odi) with the root seedname
     is found, run OptaDOS on the resulting density of states.
3. If `spectral_kpoints_path_spacing` keyword is present, create
   a bandstructure on the SeeK-path generated path.

"""


import os
import copy
import glob
import logging
from matador.workflows import Workflow
from matador.scrapers import arbitrary2dict

LOG = logging.getLogger('run3')


def castep_full_spectral(relaxer, calc_doc, seed, **kwargs):
    """ Perform a "full" spectral calculation on a system, i.e. first
    perform an SCF then interpolate to different kpoint paths/grids to
    form DOS and dispersions. Optionally use OptaDOS for post-processing
    of DOS.

    Parameters:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be calling CASTEP.
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
        relaxer (:obj:`matador.compute.ComputeTask`): the object that calls CASTEP.
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
        todo = {'scf': True, 'dos': False, 'pdos': False, 'broadening': False, 'dispersion': False, 'pdis': False}
        # definition of steps and names
        steps = {'scf': castep_spectral_scf,
                 'dos': castep_spectral_dos,
                 'pdos': optados_pdos,
                 'broadening': optados_dos_broadening,
                 'dispersion': castep_spectral_dispersion,
                 'pdis': optados_pdispersion}

        exts = {
            'scf': {
                'input': ['.cell', '.param'],
                'output': ['.castep', '.*err', '-out.cell']
            },
            'dos': {
                'input': ['.cell', '.param'],
                'output': ['.castep', '.bands', '.pdos_bin', '.dome_bin', '.*err', '-out.cell']
            },
            'dispersion': {
                'input': ['.cell', '.param'],
                'output': ['.castep', '.bands', '.pdos_bin', '.dome_bin', '.*err', '-out.cell']
            },
            'pdis': {
                'input': ['.odi', '.pdos_bin'],
                'output': ['.odo', '.*err']
            },
            'pdos': {
                'input': ['.odi', '.pdos_bin', '.dome_bin'],
                'output': ['.odo', '.*err']
            },
            'broadening': {
                'input': ['.odi', '.pdos_bin', '.dome_bin'],
                'output': ['.odo', '.*err']
            }
        }

        if os.path.isfile(self.seed + '.check'):
            LOG.info('Found {}.check, so skipping initial SCF.'.format(self.seed))
            todo['scf'] = False

        if (
                (
                    'spectral_kpoints_path' in self.calc_doc or
                    'spectral_kpoints_list' in self.calc_doc or
                    'spectral_kpoints_path_spacing' in self.calc_doc or
                    self.calc_doc.get('spectral_task', '').lower() == 'bandstructure'
                )
        ):
            todo['dispersion'] = not os.path.isfile(self.seed + '.bands_dispersion')

        if ('spectral_kpoints_mp_spacing' in self.calc_doc or
                self.calc_doc.get('spectral_task', '').lower() == 'dos'):
            todo['dos'] = not os.path.isfile(self.seed + '.bands_dos')

        odi_fname = _get_optados_fname(self.seed)
        if odi_fname is not None:
            odi_dict, _ = arbitrary2dict(odi_fname)
            if todo['dispersion']:
                todo['pdis'] = 'pdispersion' in odi_dict
            if todo['dos']:
                todo['broadening'] = 'broadening' in odi_dict
                todo['pdos'] = 'pdos' in odi_dict

        for key in todo:
            if todo[key]:
                self.add_step(steps[key], key,
                              input_exts=exts[key].get('input'),
                              output_exts=exts[key].get('output'))

        if self.relaxer.run3_settings.get('run3_settings') is not None:
            settings = self.relaxer.kwargs.get('run3_settings')
            # check that relaxer.exec was not overriden at cmd-line, then check settings file
            if settings.get('castep_executable') is not None and self.relaxer.executable == 'castep':
                self.castep_executable = settings.get('castep_executable', 'castep')
                self.relaxer.executable = self.castep_executable
            if settings.get('optados_executable') is not None:
                self.optados_executable = settings.get('optados_executable', 'optados')
                self.relaxer.optados_executable = self.optados_executable

        # if not using a user-requested path, use seekpath and spglib
        # to reduce to primitive and use consistent path
        if 'spectral_kpoints_list' not in self.calc_doc and 'spectral_kpoints_path' not in self.calc_doc:
            from matador.utils.cell_utils import cart2abc
            prim_doc, kpt_path = self.relaxer.get_seekpath_compliant_input(
                self.calc_doc, self.calc_doc.get('spectral_kpoints_path_spacing', 0.05))
            self.calc_doc.update(prim_doc)
            self.calc_doc['lattice_abc'] = cart2abc(self.calc_doc['lattice_cart'])
            if todo['dispersion']:
                self.calc_doc['spectral_kpoints_list'] = kpt_path
        elif todo['dispersion'] and 'spectral_kpoints_path' in self.calc_doc:
            self._user_defined_kpt_path = True
            LOG.warning('Using user-defined k-point path for all structures.')
            self.calc_doc['spectral_kpoints_path_spacing'] = self.calc_doc.get('spectral_kpoints_path_spacing', 0.05)

        if todo['dos']:
            self.calc_doc['spectral_kpoints_mp_spacing'] = self.calc_doc.get('spectral_kpoints_mp_spacing', 0.05)

        # always use continuation
        self.calc_doc['continuation'] = 'default'

        LOG.info('Preprocessing completed: run3 spectral options {}'.format(todo))


def castep_spectral_scf(relaxer, calc_doc, seed):
    """ Run a singleshot SCF calculation.

    Parameters:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info('Performing CASTEP spectral SCF...')
    scf_doc = copy.deepcopy(calc_doc)
    scf_doc['write_checkpoint'] = 'ALL'
    scf_doc['task'] = 'singlepoint'
    scf_doc['write_cell_structure'] = True
    if 'spectral_task' in scf_doc:
        del scf_doc['spectral_task']

    required = []
    forbidden = ['spectral_task',
                 'spectral_kpoints_list',
                 'spectral_kpoints_path',
                 'spectral_kpoints_mp_spacing',
                 'spectral_kpoints_path_spacing']

    relaxer.validate_calc_doc(scf_doc, required, forbidden)

    return relaxer.run_castep_singleshot(scf_doc, seed, keep=True, intermediate=True)


def castep_spectral_dos(relaxer, calc_doc, seed):
    """ Runs a DOS interpolation on top of a completed SCF. If a single
    .odi file is found, run OptaDOS on the resulting DOS.

    Parameters:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info('Performing CASTEP spectral DOS calculation...')

    dos_doc = copy.deepcopy(calc_doc)
    dos_doc['task'] = 'spectral'
    dos_doc['spectral_task'] = 'dos'
    # disable checkpointing for BS/DOS by default, leaving just SCF
    dos_doc['write_checkpoint'] = 'none'
    dos_doc['write_cell_structure'] = True
    dos_doc['pdos_calculate_weights'] = True

    required = ['spectral_kpoints_mp_spacing']
    forbidden = ['spectral_kpoints_list',
                 'spectral_kpoints_path',
                 'spectral_kpoints_path_spacing']

    relaxer.validate_calc_doc(dos_doc, required, forbidden)
    success = relaxer.run_castep_singleshot(dos_doc, seed, keep=True, intermediate=True)

    return success


def castep_spectral_dispersion(relaxer, calc_doc, seed):
    """ Runs a dispersion interpolation on top of a completed SCF calculation,
    optionally running orbitals2bands and OptaDOS projected dispersion.

    Parameters:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info('Performing CASTEP spectral dispersion calculation...')
    disp_doc = copy.deepcopy(calc_doc)
    disp_doc['task'] = 'spectral'
    disp_doc['spectral_task'] = 'bandstructure'

    # disable checkpointing for BS/DOS by default, leaving just SCF
    disp_doc['write_checkpoint'] = 'none'
    disp_doc['pdos_calculate_weights'] = True
    disp_doc['write_cell_structure'] = True
    disp_doc['continuation'] = 'default'

    required = []
    forbidden = ['spectral_kpoints_mp_spacing']

    relaxer.validate_calc_doc(disp_doc, required, forbidden)
    success = relaxer.run_castep_singleshot(disp_doc, seed, keep=True, intermediate=True)

    if disp_doc.get('write_orbitals'):
        LOG.info('Planning to call orbitals2bands...')

        _cache_executable = copy.deepcopy(relaxer.executable)
        _cache_core = copy.deepcopy(relaxer.ncores)
        relaxer.ncores = 1
        relaxer.executable = 'orbitals2bands'
        try:
            success = relaxer.run_generic(intermediate=True, mv_bad_on_failure=False)
        except Exception as exc:
            relaxer.executable = _cache_executable
            relaxer.ncores = _cache_core
            LOG.warning('Failed to call orbitals2bands, with error: {}'.format(exc))

        relaxer.ncores = _cache_core
        relaxer.executable = _cache_executable

    return success


def optados_pdos(relaxer, _, seed):
    """ Run an OptaDOS projected-DOS.

    Parameters:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be calling OptaDOS.
        _ : second parameter is required but ignored.
        seed (str): root filename of structure.

    """

    odi_fname = _get_optados_fname(seed)
    if odi_fname is not None:
        odi_dict, _ = arbitrary2dict(odi_fname)

        odi_dict['task'] = 'pdos'
        if 'pdispersion' in odi_dict:
            del odi_dict['pdispersion']

        LOG.info('Performing OptaDOS pDOS calculation with parameters from {}'.format(odi_fname))
        success = _run_optados(relaxer, odi_dict, seed, suffix='dos')
        return success

    return None


def optados_dos_broadening(relaxer, _, seed):
    """ Run an OptaDOS total DOS broadening.

    Parameters:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be calling OptaDOS.
        _ : second parameter is required but ignored.
        seed (str): root filename of structure.

    """

    odi_fname = _get_optados_fname(seed)
    if odi_fname is not None:
        odi_dict, _ = arbitrary2dict(odi_fname)
        # if broadening keyword is present, try to run a normal DOS
        odi_dict['task'] = 'dos'
        if 'pdos' in odi_dict:
            del odi_dict['pdos']
        if 'pdispersion' in odi_dict:
            del odi_dict['pdispersion']

        LOG.info('Performing OptaDOS DOS broadening with parameters from {}'.format(odi_fname))
        return _run_optados(relaxer, odi_dict, seed, suffix='dos')

    return None


def optados_pdispersion(relaxer, _, seed):
    """ Runs an OptaDOS projected dispersion calculation.

    Parameters:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be calling OptaDOS.
        _ : second parameter is required but ignored.
        seed (str): root filename of structure.

    """
    odi_fname = _get_optados_fname(seed)
    if odi_fname is not None:
        odi_dict, _ = arbitrary2dict(odi_fname)
        # if pdispersion keyword is present, try to run a pdis
        odi_dict['task'] = 'pdispersion'
        if 'pdos' in odi_dict:
            del odi_dict['pdos']

        LOG.info('Performing OptaDOS pDIS calculation with parameters from {}'.format(odi_fname))
        return _run_optados(relaxer, odi_dict, seed, suffix='dispersion')

    return None


def _run_optados(relaxer, odi_dict, seed, suffix=None):
    """ Run OptaDOS with given relaxer object, parameters and seed, adjusting
    the number of cores and the executable to call, then restoring them after.

    Parameters:
        relaxer (:obj:`matador.compute.ComputeTask`): the object that will be calling OptaDOS.
        odi_dict (dict): the OptaDOS parameters to write to file.
        seed (str): root filename of structure.

    Keyword arguments:
        suffix (str): either 'dos' or 'dispersion' for backup/restore.
        executable (str): optados executable to call.

    """

    from matador.export import doc2arbitrary
    odi_path = '{}.odi'.format(seed)
    if relaxer.compute_dir is not None:
        odi_path = relaxer.compute_dir + '/' + odi_path
    doc2arbitrary(odi_dict, odi_path, overwrite=True)

    if suffix is not None:
        _get_correct_files_for_optados(seed, suffix=suffix)

    _cache_executable = copy.deepcopy(relaxer.executable)
    _cache_core = copy.deepcopy(relaxer.ncores)
    _cache_nodes = copy.deepcopy(relaxer.nnodes)
    relaxer.ncores = 1
    relaxer.nnodes = 1
    relaxer.executable = relaxer.optados_executable
    success = False

    try:
        success = relaxer.run_generic(intermediate=True, mv_bad_on_failure=False)
    except Exception as exc:
        LOG.warning('Failed to call optados with error: {}'.format(exc))

    relaxer.ncores = _cache_core
    relaxer.nnodes = _cache_nodes
    relaxer.executable = _cache_executable
    if suffix is not None:
        _get_correct_files_for_optados(seed, suffix='bak')

    return success


def _get_optados_fname(seed):
    """ Get the most likely OptaDOS input file name, which here means either
    only existing one, or the shortest one.

    Parameters:
        seed (str): seedname for the calculation.

    Returns:
        str: the OptaDOS input filename.

    """
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

    return odi_fname


def _get_correct_files_for_optados(seed, suffix=None):
    """ If e.g. dispersion and DOS calculations were run previously, but
    it is unclear which exists as the current <seed>.<ext>, try to copy
    the old <seed>.<ext>_dispersion/dos files to the correct place.

    Parameters:
        seed (str): seedname for the calculation.

    Keyword arguments:
        suffix (str): either 'dos' or 'dispersion', or 'bak' to get old files back.

    """
    import shutil
    LOG.debug('Getting files for OptaDOS: {} {}'.format(seed, suffix))
    if suffix is None:
        return
    files_to_cache = ['.pdos_bin', '.dome_bin', '-out.cell', '.bands', '.cell', '.param']
    for ext in files_to_cache:
        old_file = '{}{}_{}'.format(seed, ext, suffix)
        current_file = '{}{}'.format(seed, ext)
        if os.path.isfile(current_file):
            if suffix != 'bak':
                backup_file = '{}{}_{}'.format(seed, ext, 'bak')
                shutil.copy2(current_file, backup_file)
            os.remove(current_file)
        if os.path.isfile(old_file):
            LOG.debug('Copying {} to {}'.format(old_file, current_file))
            shutil.copy2(old_file, current_file)
            if suffix == 'bak':
                os.remove(old_file)
