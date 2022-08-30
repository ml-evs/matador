# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the :class:`CastepPhononWorkflow`
class, which performs phonon calculations with CASTEP in
multiple steps (only when necessary):

    1. Try to pre-relax structure (skipped if check file
           is already present).
    2. Calculate dynamical matrix.
    3. If ``phonon_fine_kpoint_mp_spacing`` keyword is found,
       interpolate dynamical matrix to form phonon DOS.
    4. If ``phonon_fine_kpoint_path_spacing`` key word is found,
       interpolate dynamical matrix to form phonon dispersion
       on the path given by seekpath.
    5. If ``task=thermodynamics``, perform a CASTEP thermodynamics
       calculation for the temperature-dependence of the free energy.

"""


import os
import copy
import logging
from matador.workflows.workflows import Workflow

LOG = logging.getLogger("run3")


def castep_full_phonon(computer, calc_doc, seed, **kwargs):
    """Perform a "full" phonon calculation on a system, i.e.
    first perform a relaxation in a standardised unit cell,
    then compute the dynamical matrix, then finally interpolate
    that dynamical matrix into dispersion curves and DOS. This function
    is a wrapper for the :class:`CastepPhononWorkflow` class.

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
    workflow = CastepPhononWorkflow(computer, calc_doc, seed, **kwargs)
    return workflow.success


class CastepPhononWorkflow(Workflow):
    """Perform a "full" phonon calculation on a system, i.e.
    first perform a relaxation in a standardised unit cell,
    then compute the dynamical matrix, then finally interpolate
    that dynamical matrix into dispersion curves and DOS.

    Attributes:
        computer (:obj:`ComputeTask`): the object that calls CASTEP.
        calc_doc (dict): the interim dictionary of structural and
            calculation parameters.
        seed (str): the root seed for the calculation.
        success (bool): the status of the Workflow: only set to True after
            post-processing method completes.

    """

    def preprocess(self):
        """Decide which parts of the Workflow need to be performed,
        and set the appropriate CASTEP parameters.

        """
        # default todo
        todo = {
            "relax": True,
            "dynmat": True,
            "vdos": False,
            "dispersion": False,
            "thermodynamics": False,
        }
        # definition of steps and names
        steps = {
            "relax": castep_phonon_prerelax,
            "dynmat": castep_phonon_dynmat,
            "vdos": castep_phonon_dos,
            "dispersion": castep_phonon_dispersion,
            "thermodynamics": castep_phonon_thermodynamics,
        }

        exts = {
            "relax": {
                "input": [".cell", ".param"],
                "output": [".castep", "-out.cell", ".*err"],
            },
            "dynmat": {"input": [".cell", ".param"], "output": [".castep", ".*err"]},
            "vdos": {
                "input": [".cell", ".param"],
                "output": [".castep", ".phonon", ".phonon_dos", ".*err"],
            },
            "dispersion": {
                "input": [".cell", ".param"],
                "output": [".castep", ".phonon", ".*err"],
            },
            "thermodynamics": {
                "input": [".cell", ".param"],
                "output": [".castep", ".*err"],
            },
        }

        if self.calc_doc.get("task").lower() in [
            "phonon",
            "thermodynamics",
            "phonon+efield",
        ]:
            if (
                "phonon_fine_kpoint_path" in self.calc_doc
                or "phonon_fine_kpoint_list" in self.calc_doc
                or "phonon_fine_kpoint_path_spacing" in self.calc_doc
            ):
                todo["dispersion"] = True
            if "phonon_fine_kpoint_mp_spacing" in self.calc_doc:
                todo["vdos"] = True
            if self.calc_doc["task"].lower() == "thermodynamics":
                todo["thermodynamics"] = True

        # prepare to do pre-relax if there's no check file
        if os.path.isfile(self.seed + ".check"):
            todo["relax"] = False
            LOG.info(
                "Restarting from {}.check, so not performing re-relaxation".format(
                    self.seed
                )
            )

        for key in todo:
            if todo[key]:
                self.add_step(
                    steps[key],
                    key,
                    input_exts=exts[key].get("input"),
                    output_exts=exts[key].get("output"),
                )

        # always standardise the cell so that any phonon calculation can have
        # post-processing performed after the fact, unless a path has been provided
        if (
            "phonon_fine_kpoint_list" not in self.calc_doc
            and "phonon_fine_kpoint_path" not in self.calc_doc
        ):
            from matador.utils.cell_utils import cart2abc

            prim_doc, kpt_path = self.computer.get_seekpath_compliant_input(
                self.calc_doc,
                self.calc_doc.get("phonon_fine_kpoint_path_spacing", 0.02),
            )
            self.calc_doc.update(prim_doc)
            self.calc_doc["lattice_abc"] = cart2abc(self.calc_doc["lattice_cart"])
            if todo["dispersion"]:
                self.calc_doc["phonon_fine_kpoint_list"] = kpt_path

        elif todo["dispersion"] and "phonon_fine_kpoint_path" in self.calc_doc:
            self._user_defined_kpt_path = True
            LOG.warning("Using user-defined k-point path for all structures.")
            self.calc_doc["phonon_fine_kpoint_spacing"] = self.calc_doc.get(
                "phonon_fine_kpoint_path_spacing", 0.05
            )

        # always shift phonon grid to include Gamma
        if "phonon_kpoint_mp_spacing" in self.calc_doc:
            from matador.utils.cell_utils import calc_mp_grid, shift_to_include_gamma

            grid = calc_mp_grid(
                self.calc_doc["lattice_cart"], self.calc_doc["phonon_kpoint_mp_spacing"]
            )
            offset = shift_to_include_gamma(grid)
            if offset != [0, 0, 0]:
                self.calc_doc["phonon_kpoint_mp_offset"] = offset
                LOG.debug("Set phonon MP grid offset to {}".format(offset))

        LOG.info("Preprocessing completed: run3 phonon options {}".format(todo))


def castep_phonon_prerelax(computer, calc_doc, seed):
    """Run a singleshot geometry optimisation before an SCF-style calculation.
    This is typically used to ensure phonon calculations start successfully.
    The phonon calculation will then be restarted from the .check file produced here.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    from matador.workflows.castep.common import castep_prerelax

    LOG.info("Performing CASTEP phonon pre-relax...")
    required = ["write_checkpoint"]
    forbidden = [
        "phonon_fine_kpoint_list",
        "phonon_fine_kpoint_path",
        "phonon_fine_kpoint_mp_spacing",
        "phonon_fine_kpoint_path_spacing",
    ]

    return castep_prerelax(
        computer, calc_doc, seed, required_keys=required, forbidden_keys=forbidden
    )


def castep_phonon_dynmat(computer, calc_doc, seed):
    """Runs a singleshot phonon dynmat calculation, with no "fine_method" interpolation.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info("Performing CASTEP dynmat calculation...")
    dynmat_doc = copy.deepcopy(calc_doc)
    dynmat_doc["write_checkpoint"] = "ALL"
    if calc_doc["task"].lower() == "phonon+efield":
        dynmat_doc["task"] = "phonon+efield"
    else:
        dynmat_doc["task"] = "phonon"

    dynmat_doc["continuation"] = "default"

    required = ["continuation", "write_checkpoint"]
    forbidden = [
        "phonon_fine_kpoint_list",
        "phonon_fine_kpoint_path",
        "phonon_fine_kpoint_mp_spacing",
        "phonon_fine_kpoint_path_spacing",
    ]

    computer.validate_calc_doc(dynmat_doc, required, forbidden)
    return computer.run_castep_singleshot(
        dynmat_doc, seed, keep=True, intermediate=True
    )


def castep_phonon_dos(computer, calc_doc, seed):
    """Runs a DOS interpolation on top of a completed
    phonon calculation.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info("Performing CASTEP phonon DOS calculation...")
    dos_doc = copy.deepcopy(calc_doc)
    dos_doc["task"] = "phonon"
    dos_doc["phonon_calculate_dos"] = True
    dos_doc["continuation"] = "default"

    required = ["phonon_fine_kpoint_mp_spacing"]
    forbidden = [
        "phonon_fine_kpoint_list",
        "phonon_fine_kpoint_path",
        "phonon_fine_kpoint_path_spacing",
    ]

    computer.validate_calc_doc(dos_doc, required, forbidden)

    return computer.run_castep_singleshot(dos_doc, seed, keep=True, intermediate=True)


def castep_phonon_dispersion(computer, calc_doc, seed):
    """Runs a dispersion interpolation on top of a completed
    phonon calculation.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info("Performing CASTEP phonon dispersion calculation...")
    disp_doc = copy.deepcopy(calc_doc)
    disp_doc["task"] = "phonon"
    disp_doc["phonon_calculate_dos"] = False
    disp_doc["continuation"] = "default"

    required = []
    forbidden = ["phonon_fine_kpoint_mp_spacing"]

    computer.validate_calc_doc(disp_doc, required, forbidden)

    return computer.run_castep_singleshot(disp_doc, seed, keep=True, intermediate=True)


def castep_phonon_thermodynamics(computer, calc_doc, seed):
    """Runs a "thermodynamics" interpolation on top of a completed
    phonon calculation, using the phonon_fine_kpoint_mp_grid.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info("Performing CASTEP thermodynamics calculation...")
    thermo_doc = copy.deepcopy(calc_doc)
    thermo_doc["continuation"] = "default"
    thermo_doc["task"] = "thermodynamics"
    thermo_doc["phonon_calculate_dos"] = False

    required = ["phonon_fine_kpoint_mp_spacing"]
    forbidden = [
        "phonon_fine_kpoint_list",
        "phonon_fine_kpoint_path",
        "phonon_fine_kpoint_path_spacing",
    ]

    computer.validate_calc_doc(thermo_doc, required, forbidden)

    return computer.run_castep_singleshot(
        thermo_doc, seed, keep=True, intermediate=True
    )
