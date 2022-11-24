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


import copy
import logging
import os
from functools import partial

from matador.scrapers import arbitrary2dict
from matador.workflows.castep.common import castep_prerelax, castep_scf
from matador.workflows.castep.spectral import (
    _get_optados_fname,
    castep_spectral_dos,
    optados_dos_broadening,
    optados_pdos,
)
from matador.workflows.workflows import Workflow

LOG = logging.getLogger("run3")

__all__ = "castep_full_magres"


def castep_full_magres(computer, calc_doc, seed, final_elec_energy_tol=1e-11, **kwargs):
    """Perform a "full" magres calculation on a system, i.e.
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
    workflow = CastepMagresWorkflow(
        computer, calc_doc, seed, final_elec_energy_tol=final_elec_energy_tol, **kwargs
    )

    return workflow.success


class CastepMagresWorkflow(Workflow):
    """Perform a "full" magres calculation on a system, i.e.
    first perform a relaxation in a standardised unit cell,
    then do a high quality SCF, then compute NMR properties.

    Attributes:
        computer (:obj:`ComputeTask`): the object that calls CASTEP.
        calc_doc (dict): the interim dictionary of structural and
            calculation parameters.
        seed (str): the root seed for the calculation.
        success (bool): the status of the Workflow: only set to True after
            post-processing method completes.
        final_elec_energy_tol (float): the electronic energy tolerance to use
            in the high-quality SCF calculation.

    """

    def preprocess(self):
        """Decide which parts of the Workflow need to be performed,
        and set the appropriate CASTEP parameters.

        """

        self.final_elec_energy_tol = self.workflow_params.get(
            "final_elec_energy_tol", 1e-11
        )
        # default todo
        todo = {
            "relax": True,
            "scf": True,
            "dos": True,
            "pdos": True,
            "broadening": True,
            "magres": True,
        }
        # definition of steps and names
        steps = {
            "relax": castep_prerelax,
            "scf": partial(
                castep_magres_scf, elec_energy_tol=self.final_elec_energy_tol
            ),
            "dos": castep_spectral_dos,
            "pdos": optados_pdos,
            "broadening": optados_dos_broadening,
            "magres": partial(
                castep_magres, elec_energy_tol=self.final_elec_energy_tol
            ),
        }

        exts = {
            "relax": {
                "input": [".cell", ".param"],
                "output": [".castep", "-out.cell", ".*err"],
            },
            "scf": {"input": [".cell", ".param"], "output": [".castep", ".bands"]},
            "magres": {"input": [".cell", ".param"], "output": [".castep", ".magres"]},
            "dos": {
                "input": [".cell", ".param"],
                "output": [
                    ".castep",
                    ".bands",
                    ".pdos_bin",
                    ".dome_bin",
                    ".*err",
                    "-out.cell",
                ],
            },
            "pdos": {
                "input": [".odi", ".pdos_bin", ".dome_bin"],
                "output": [".odo", ".*err"],
            },
            "broadening": {
                "input": [".odi", ".pdos_bin", ".dome_bin"],
                "output": [".odo", ".*err"],
            },
        }

        odi_fname = _get_optados_fname(self.seed)
        if odi_fname is not None:
            odi_dict, _ = arbitrary2dict(odi_fname)
            if todo["dos"]:
                todo["broadening"] = "broadening" in odi_dict
                todo["pdos"] = "pdos" in odi_dict
        else:
            todo["dos"] = False
            todo["pdos"] = False
            todo["broadening"] = False

        # prepare to do pre-relax if there's no check file
        if os.path.isfile(self.seed + ".check"):
            todo["scf"] = True
            todo["relax"] = False
            LOG.info(
                "Restarting from {}.check, so not performing re-relaxation".format(
                    self.seed
                )
            )

        # If geom force tol is not set, do not perform a relaxation
        if self.calc_doc.get("geom_force_tol") is None:
            todo["relax"] = False

        for key in todo:
            if todo[key]:
                self.add_step(
                    steps[key],
                    key,
                    input_exts=exts[key].get("input"),
                    output_exts=exts[key].get("output"),
                )


def castep_magres_scf(computer, calc_doc, seed, elec_energy_tol=1e-11):
    """Run a singleshot SCF calculation with a high elec_energy_tol.

    Parameters:
        computer (:obj:`matador.compute.ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    Returns:
        bool: whether or not the SCF was successful.

    """
    calc_doc["write_checkpoint"] = "ALL"
    calc_doc["continuation"] = "default"

    required = ["write_checkpoint", "continuation"]

    return castep_scf(
        computer,
        calc_doc,
        seed,
        elec_energy_tol=elec_energy_tol,
        required_keys=required,
    )


def castep_magres(computer, calc_doc, seed, elec_energy_tol=1e-11):
    """Runs a NMR properties calculation on top of a completed
    SCF calculation.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info("Performing CASTEP Magres calculation...")
    magres_doc = copy.deepcopy(calc_doc)
    magres_doc["task"] = "magres"
    magres_doc["magres_task"] = calc_doc.get("magres_task", "NMR")
    magres_doc["continuation"] = "default"
    # this is just to suppress a warning that elec_energy_tol has changed
    magres_doc["elec_energy_tol"] = elec_energy_tol

    required = []
    forbidden = []

    computer.validate_calc_doc(magres_doc, required, forbidden)

    return computer.run_castep_singleshot(
        magres_doc, seed, keep=True, intermediate=True
    )
