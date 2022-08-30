# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements some common workflow steps for use in
more complicated workflows.

"""

import copy
import logging

__all__ = ("castep_prerelax", "castep_scf")

LOG = logging.getLogger("run3")


def castep_scf(
    computer,
    calc_doc,
    seed,
    elec_energy_tol=None,
    write_checkpoint="ALL",
    required_keys=None,
    forbidden_keys=None,
):
    """Run a singleshot SCF calculation.

    Parameters:
        computer (:obj:`matador.compute.ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    Keyword arguments:
        elec_energy_tol (float or str): keyword to pass to ``elec_energy_tol``.
        write_checkpoint (bool or str): keyword to pass to CASTEP's ``write_checkpoint`` parameter. If
            ``True`` (``False``), CASTEP parameter set to ``ALL`` (``NONE``).
        required_keys (:obj:`list` of :obj:`str`): list of keys required in calc doc to perform
            the calculation.
        forbidden_keys (:obj:`list` of :obj:`str`): list of keys to scrub from calc doc to perform
            the calculation.

    Returns:
        bool: whether or not the SCF was successful.

    """
    LOG.info("Performing singleshot CASTEP SCF...")
    scf_doc = copy.deepcopy(calc_doc)

    scf_doc["write_checkpoint"] = _parse_write_checkpoint(write_checkpoint)
    scf_doc["task"] = "singlepoint"
    if elec_energy_tol is not None:
        scf_doc["elec_energy_tol"] = elec_energy_tol

    required = []
    forbidden = [
        "spectral_task",
        "spectral_kpoints_list",
        "spectral_kpoints_path",
        "spectral_kpoints_mp_spacing",
        "spectral_kpoints_path_spacing",
    ]

    computer.validate_calc_doc(scf_doc, required, forbidden)

    return computer.run_castep_singleshot(scf_doc, seed, keep=True, intermediate=True)


def castep_prerelax(
    computer,
    calc_doc,
    seed,
    write_checkpoint="all",
    required_keys=None,
    forbidden_keys=None,
):
    """Run a self-consistent (i.e. restarted) geometry optimisation.
    Optionally write a check file containing the final structure and density.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    Keyword arguments:
        write_checkpoint (bool or str): keyword to pass to CASTEP's ``write_checkpoint`` parameter. If
            ``True`` (``False``), CASTEP parameter set to ``ALL`` (``NONE``).
        required_keys (:obj:`list` of :obj:`str`): list of keys required in calc doc to perform
            the calculation.
        forbidden_keys (:obj:`list` of :obj:`str`): list of keys to scrub from calc doc to perform
            the calculation.

    Returns:
        bool: whether or not the relaxation was successful.

    """
    LOG.info("Performing CASTEP pre-relax...")

    relax_doc = copy.deepcopy(calc_doc)

    relax_doc["write_checkpoint"] = _parse_write_checkpoint(write_checkpoint)
    if "geom_max_iter" not in relax_doc:
        relax_doc["geom_max_iter"] = 100
    relax_doc["task"] = "geometryoptimisation"

    computer.validate_calc_doc(relax_doc, required_keys, forbidden_keys)
    computer.calc_doc = relax_doc

    return computer.run_castep_relaxation(intermediate=True)


def _parse_write_checkpoint(write_checkpoint):
    """Returns the appropriate value of ``write_checkpoint``."""
    if isinstance(write_checkpoint, bool):
        if not write_checkpoint:
            write_checkpoint = "NONE"
        else:
            write_checkpoint = "ALL"
    if write_checkpoint.upper() not in ("NONE", "MINIMAL", "ALL", "BOTH", "FULL"):
        LOG.warning(
            f"Invalid value of `write_checkpoint` provided: {write_checkpoint}, using 'ALL'"
        )
        write_checkpoint = "ALL"

    return write_checkpoint
