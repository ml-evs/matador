# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module implements the CastepElasticWorkflow, which performs
elastic tensor calculations with CASTEP. So far, only bulk modulus
calculations (i.e. varying volume from equilibrium) have been implemented.

"""


import copy
import os
import logging
import numpy as np
from matador.workflows import Workflow
from matador.crystal.elastic import get_equation_of_state
from matador.workflows.castep.common import castep_prerelax
from matador.scrapers import cell2dict

LOG = logging.getLogger("run3")


def castep_elastic(computer, calc_doc, seed, **kwargs):
    """Perform a calculation of the elastic tensor on a system.
    Currently only calculation of the bulk modulus is implemented.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): dictionary of structure and calculation
            parameters.
        seed (str): root seed for the calculation.

    Raises:
        RuntimeError: if any part of the calculation fails.

    Returns:
        bool: True if Workflow completed successfully, or False otherwise.

    """
    workflow = CastepElasticWorkflow(computer, calc_doc, seed, **kwargs)
    return workflow.success


class CastepElasticWorkflow(Workflow):
    """Perform a "full" spectral calculation on a system, i.e. first
    perform an SCF then interpolate to different kpoint paths/grids to
    form DOS and dispersions. Optionally use OptaDOS for post-processing
    of DOS.

    Attributes:
        computer (:obj:`ComputeTask`): the object that calls CASTEP.
        calc_doc (dict): the interim dictionary of structural and
            calculation parameters.
        seed (str): the root seed for the calculation.
        volume_rescale (list): a list of floats containing the amounts
            by which to rescale volumes when performing bulk modulus calculation.
        success (bool): the status of the Workflow: only set to True after
            post-processing method completes.

    """

    def preprocess(self):
        """Decide which parts of the Workflow need to be performed,
        and set the appropriate CASTEP parameters.

        """
        num_volumes = self.workflow_params.get("num_volumes", 9)
        self.plot = self.workflow_params.get("plot", True)

        self.volume_rescale = np.cbrt(
            np.geomspace(0.7, 1.2, num=num_volumes, endpoint=True)
        ).tolist()
        self.volume_rescale.append(1.0)
        self._completed_volumes = []
        LOG.info(
            "Preprocessing completed: run3 bulk modulus calculation options {}".format(
                self.volume_rescale
            )
        )

        # if there is not geom_force_tol parameter set, then don't pre-relax
        if self.calc_doc.get("geom_force_tol"):
            self.add_step(
                castep_elastic_prerelax,
                "relax",
                clean_after=True,  # clean up after geometry step as seed is going to change
                output_exts=["-out.cell"],
            )

        for volume in self.volume_rescale:
            self.add_step(
                castep_rescaled_volume_scf, "rescaled_volume_scf", rescale=volume
            )

    def postprocess(self):
        """Fit some equations of state, then save a plot and a datafile."""
        results = get_equation_of_state(self.seed + ".bulk_mod", plot=self.plot)
        if "summary" in results:
            for line in results["summary"]:
                print(line)
            print("Writing summary to {seed}.bulk_mod.results".format(seed=self.seed))
            with open(self.seed + ".bulk_mod.results", "w") as f:
                f.writelines(results["summary"])


def castep_rescaled_volume_scf(computer, calc_doc, seed, rescale=1):
    """Run a singleshot SCF calculation.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    Keyword arguments:
        rescale (float): amount to rescale volume by.

    """
    assert rescale > 0
    LOG.info("Performing CASTEP SCF on volume rescaled by {:.2f}.".format(rescale**3))
    scf_doc = copy.deepcopy(calc_doc)

    # get relaxed cell from previous step
    fname = f"completed/{seed}-out.cell"
    if not os.path.isfile(fname):
        raise RuntimeError(f"Could not find relaxed geometry in {fname}.")
    relaxed_cell, s = cell2dict(fname, positions=True, lattice=True)

    scf_doc["lattice_cart"] = relaxed_cell["lattice_cart"]
    scf_doc["positions_frac"] = relaxed_cell["positions_frac"]

    LOG.debug("Found equilibrium geometry: {scf_doc['lattice_cart']}")

    for i in range(3):
        for k in range(3):
            scf_doc["lattice_cart"][i][k] *= rescale

    scf_doc["task"] = "singlepoint"
    bulk_mod_seed = seed + ".bulk_mod"
    computer.seed = bulk_mod_seed

    return computer.run_castep_singleshot(
        scf_doc, bulk_mod_seed, keep=True, intermediate=True
    )


def castep_elastic_prerelax(computer, calc_doc, seed):
    """Run a geometry optimisation before re-scaling volumes SCF-style calculation.

    Parameters:
        computer (:obj:`ComputeTask`): the object that will be calling CASTEP.
        calc_doc (dict): the structure to run on.
        seed (str): root filename of structure.

    """
    LOG.info("Performing CASTEP elastic pre-relax...")
    relax_doc = copy.deepcopy(calc_doc)
    relax_doc["write_cell_structure"] = True

    return castep_prerelax(computer, relax_doc, seed)
