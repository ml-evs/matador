# coding: utf-8
""" This file implements the PDFFitter class,
which attempts fitting of structures to an
experimental diffraction input to structures
in the database with the diffpy package.
"""

from __future__ import print_function

# matador modules
from matador.utils.print_utils import print_failure, print_notify
from matador.utils.cell_utils import abc2cart
from matador.utils.chem_utils import get_atomic_number

# external libraries
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import spglib as spg
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe, FitResults
from diffpy.srfit.structure import constrainAsSpaceGroup

# standard library
import multiprocessing as mp
from traceback import print_exc
from os.path import isfile


class PDFFitter(object):
    """Take in a completed query and find the best fit of
    one or two structures to the experimental PDF.
    """

    def __init__(self, cursor, dx=0.01, xmin=1, xmax=50, two_phase=False, **kwargs):
        """Read in args and begin computing PDFs,
        call fitting functions then display best
        candidates.
        """
        self.args = kwargs
        self.nprocesses = self.args.get("num_processes")
        if self.nprocesses is None:
            self.nprocesses = 2
        self.input_file = self.args.get("file")
        self.cursor = list(cursor)
        self.dx = dx
        self.xmin = xmin
        self.xmax = xmax
        self.two_phase = two_phase
        self.total_num_candidates = len(self.cursor)
        self.structures = []
        self.space_groups = []
        # collect structures and space groups
        for ind, doc in enumerate(self.cursor):
            # cast db document as diffpy Structure
            self.structures.append(doc2diffpy(doc))
            self.space_groups.append(doc["space_group"])

        self.jobs_fname = "pdf_jobs.txt"
        if not isfile(self.jobs_fname):
            with open(self.jobs_fname, "a"):
                pass
        self.completed_fname = "pdf_complete.txt"
        if not isfile(self.completed_fname):
            with open(self.completed_fname, "a"):
                pass
        self.failed_fname = "pdf_failed.txt"
        if not isfile(self.failed_fname):
            with open(self.failed_fname, "a"):
                pass
        self.spawn()

    def perform_fits(self):
        """Scan for fit that has not yet been
        done, then do it.
        """
        for ind, structure in enumerate(self.structures):
            running = False
            with open(self.jobs_fname, "rw") as job_file:
                flines = job_file.readlines()
                for line in flines:
                    if structure.title in line:
                        running = True
                        break
            if not running:
                sg = self.space_groups[ind]
                with open(self.jobs_fname, "a") as job_file:
                    job_file.write(structure.title + "\n")
                try:
                    self.fit(structure, sg)
                    with open(self.completed_fname, "a") as job_file:
                        job_file.write(structure.title + "\n")
                except (KeyboardInterrupt, SystemExit, RuntimeError):
                    print_exc()
                    raise SystemExit
                except:
                    with open(self.failed_fname, "a") as job_file:
                        job_file.write(structure.title + "\n")
                    # print_exc()
                    pass
        return

    def fit(self, structure, sg):
        """Prepare to fit, fit, then plot."""
        fit = self.make_recipe(structure, sg)
        try:
            self.regression(fit, structure.title)
            self.write_to_file(fit, structure.title)
            self.plot_fit(fit, structure)
        except:
            print_exc()
        return

    def plot_fit(self, fit, title, num=None):
        """Plot results."""
        if num is None:
            fig_name = title + ".pdf"
        else:
            fig_name = title + str(num) + ".pdf"
        try:
            plt.style.use("bmh")
        except:
            pass
        # Plot the observed and refined PDF.
        # Get the experimental data from the recipe
        r_expt = fit.Contribution.profile.x
        g_expt = fit.Contribution.profile.y
        # Get the calculated PDF and compute the difference between the calculated and
        # measured PDF
        g_calc = fit.Contribution.evaluate()
        baseline = 1.2 * g_expt.min()
        g_diff = g_expt - g_calc

        fig = plt.figure(figsize=(5, 3))
        ax = fig.add_subplot(111)
        ax.plot(r_expt, g_expt, "o", label="G(r) expt.", markerfacecolor="none", lw=1)
        ax.plot(r_expt, g_calc, "-", label="G(r) calc. fit")
        ax.plot(r_expt, g_diff + baseline, "-", label="G(r) diff")
        ax.axhline(baseline, ls="--", c="k", lw=1)
        ax.set_xlabel(r"r ($\mathrm{\AA}$)")
        ax.set_xticklabels(ax.get_xticks())
        ax.set_ylabel(r"G(r) ($\mathrm{\AA}^{-2}$)")
        ax.set_yticklabels(ax.get_yticks())
        ax.grid("off")
        ax.legend(
            loc="upper center",
            ncol=3,
            fontsize=12,
            shadow=True,
            bbox_to_anchor=(0.5, 1.2),
        )
        plt.savefig(fig_name, bbox_inches="tight", dpi=300)
        return

    def make_recipe(self, structure, sg):
        """Construct PDF with diffpy."""

        # construct a PDFContribution object
        pdf = PDFContribution("Contribution")
        # read experimental data
        try:
            pdf.loadData(self.input_file)
        except:
            print_failure("Failed to parse " + self.input_file + ". Exiting...")
            exit()

        print("Constructing PDF object for", structure.title)

        pdf.setCalculationRange(self.xmin, self.xmax, self.dx)
        pdf.addStructure("Contribution", structure)

        # create FitRecipe to calculate PDF with chosen fit variable
        fit = FitRecipe()
        fit.addContribution(pdf)
        # configure variables and add to recipe
        if sg != "xxx" and sg is not None:
            print(sg)
            spacegroup_params = constrainAsSpaceGroup(pdf.Contribution.phase, sg)
        else:
            cart_lat = abc2cart(
                [
                    [structure.lattice.a, structure.lattice.b, structure.lattice.c],
                    [
                        structure.lattice.alpha,
                        structure.lattice.beta,
                        structure.lattice.gamma,
                    ],
                ]
            )
            positions_frac = structure.xyz
            atomic_numbers = []
            for atom in structure.element:
                atomic_numbers.append(get_atomic_number(atom))
            cell = (cart_lat, positions_frac, atomic_numbers)
            sg = int(
                spg.get_spacegroup(cell, symprec=1e-2)
                .split(" ")[1]
                .replace("(", "")
                .replace(")", "")
            )
            spacegroup_params = constrainAsSpaceGroup(pdf.Contribution.phase, sg)
        # print('Space group parameters:')
        # print(', '.join([param.name for param in spacegroup_params]))
        # iterate through spacegroup params and activate them
        for param in spacegroup_params.latpars:
            fit.addVar(param)
        for param in spacegroup_params.xyzpars:
            fit.addVar(param, fixed=True)
        # these next parameters are taken from Martin's PDFht.py,
        # though I have a feeling that was not their origin...
        # set initial ADP parameters
        for param in spacegroup_params.adppars:
            fit.addVar(param, value=0.03, fixed=True)
        # overall scale of PDF and delta2 parameter for correlated motion - from PDFht.py
        fit.addVar(pdf.scale, 1, fixed=True)
        fit.restrain(pdf.scale, lb=0, ub=0.1, scaled=True)
        fit.addVar(pdf.Contribution.delta2, 5, fixed=True)
        fit.restrain(pdf.Contribution.delta2, lb=1, ub=10, scaled=True)
        # fix Qdamp based on information about "our beamline": yep, no idea
        fit.addVar(pdf.qdamp, 0.03, fixed=True)
        fit.restrain(pdf.qdamp, lb=0, ub=0.1, scaled=True)

        return fit

    def regression(self, fit, title):
        """Apply least squares to the free parameters."""
        print("Fitting PDF of", title, "to expt. data.")
        fit.fithooks[0].verbose = 0

        # We can now execute the fit using scipy's least square optimizer.
        # free parameters one-by-one and fit
        self.plot_fit(fit, title, num=0)
        print(title, "[1/3]")
        fit.free("scale")
        leastsq(fit.residual, fit.values)
        self.plot_fit(fit, title, num=1)
        print(title, "[2/3]")
        fit.free("delta2")
        leastsq(fit.residual, fit.values)
        self.plot_fit(fit, title, num=2)
        print(title, "[3/3]")
        fit.free("all")
        leastsq(fit.residual, fit.values)
        self.plot_fit(fit, title, num=3)
        print(title, "Done!")

        ContributionResult = FitResults(fit)
        ContributionResult.saveResults(title + ".results")

        return

    def write_to_file(self, fit, title):
        """Write final fit to file."""
        from numpy import savetxt, zeros

        r_expt = fit.Contribution.profile.x
        g_expt = fit.Contribution.profile.y
        g_calc = fit.Contribution.evaluate()
        print(g_calc)
        fit_data = zeros((3, len(r_expt)))
        fit_data[0] = r_expt
        fit_data[1] = g_expt
        fit_data[2] = g_calc
        savetxt(
            title + ".fit",
            fit_data.T,
            header="\tr (Angstrom)\t\t\tG_expt(r)\t\t\tG_calc(r)",
        )
        return

    def spawn(self):
        """Spawn processes to perform PDF fitting."""
        print_notify("Performing " + str(self.nprocesses) + " concurrent fits.")
        procs = []
        for ind in range(self.nprocesses):
            procs.append(mp.Process(target=self.perform_fits))
        try:
            for proc in procs:
                proc.start()
        except (KeyboardInterrupt, SystemExit, RuntimeError):
            for proc in procs:
                proc.terminate()
            exit("Killing running jobs and exiting...")


def doc2diffpy(doc):
    """Convert doc into diffpy Structure object."""
    from numpy import asarray
    from diffpy.Structure.atom import Atom
    from diffpy.Structure.lattice import Lattice
    from diffpy.Structure.structure import Structure

    lattice = Lattice(
        a=doc["lattice_abc"][0][0],
        b=doc["lattice_abc"][0][1],
        c=doc["lattice_abc"][0][2],
        alpha=doc["lattice_abc"][1][0],
        beta=doc["lattice_abc"][1][1],
        gamma=doc["lattice_abc"][1][2],
    )
    atoms = []
    for ind, atom in enumerate(doc["atom_types"]):
        # encode atype as utf-8 or you will waste hours of your life
        atoms.append(
            Atom(atype=atom.encode("utf-8"), xyz=asarray(doc["positions_frac"][ind]))
        )
    title = None
    for sources in doc["source"]:
        if sources.endswith(".res") or sources.endswith(".castep"):
            title = sources.split("/")[-1].split(".")[0].encode("utf-8")

    return Structure(atoms, lattice, title=title)
