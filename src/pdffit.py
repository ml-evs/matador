#!/usr/bin/python
# coding: utf-8
""" This file implements the PDFFitter class,
which attempts fitting of structures to an
experimental diffraction input to structures
in the database with the diffpy package.
"""
from __future__ import print_function

# matador functionality
from print_utils import print_failure, print_warning

# standard library
from traceback import print_exc
# external libraries
from scipy.optimize import leastsq
from numpy import zeros_like
import matplotlib.pyplot as plt
# diffpy
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe, FitResults
from diffpy.srfit.structure import constrainAsSpaceGroup


class PDFFitter:
    """ Take in a completed query and find the best fit of
    one or two structures to the experimental PDF.
    """
    def __init__(self, cursor,
                 dx=0.01, xmin=1, xmax=50, two_phase=False,
                 **kwargs):
        """ Read in args and begin computing PDFs,
        call fitting functions then display best
        candidates.
        """
        self.args = kwargs
        self.input_file = self.args.get('file')
        self.cursor = list(cursor)
        self.dx = dx
        self.xmin = xmin
        self.xmax = xmax
        self.two_phase = two_phase
        self.total_num_candidates = len(self.cursor)

        for ind, doc in enumerate(self.cursor):
            print(doc['text_id'][0], doc['text_id'][1])
            # cast db document as diffpy Structure
            structure = doc2diffpy(doc)
            # make recipe and fit
            fit = self.make_recipe(structure, doc)
            # assess quality of fit
            try:
                self.regression(fit, structure.title)
            except:
                print_exc()
            # self.plot_fit(fit, structure)

    def plot_fit(self, fit, structure):
        """ Plot results. """
        fig_name = structure.title + ".png"
        # Plot the observed and refined PDF.
        # Get the experimental data from the recipe
        r = fit.Contribution.profile.x
        gobs = fit.Contribution.profile.y
        # Get the calculated PDF and compute the difference between the calculated and
        # measured PDF
        gcalc = fit.Contribution.evaluate()
        baseline = 1.1 * gobs.min()
        gdiff = gobs - gcalc
        # Plot!
        plt.figure(figsize=(5, 5))
        plt.plot(r, gobs, 'bo', label="G(r) data",
                 markerfacecolor='none', markeredgecolor='b')
        plt.plot(r, gcalc, 'r-', label="G(r) fit")
        plt.plot(r, gdiff + baseline, 'g-', label="G(r) diff")
        plt.plot(r, zeros_like(r) + baseline, 'k:')
        plt.xlabel(r"r ($\AA$)")
        plt.ylabel(r"G ($\AA^{-2}$)")
        plt.legend()
        plt.savefig(fig_name, bbox_inches='tight')
        plt.show()

    def make_recipe(self, structure, doc):
        """ Construct PDF with diffpy. """
        # construct a PDFContribution object
        pdf = PDFContribution("structure")

        # read experimental data
        try:
            pdf.loadData(self.input_file)
        except:
            print_failure('Failed to parse ' + self.input_file + '. Exiting...')
            exit()

        pdf.setCalculationRange(self.xmin, self.xmax, self.dx)
        pdf.addStructure("Contribution", structure)

        # create FitRecipe to calculate PDF with chosen fit variable
        fit = FitRecipe()
        fit.addContribution(pdf)
        # configure variables and add to recipe
        if doc['space_group'] != 'xxx':
            spacegroup_params = constrainAsSpaceGroup(pdf.Contribution.phase,
                                                      str(doc['space_group']))
        else:
            print_warning('Invalid space group... skipping')
            raise RuntimeError('Invalid space group for', doc['text_id'][0], doc['text_id'][1])
        print('Space group parameters:')
        print(', '.join([param.name for param in spacegroup_params]))
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
        """ Assess quality of fit. """

        fit.fithooks[0].verbose = 0

        # We can now execute the fit using scipy's least square optimizer.
        print("Refine PDF using scipy's least-squares optimizer:")
        print("  variables:", fit.names)
        print("  initial values:", fit.values)
        # free parameters one-by-one and fit
        fit.free("scale")
        leastsq(fit.residual, fit.values)
        print("  freed scale:", fit.values)
        fit.free("delta2")
        leastsq(fit.residual, fit.values)
        print("  freed delta2:", fit.values)
        fit.free("all")
        leastsq(fit.residual, fit.values)
        print("  freed all:", fit.values)

        ContributionResult = FitResults(fit)
        ContributionResult.saveResults(title+'.results')

        return


def doc2diffpy(doc):
    """ Convert doc into diffpy Structure object. """
    from numpy import asarray
    from diffpy.Structure.atom import Atom
    from diffpy.Structure.lattice import Lattice
    from diffpy.Structure.structure import Structure
    
    lattice = Lattice(a=doc['lattice_abc'][0][0],
                      b=doc['lattice_abc'][0][1],
                      c=doc['lattice_abc'][0][2],
                      alpha=doc['lattice_abc'][1][0],
                      beta=doc['lattice_abc'][1][1],
                      gamma=doc['lattice_abc'][1][2]
                      )
    atoms = []
    for ind, atom in enumerate(doc['atom_types']):
        # encode atype as utf-8 or you will waste hours of your life
        atoms.append(Atom(atype=atom.encode('utf-8'),
                          xyz=asarray(doc['positions_frac'][ind])))
    title = None
    for sources in doc['source']:
        if sources.endswith('.res') or sources.endswith('.castep'):
            title = sources.split('/')[-1].split('.')[0].encode('utf-8')

    return Structure(atoms, lattice, title=title)
