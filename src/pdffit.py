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
from scipy.optimize.minpack import leastsq
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
            # cast db document as diffpy Structure
            structure = doc2diffpy(doc)
            # make recipe and fit
            fit = self.make_recipe(structure, doc)
            # assess quality of fit
            try:
                result = self.regression(fit)
            except:
                print_exc()


    def make_recipe(self, structure, doc):
        """ Construct PDF with diffpy. """
        # construct a PDFContribution object
        pdf = PDFContribution("structure")

        # read experimental data
        try:
            pdf.loadData(self.input_file)
        except:
            print_failure('Failed to parse ' + data_file + '. Exiting...')
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
        print (', '.join([p.name for p in spacegroup_params]))
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


    def regression(self, fit):
        """ Assess quality of fit. """
        fit.fithooks[0].verbose = 0
        
        # free parameters one-by-one and fit
        # fit.free("scale")
        # We can now execute the fit using scipy's least square optimizer.
        print("Refine PDF using scipy's least-squares optimizer:")
        print("  variables:", fit.names)
        print("  initial values:", fit.values)
        leastsq(fit.residual, fit.values)
        print("  final values:", fit.values)
        print

        # leastsq(fit.residual, fit.values)
        # fit.free("delta2")
        # leastsq(fit.residual, fit.values)
        # fit.free("all")
        leastsq(fit.residual, fit.values)

        ContributionResult = FitResults(fit)
        ContributionResult.saveResults('test.dat')

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
        atoms.append(Atom(atype=atom,
                          xyz=asarray(doc['positions_frac'][ind])))
    title = None
    for sources in doc['source']:
        if sources.endswith('.res') or sources.endswith('.castep'):
            title = sources.split('/')[-1].split('.')[0].encode('utf-8')

    return Structure(atoms, lattice)#, title=title)
