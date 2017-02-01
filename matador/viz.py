# encoding: utf-8
""" This file is a dirty wrapper of ase-gui for
quick visualisation of structures.
"""


def viz(doc):
    from ase.visualize import view
    view(doc2ase(doc))
    return


def doc2ase(doc):
    from ase import Atoms
    return Atoms(symbols=doc['atom_types'],
                 scaled_positions=doc['positions_frac'],
                 cell=doc['lattice_cart'],
                 pbc=True)
