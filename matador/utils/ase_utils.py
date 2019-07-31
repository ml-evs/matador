# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements some light wrappers to
the Atomic Simulation Environment (ASE).

"""
import copy


def ase2dict(atoms):
    """ Return a simple matador-style dictionary from
    an ase.Atoms object.

    Parameters:
        atoms (ase.Atoms): input structure.

    Returns:
        dict: matador output.

    """
    from matador.utils.cell_utils import cart2abc
    doc = {}
    doc['atom_types'] = atoms.get_chemical_symbols()
    doc['lattice_cart'] = atoms.get_cell()
    doc['lattice_abc'] = cart2abc(doc['lattice_cart'])
    doc['positions_frac'] = atoms.get_scaled_positions()
    doc['num_atoms'] = len(doc['atom_types'])
    return doc


def doc2ase(doc):
    """ Convert matador document to simple ASE object. """
    from ase import Atoms

    atoms = Atoms(symbols=doc['atom_types'],
                  scaled_positions=doc['positions_frac'],
                  cell=doc['lattice_cart'],
                  pbc=True)

    atoms.info['matador'] = copy.deepcopy(doc)
    if '_id' in atoms.info['matador']:
        atoms.info['matador']['_id'] = str(atoms.info['matador']['_id'])

    return atoms
