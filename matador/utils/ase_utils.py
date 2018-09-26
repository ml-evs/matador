# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements some light wrappers to
the Atomic Simulation Environment (ASE).

"""


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
