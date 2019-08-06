# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements some light wrappers to
the Atomic Simulation Environment (ASE).

"""
import copy
import matador.crystal


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


def doc2ase(doc, add_keys_to_info=True):
    """ Convert matador document to simple ASE object.

    Parameters:
        doc (dict/:obj:`matador.crystal.Crystal`): matador document  or
            `Crystal` containing the structure.

    Keyword arguments:
        add_keys_to_info (bool): whether or not to add the keys from the
            matador document to the info section of the Atoms object.

    """
    from ase import Atoms

    atoms = Atoms(symbols=doc['atom_types'],
                  scaled_positions=doc['positions_frac'],
                  cell=doc['lattice_cart'],
                  pbc=True)

    if add_keys_to_info:
        if isinstance(doc, matador.crystal.Crystal):
            atoms.info['matador'] = doc._data
        else:
            atoms.info['matador'] = copy.deepcopy(doc)
        if '_id' in atoms.info['matador']:
            atoms.info['matador']['_id'] = str(atoms.info['matador']['_id'])

    return atoms
