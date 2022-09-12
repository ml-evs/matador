# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements some light wrappers to
the Atomic Simulation Environment (ASE).

"""

import copy
from typing import Union
from matador.utils.chem_utils import get_stoich
from matador.utils.cell_utils import get_spacegroup_spg
from matador.crystal import Crystal

__all__ = ["ase2dict", "doc2ase"]


def ase2dict(atoms, as_model=False) -> Union[dict, Crystal]:
    """Return a matador document (dictionary or :obj:`Crystal`)
    from an `ase.Atoms` object.

    Parameters:
        atoms (ase.Atoms): input structure.

    Keyword arguments:
        as_model (bool): if `True`, return a
            Crystal instead of a dictionary.

    Returns:
        Union[dict, Crystal]: matador output.

    """
    from matador.utils.cell_utils import cart2abc

    doc = {}

    # sort atoms, then their positions
    doc["atom_types"] = atoms.get_chemical_symbols()
    inds = [i[0] for i in sorted(enumerate(doc["atom_types"]), key=lambda x: x[1])]
    doc["positions_frac"] = atoms.get_scaled_positions().tolist()
    doc["positions_frac"] = [doc["positions_frac"][ind] for ind in inds]
    doc["atom_types"] = [doc["atom_types"][ind] for ind in inds]
    try:
        doc["lattice_cart"] = atoms.get_cell().tolist()
    except AttributeError:
        doc["lattice_cart"] = atoms.get_cell().array.tolist()
    doc["lattice_abc"] = cart2abc(doc["lattice_cart"])
    doc["num_atoms"] = len(doc["atom_types"])
    doc["stoichiometry"] = get_stoich(doc["atom_types"])
    doc["cell_volume"] = atoms.get_volume()
    doc["elems"] = {atom for atom in doc["atom_types"]}
    doc["num_fu"] = doc["num_atoms"] / int(
        sum(doc["stoichiometry"][i][1] for i in range(len(doc["stoichiometry"])))
    )
    doc["space_group"] = get_spacegroup_spg(doc, symprec=0.001)

    if hasattr(atoms, "info"):
        doc["ase_info"] = copy.deepcopy(atoms.info)

    if as_model:
        doc = Crystal(doc)

    return doc


def doc2ase(doc: Union[dict, Crystal], add_keys_to_info=True):
    """Convert matador document to simple ASE object.

    Parameters:
        doc (dict/:obj:`Crystal`): matador document or
            `Crystal` containing the structure.

    Keyword arguments:
        add_keys_to_info (bool): whether or not to add the keys from the
            matador document to the info section of the Atoms object.

    """
    from ase import Atoms

    atoms = Atoms(
        symbols=doc["atom_types"],
        scaled_positions=doc["positions_frac"],
        cell=doc["lattice_cart"],
        pbc=True,
    )

    if add_keys_to_info:
        if isinstance(doc, Crystal):
            atoms.info["matador"] = doc._data
        else:
            atoms.info["matador"] = copy.deepcopy(doc)
        if "_id" in atoms.info["matador"]:
            atoms.info["matador"]["_id"] = str(atoms.info["matador"]["_id"])

    return atoms
