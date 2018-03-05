""" Simple ASE wrapper to read in cif files as
matador documents.
"""


def cif2dict(fname: str):
    """ Read cif file into ASE object,
    then convert ASE Atoms into matador document.

    Input:

        | fname: str, cif filename
    Returns:

        | doc: dict, simple matador document

    """
    import ase.io
    from matador.utils.cell_utils import cart2abc
    fname = fname.replace('.cif', '')

    try:
        atoms = ase.io.read(fname + '.cif')
        doc = dict()
        doc['atom_types'] = atoms.get_chemical_symbols()
        doc['lattice_cart'] = atoms.get_cell()
        doc['lattice_abc'] = cart2abc(doc['lattice_cart'])
        doc['positions_frac'] = atoms.get_scaled_positions()
        doc['num_atoms'] = len(doc['atom_types'])
        return doc, True

    except Exception as oops:
        print(oops)
        return oops, False
