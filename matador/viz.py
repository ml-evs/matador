# encoding: utf-8
""" This file is a dirty wrapper of ase-gui for
quick visualisation of structures.
"""
ELEMENT_COLOURS = {
    'K': '#66236D',
    'Sn': '#938CAF',
    'P': '#D66814'
}


def viz(doc):
    """ Quick and dirty ase-gui visualisation from matador doc. """
    from ase.visualize import view
    view(doc2ase(doc))
    return


def nb_viz(doc, repeat=1):
    """ Return an ipywidget for nglview visualisation in
    a Jupyter notebook or otherwise.
    """
    import nglview
    atoms = doc2ase(doc)
    atoms = atoms.repeat((repeat, repeat, repeat))
    view = nglview.show_ase(atoms)
    view.add_unitcell()
    view.add_spacefill(radius_type='vdw', scale=0.3)
    for elements in set(doc['atom_types']):
        view.add_ball_and_stick(selection='#{}'.format(elements).upper())
    view.background = '#FFFFFF'
    view.camera = 'orthographic'
    view.center()
    return view


def doc2ase(doc):
    """ Convert matador document to simple ASE object. """
    from ase import Atoms
    return Atoms(symbols=doc['atom_types'],
                 scaled_positions=doc['positions_frac'],
                 cell=doc['lattice_cart'],
                 pbc=True)
