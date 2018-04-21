# encoding: utf-8
""" This file is a dirty wrapper of ase-gui for
quick visualisation of structures.
"""


def viz(doc):
    """ Quick and dirty ase-gui visualisation from matador doc. """
    from ase.visualize import view
    view(doc2ase(doc))
    return


def get_element_colours():
    """ Read element colours from VESTA file. """
    import os
    colours_fname = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/../config/vesta_elements.ini'
    with open(colours_fname, 'r') as f:
        flines = f.readlines()
    element_colours = dict()
    for line in flines:
        line = line.split()
        elem = line[1]
        colour = list(map(float, line[-3:]))
        element_colours[elem] = colour
    return element_colours


def nb_viz(doc, repeat=1, bonds=None):
    """ Return an ipywidget for nglview visualisation in
    a Jupyter notebook or otherwise.
    """
    import nglview
    atoms = doc2ase(doc)
    atoms = atoms.repeat((repeat, repeat, repeat))
    view = nglview.show_ase(atoms)
    view.add_unitcell()
    view.remove_ball_and_stick()
    view.add_spacefill(radius_type='vdw', scale=0.3)
    if bonds is None:
        for elements in set(doc['atom_types']):
            view.add_ball_and_stick(selection='#{}'.format(elements.upper()))
    else:
        if not isinstance(bonds, list):
            bonds = [bonds]
        for bond in bonds:
            view.add_ball_and_stick(selection='#{}'.format(bond.upper()))
    view.parameters = {'clipDist': 0}
    view.camera = 'orthographic'
    view.background = '#FFFFFF'
    view.center()
    return view


def doc2ase(doc):
    """ Convert matador document to simple ASE object. """
    from ase import Atoms
    return Atoms(symbols=doc['atom_types'],
                 scaled_positions=doc['positions_frac'],
                 cell=doc['lattice_cart'],
                 pbc=True)


ELEMENT_COLOURS = get_element_colours()
