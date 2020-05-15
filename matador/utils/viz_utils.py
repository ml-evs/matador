# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains a dirty wrapper of ase-gui for quick
visualisation, and nglview wrapper for JupyterNotebook visualisation,
and some colour definitions scraped from VESTA configs.

"""

from matador.utils.ase_utils import doc2ase


def viz(doc):
    """ Quick and dirty ase-gui visualisation from matador doc. """
    from ase.visualize import view
    view(doc2ase(doc))
    return


def get_element_colours():
    """ Read element colours from VESTA file. The colours file can be
    specified in the matadorrc. If unspecified, the default
    ../config/vesta_elements.ini will be used.

    """
    import os
    from matador.config import load_custom_settings, SETTINGS
    # check if element_colours has been given as an absolute path
    if SETTINGS:
        colours_fname = SETTINGS.get('plotting', {}).get('element_colours')
    else:
        colours_fname = load_custom_settings().get('plotting', {}).get('element_colours')

    if colours_fname is None:
        colours_fname = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/../config/vesta_elements.ini'
    elif not os.path.isfile(colours_fname):
        print('Could not find {}, please specify an absolute path. Falling back to default...'.format(colours_fname))
        # otherwise fallback to ../config/vesta_elements.ini
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

    Parameters:
        doc (matador.crystal.Crystal / dict): matador document to show.

    Keyword arguments:
        repeat (int): number of periodic images to include.
        bonds (str): custom bond selection.

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
