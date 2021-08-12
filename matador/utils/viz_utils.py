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


def get_element_radii():
    """ Read element radii from VESTA config file. The colours file can be
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
    element_radii = dict()
    for line in flines:
        line = line.split()
        elem = line[1]
        r = float(line[2])
        element_radii[elem] = r

    return element_radii


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


def fresnel_view(doc, standardize=True, extension=None, show_bonds=True, show_cell=True, bond_dict=None, images=True):
    try:
        import fresnel
    except ImportError as exc:
        raise ImportError(
            "Optional dependency 'fresnel' is missing, please install it from conda forge."
        ) from exc

    import numpy as np

    scene = fresnel.Scene()
    if standardize:
        doc = doc.standardized()

    cell = None
    if show_cell:
        cell = _draw_cell(scene, doc)

    if extension is not None:
        doc = doc.supercell(extension)

    atoms = _draw_atoms(scene, doc, images=images)

    bonds = None
    if show_bonds:
        bonds = _draw_bonds(scene, doc, bond_dict=bond_dict)

    geometry = (cell, atoms, bonds)

    a, b, c = doc.lattice_cart
    middle = np.sum(doc.lattice_cart, axis=-1) / 2
    scene.camera = fresnel.iamera.Orthographic(
        position=1.3 * np.asarray(a) - 1.2 * np.asarray(b) + np.asarray(c),
        look_at=middle,
        up=c,
        height=np.sum(doc.lattice_abc[0])
    )
    scene.lights = fresnel.light.lightbox()
    return scene, geometry


def _draw_cell(scene, doc, thickness=0.05, compass=True):
    import numpy as np
    import fresnel
    vertices = np.zeros((8, 3), dtype=np.float64)
    vertices[1:4, :] = doc.lattice_cart
    vertices[4, :] = vertices[1, :] + vertices[2, :]
    vertices[5, :] = vertices[1, :] + vertices[3, :]
    vertices[6, :] = vertices[2, :] + vertices[3, :]
    vertices[7, :] = np.sum(vertices[1:4, :], axis=0)
    unit_cell = fresnel.geometry.Cylinder(scene, N=12)
    unit_cell.points[:] = [
        [vertices[0, :], vertices[1, :]],
        [vertices[0, :], vertices[2, :]],
        [vertices[0, :], vertices[3, :]],
        [vertices[1, :], vertices[4, :]],
        [vertices[2, :], vertices[4, :]],
        [vertices[1, :], vertices[5, :]],
        [vertices[3, :], vertices[5, :]],
        [vertices[2, :], vertices[6, :]],
        [vertices[3, :], vertices[6, :]],
        [vertices[4, :], vertices[7, :]],
        [vertices[5, :], vertices[7, :]],
        [vertices[6, :], vertices[7, :]],
    ]
    unit_cell.radius[:] = thickness * np.ones(12)
    unit_cell.material = fresnel.material.Material(color=fresnel.color.linear([0, 0, 0]))
    unit_cell.material.primitive_color_mix = 1
    unit_cell.material.solid = 1.

    if compass:
        compass = fresnel.geometry.Cylinder(scene, N=3)
        compass.points[:] = [[[0, 0, 0], doc.lattice_cart[ind] / doc.lattice_abc[0][ind]] for ind in range(3)]
        compass.color[0] = ([1.0, 0, 0], [1, 0, 0])
        compass.color[1] = ([0, 1.0, 0], [0, 1.0, 0])
        compass.color[2] = ([0, 0, 1.0], [0, 0, 1.0])
        compass.radius[:] = 0.1
        compass.material = fresnel.material.Material(color=fresnel.color.linear([0, 0, 0]))
        compass.material.primitive_color_mix = 1
        compass.material.solid = 0.
        compass.material.roughness = 0.5
        compass.material.specular = 0.7
        compass.material.spec_trans = 0.0
        compass.material.metal = 0.0

    return unit_cell


def _draw_bonds(scene, doc, bond_dict=None, bond_thickness=0.2):
    import numpy as np
    import fresnel
    bonds = []
    bond_colours = []
    rendered_bonds = set()
    colours = get_element_colours()
    for bond in doc.network.edges:
        i, j, k = bond
        indices = tuple(sorted([i, j]))
        if indices in rendered_bonds:
            continue
        rendered_bonds.add(indices)
        if bond_dict is not None:
            species = (doc[i].species, doc[j].species)
            dist = doc.network[i][j][k]["dist"]
            if bond_dict.get(tuple(sorted(species)), 0.0) < dist:
                continue
        bonds.append((doc[i].coords_cartesian.tolist(), doc[j].coords_cartesian.tolist()))
        bond_colours.append((
            fresnel.color.linear(colours[doc[i].species]).tolist(),
            fresnel.color.linear(colours[doc[j].species]).tolist()
        ))
    if bonds:
        bonds_geom = fresnel.geometry.Cylinder(
            scene,
            N=len(bonds),
            color=bond_colours,
        )
        bonds_geom.points[:] = bonds
        bonds_geom.radius[:] = bond_thickness * np.ones(len(bonds))
        bonds_geom.material.primitive_color_mix = 1
        bonds_geom.outline_width = 0.05

    return bonds_geom


def _draw_atoms(scene, doc, images=True):
    import fresnel
    import numpy as np

    colours = get_element_colours()
    radii = get_element_radii()
    positions = []
    species = []
    for ind, (pos_cart, pos_frac) in enumerate(zip(doc.positions_abs, doc.positions_frac)):
        positions.append(pos_cart)
        EPS = 1e-3
        species.append(doc[ind].species)
        if images:
            marked = []
            for j, e in enumerate(pos_frac):
                if e < EPS:
                    marked.append(j)
            for j in marked:
                pos_image = np.array(pos_cart) + doc.lattice_cart[j]
                pos_frac_image = np.array(pos_frac)
                pos_frac_image[j] += 1
                positions.append(pos_image)
                species.append(doc[ind].species)

    atoms = fresnel.geometry.Sphere(scene, position=positions, radius=1.0, outline_width=0.1)
    atoms.material = fresnel.material.Material()
    atoms.material.solid = 0.
    atoms.material.primitive_color_mix = 1
    atoms.material.roughness = 0.8
    atoms.material.specular = 0.9
    atoms.material.spec_trans = 0.0
    atoms.material.metal = 0
    atoms.color[:] = fresnel.color.linear([colours[s] for s in species])
    atoms.radius[:] = [0.5 * radii[s] for s in species]

    return atoms
