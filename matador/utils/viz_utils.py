# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains a dirty wrapper of ase-gui for quick
visualisation, and nglview wrapper for JupyterNotebook visualisation,
and some colour definitions scraped from VESTA configs.

"""

from typing import Tuple, Optional, Dict, List, Union, Callable, Set, TYPE_CHECKING
from pathlib import Path
from functools import lru_cache

import tqdm
import numpy as np

from matador.utils.ase_utils import doc2ase
from matador.crystal import Crystal

if TYPE_CHECKING:
    import fresnel


def viz(doc):
    """Quick and dirty ase-gui visualisation from matador doc."""
    try:
        from ase.visualize import view
    except ImportError as exc:
        raise ImportError(
            "Optional dependency 'ase' is missing, please install it from PyPI with `pip install ase"
        ) from exc
    view(doc2ase(doc))
    return


@lru_cache(maxsize=1)
def get_element_colours() -> Dict[str, List[float]]:
    """Return RGB element colours from VESTA file.

    The colours file can be specified in the matadorrc.
    If unspecified, the default `../config/vesta_elements.ini` will be used.

    """
    import os
    from matador.config import load_custom_settings, SETTINGS

    default_colours_fname = (
        Path(os.path.realpath(__file__)).parent / "../config/vesta_elements.ini"
    )
    # check if element_colours has been given as an absolute path
    if SETTINGS:
        colours_fname = SETTINGS.get("plotting", {}).get("element_colours")
    else:
        colours_fname = (
            load_custom_settings().get("plotting", {}).get("element_colours")
        )

    if colours_fname is None or not os.path.isfile(colours_fname):
        colours_fname = default_colours_fname

    with open(colours_fname, "r") as f:
        flines = f.readlines()
    element_colours = dict()
    for line in flines:
        line = line.split()
        elem = line[1]
        colour = list(map(float, line[-3:]))
        element_colours[elem] = colour
    return element_colours


@lru_cache(maxsize=1)
def get_element_radii():
    """Read element radii from VESTA config file. The colours file can be
    specified in the matadorrc. If unspecified, the default
    ../config/vesta_elements.ini (relative to this file) will be used.

    """
    import os
    from matador.config import load_custom_settings, SETTINGS

    default_colours_fname = (
        Path(os.path.realpath(__file__)).parent / "../config/vesta_elements.ini"
    )
    # check if element_colours has been given as an absolute path
    if SETTINGS:
        colours_fname = SETTINGS.get("plotting", {}).get("element_colours")
    else:
        colours_fname = (
            load_custom_settings().get("plotting", {}).get("element_colours")
        )

    if colours_fname is None or not os.path.isfile(colours_fname):
        colours_fname = default_colours_fname

    with open(colours_fname, "r") as f:
        flines = f.readlines()
    element_radii = dict()
    for line in flines:
        line = line.split()
        elem = line[1]
        r = float(line[2])
        element_radii[elem] = r

    return element_radii


def nb_viz(doc, repeat=1, bonds=None):
    """Return an ipywidget for nglview visualisation in
    a Jupyter notebook or otherwise.

    Parameters:
        doc (matador.crystal.Crystal / dict): matador document to show.

    Keyword arguments:
        repeat (int): number of periodic images to include.
        bonds (str): custom bond selection.

    """
    try:
        import nglview
    except ImportError as exc:
        raise ImportError(
            "Optional dependency 'nglview' is missing, please install it from PyPI with `pip install nglview"
        ) from exc
    atoms = doc2ase(doc)
    atoms = atoms.repeat((repeat, repeat, repeat))
    view = nglview.show_ase(atoms)
    view.add_unitcell()
    view.remove_ball_and_stick()
    view.add_spacefill(radius_type="vdw", scale=0.3)
    if bonds is None:
        for elements in set(doc["atom_types"]):
            view.add_ball_and_stick(selection="#{}".format(elements.upper()))
    else:
        if not isinstance(bonds, list):
            bonds = [bonds]
        for bond in bonds:
            view.add_ball_and_stick(selection="#{}".format(bond.upper()))
    view.parameters = {"clipDist": 0}
    view.camera = "orthographic"
    view.background = "#FFFFFF"
    view.center()
    return view


def fresnel_view(
    doc: Crystal,
    standardize: Union[bool, float] = True,
    extension: Optional[Tuple[int, int, int]] = None,
    show_bonds: bool = True,
    show_cell: bool = True,
    show_compass: bool = True,
    bond_dict: Optional[Dict[Tuple[str, str], float]] = None,
    images: Union[bool, float] = True,
    pad_cell: bool = True,
    lights: Optional[Callable] = None,
    **camera_kwargs,
) -> "fresnel.Scene":
    """Return a fresnel scene visualising the input crystal.

    The scene can then be rendered/ray-traced with `fresnel.pathtrace` or
    `fresnel.preview`.

    Notes:
        Periodic images of sites will be included in the visualisation if they are
        within the desired distance tolerance of another site in the crystal, and
        are only rendered if they are bonded within the definitions of the `bond_dict`.
        This process requires the construction of two kd-trees, so disabling the `images`
        flag will yield a significant speed-up for large crystal structures.

    Parameters:
        doc: The crystal to visualize.
        standardize: Whether to standardize the crystal before visualising. If a float
            is provided, this will be used for the symmetry tolerance.
        extension: Visualize a supercell with this lattice extension.
        show_bonds: Whether to show bonds in the visualisation.
        show_cell: Whether to show the unit cell in the visualisation.
        show_compass: Whether to show a compass in the visualisation.
        bond_dict: A dictionary of bond lengths to override the default bonds, e.g.,
            `{('C', 'C'): 1.5, ('C', 'H'): 1.0}`.
        images: Whether to show periodic images in the visualisation (if given a float, this will be used
            as the clipping distance).
        pad_cell: Size to aim for all directions to be visualised, e.g. a supercell with minimum side length
            `pad_cell` will be constructed.
        lights: An optional callable that sets the lights for the scene.

    Returns:
        The fresnel scene to render.
    """
    try:
        import fresnel
    except ImportError as exc:
        raise ImportError(
            "Optional dependency 'fresnel' is missing, please install it "
            "from conda forge with `conda install -c conda-forge fresnel`"
        ) from exc

    scene = fresnel.Scene()
    if standardize:
        doc = doc.standardized()

    if show_cell:
        _draw_cell(scene, doc, compass=show_compass)

    if extension is not None:
        if not extension == (1, 1, 1):
            doc = doc.supercell(extension)
    elif pad_cell:
        doc = doc.supercell(target=pad_cell)

    clip = None
    if isinstance(images, float):
        images = True
        clip = images
    elif images and bond_dict:
        clip = max(bond_dict[d] for d in bond_dict)

    species, positions = _get_sites(doc, images=images, clip=clip)

    to_remove = set()
    if show_bonds:
        bond_sets = _draw_bonds(
            scene,
            list(zip(species, positions)),
            bond_dict=bond_dict,
            nsites=len(doc.positions_abs),
            bond_images=False,
        )

        if images:
            image_ind = len(doc.positions_abs)
            for i in range(image_ind, len(positions)):
                if i not in bond_sets:
                    to_remove.add(i)

    positions = [positions[i] for i, p in enumerate(positions) if i not in to_remove]
    species = [species[i] for i, p in enumerate(species) if i not in to_remove]

    _draw_atoms(scene, positions, species)

    a, b, c = doc.lattice_cart
    fit_orthographic_camera(scene, **camera_kwargs)
    if lights is None:
        scene.lights = fresnel.light.lightbox()
    else:
        scene.lights = lights()

    return scene


def _get_sites(
    doc, images: True, clip: float = 3.5
) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """Return the species and positions of the sites in the crystal, including
    any periodic images within the distance tolerance.

    Parameters:
        doc: The crystal to visualize.
        images: Whether to include periodic images.
        clip: The distance tolerance to other sites for periodic image inclusion.

    Returns:
        A list of species and a list of Cartesian positions.

    """
    positions = []
    species = []
    clip = clip or 3.5
    for ind, pos_cart in enumerate(doc.positions_abs):
        positions.append(pos_cart)
        species.append(doc[ind].species)
    if images:
        image_positions, image_species = get_image_atoms_near_cell_boundaries(
            doc, tolerance=clip
        )
        positions.extend(image_positions)
        species.extend(image_species)

    return species, positions


def _draw_cell(scene, doc, thickness=0.05, compass=True) -> None:
    """Add the unit cell to a fresnel scene.

    Parameters:
        scene: The fresnel scene.
        doc: The crystal to visualize.
        thickness: The line thickness to use.
        compass: Whether to draw the compass of lattice vectors.

    """
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
    unit_cell.material = fresnel.material.Material(
        color=fresnel.color.linear([0, 0, 0])
    )
    unit_cell.material.primitive_color_mix = 1
    unit_cell.material.solid = 1.0

    if compass:
        compass = fresnel.geometry.Cylinder(scene, N=3)
        compass.points[:] = [
            [[0, 0, 0], doc.lattice_cart[ind] / doc.lattice_abc[0][ind]]
            for ind in range(3)
        ]
        compass.color[0] = ([1.0, 0, 0], [1, 0, 0])
        compass.color[1] = ([0, 1.0, 0], [0, 1.0, 0])
        compass.color[2] = ([0, 0, 1.0], [0, 0, 1.0])
        compass.radius[:] = 0.1
        compass.material = fresnel.material.Material(
            color=fresnel.color.linear([0, 0, 0])
        )
        compass.material.primitive_color_mix = 1
        compass.material.solid = 0.0
        compass.material.roughness = 0.5
        compass.material.specular = 0.7
        compass.material.spec_trans = 0.0
        compass.material.metal = 0.0


def _draw_bonds(
    scene: "fresnel.Scene",
    sites: List[Tuple[str, List[List[float]]]],
    bond_dict: Dict[Tuple[str, str], float] = None,
    bond_thickness: float = 0.15,
    nsites: int = None,
    bond_images: bool = False,
) -> Dict[int, Set[int]]:
    """Add bonds to a fresnel scene.

    Parameters:
        scene: The fresnel scene.
        sites: A zipped list of species and positions.
        bond_dict: The dictionary defining bond lengths between species.
        bond_thickness: The thickness of the bonds in the visualisation.
        nsites: The number of "real" sites in the structure.
        bond_images: Whether to draw bonds between images.

    Returns:
        Indices of bonded atoms as a dictionary of sets.

    """
    import fresnel
    from collections import defaultdict

    bond_colours = []
    bonds = []
    bond_sets = defaultdict(set)
    bond_dict = {} if bond_dict is None else bond_dict
    for k, v in list(bond_dict.items()):
        bond_dict[tuple(sorted(k))] = v

    colours = get_element_colours()
    for i, i_atom in enumerate(sites):
        for j, j_atom in enumerate(sites):
            if j <= i:
                continue
            species = (sites[i][0], sites[j][0])
            dist = np.linalg.norm(np.array(sites[i][1]) - np.array(sites[j][1]))
            if bond_dict.get(tuple(sorted(species)), -1) < dist:
                continue
            # If these are both image atoms
            if not bond_images and nsites and i >= nsites and j >= nsites:
                continue
            bonds.append((sites[i][1], sites[j][1]))
            bond_sets[i].add(j)
            bond_sets[j].add(i)
            bond_colours.append(
                (
                    fresnel.color.linear(colours[sites[i][0]]),
                    fresnel.color.linear(colours[sites[j][0]]),
                )
            )
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

    return bond_sets


def _draw_atoms(
    scene: "fresnel.Scene",
    positions: List[Tuple[float, float, float]],
    species: List[str],
):
    """Add atoms to the scene, coloured by species.

    Parameters:
        scene: The fresnel scene to add the atoms to.
        positions: List of positions of the atoms to draw.
        species: List of species of the atoms to draw.

    Returns:
        atoms: The atoms geometry.
        sites: The list of sites (as species, position tuples) that were added to the scene.

    """
    import fresnel

    colours = get_element_colours()
    radii = get_element_radii()

    atoms = fresnel.geometry.Sphere(
        scene, position=positions, radius=1.0, outline_width=0.1
    )
    atoms.material = fresnel.material.Material()
    atoms.material.solid = 0.0
    atoms.material.primitive_color_mix = 1
    atoms.material.roughness = 0.8
    atoms.material.specular = 0.9
    atoms.material.spec_trans = 0.0
    atoms.material.metal = 0
    atoms.color[:] = fresnel.color.linear([colours[s] for s in species])
    atoms.radius[:] = [0.5 * radii[s] for s in species]


def rerender_scenes_to_axes(scenes, axes, renderer=None):
    """(Re)render the scenes to the axes.

    Parameters:
        scenes: List of fresnel scenes to render.
        axes: List or grid of matplolib axes to render to.
        renderer: The renderer to use (default pathtrace).

    """
    import fresnel

    _axes = axes.flatten() if not isinstance(axes, list) else axes
    if renderer is None:
        from functools import partial

        renderer = partial(fresnel.pathtrace, light_samples=32)

    for i, ax in tqdm.tqdm(enumerate(_axes), desc="Rendering scenes"):
        old_title = ax.title.get_text()
        ax.clear()
        ax.axis("off")
        ax.set_aspect("equal")
        ax.set_title(old_title)
        if i < len(scenes):
            render = renderer(scenes[i])
            ax.imshow(render[:])


def fit_orthographic_camera(
    scene: "fresnel.Scene", scale=1.1, margin=0.05, view="isometric", vectors=None
) -> None:
    """Alternative orthographic camera fitting, based on the equivalent fresnel function.

    Parameters:
        scene: The fresnel scene to visualize.
        scale: The scale of the camera.
        margin: The margin around the viewbox.
        view: Either 'front', 'isometric', 'a', 'b', 'c', 'vector' or an actual vector.
        vectors: A custom dictionary of vectors with keys 'v', 'up' and 'right' that defines
            the camera basis.

    """
    import fresnel

    if not scene.camera:
        scene.camera = fresnel.camera.Orthographic()

    if view == "front":
        view = "c"

    _default_vectors = {
        "a": dict(
            v=np.array([1, 0, 0]), up=np.array([0, 0, 1]), right=np.array([0, 1, 0])
        ),
        "b": dict(
            v=np.array([0, 1, 0]), up=np.array([1, 0, 0]), right=np.array([0, 0, 1])
        ),
        "c": dict(
            v=np.array([0, 0, 1]), up=np.array([0, 1, 0]), right=np.array([1, 0, 0])
        ),
        "isometric": dict(
            v=np.array([1, 1, 1]) / np.sqrt(3),
            up=np.array([-1, 2, -1]) / np.sqrt(6),
            right=np.array([1, 0, -1]) / np.sqrt(2),
        ),
    }

    if view == "vector" or vectors is not None:
        view = "vector"
        _default_vectors["vector"] = vectors
        if "angle" in vectors:
            angle = np.deg2rad(vectors["angle"])
            vectors["v"] = np.array([np.cos(angle), np.sin(angle), 0])
            vectors["right"] = np.array([-np.sin(angle), np.cos(angle), 0])
        if "right" not in vectors:
            vectors["right"] = np.cross(vectors["up"], vectors["v"])
        for k in ("up", "v", "right"):
            vectors[k] = np.array(vectors[k]) / np.linalg.norm(vectors[k])

    if isinstance(view, list):
        _view = np.array(view)
        _view = _view / np.linalg.norm(_view)
        up = [0, 0, 1]
        right = np.cross(_view, up)
        up = np.cross(_view, right)
        vector_view = {"v": _view, "up": up, "right": right}
        _default_vectors["vector"] = vector_view
        view = "vector"

    # raise error if the scene is empty
    if len(scene.geometry) == 0:
        raise ValueError("The scene is empty")

    # find the center of the scene
    extents = scene.get_extents()

    # choose an appropriate view automatically
    if view == "auto":
        xw = extents[1, 0] - extents[0, 0]
        yw = extents[1, 1] - extents[0, 1]
        zw = extents[1, 2] - extents[0, 2]

        if zw < 0.51 * max(xw, yw):
            view = "c"
        else:
            view = "isometric"

    v = _default_vectors[view]["v"]
    up = _default_vectors[view]["up"]

    # make a list of points of the cube surrounding the scene
    points = np.array(
        [
            [extents[0, 0], extents[0, 1], extents[0, 2]],
            [extents[0, 0], extents[0, 1], extents[1, 2]],
            [extents[0, 0], extents[1, 1], extents[0, 2]],
            [extents[0, 0], extents[1, 1], extents[1, 2]],
            [extents[1, 0], extents[0, 1], extents[0, 2]],
            [extents[1, 0], extents[0, 1], extents[1, 2]],
            [extents[1, 0], extents[1, 1], extents[0, 2]],
            [extents[1, 0], extents[1, 1], extents[1, 2]],
        ]
    )

    # find the center of the box
    center = (extents[0, :] + extents[1, :]) / 2
    points = points - center

    # determine the extent of the scene box in the up direction
    up_projection = np.dot(points, up)
    height = (1 + margin) * np.max(np.abs(up_projection)) * 2

    # determine the extent of the scene box in the view direction
    view_projection = np.dot(points, v)
    view_distance = np.max(view_projection) * scale

    # build the camera
    scene.camera.position = center + view_distance * v
    scene.camera.look_at = center
    scene.camera.up = up
    scene.camera.height = height


def get_image_atoms_near_cell_boundaries(
    doc: Crystal, tolerance: float, bond_dict: Dict[Tuple[str, str], float] = None
) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """Returns a list of atom positions and species that are within the
    distance `tolerance` of the non-image atoms.

    Parameters:
        doc: The crystal object to consider.
        tolerance: The distance tolerance in Å.
        bond_dict: If provided, use asspecies-specific tolerances.

    Returns:
        A list of species and a list of positions in Cartesian coordinates.

    """
    from itertools import product
    from scipy.spatial import KDTree

    n_iter = 1

    tree = KDTree(doc.positions_abs, compact_nodes=False, copy_data=True)

    test_displacements = [
        np.sum(np.array(doc.lattice_cart) * np.array(p), axis=-1)
        for p in product(range(-1, 2), repeat=3)
        if p != (0, 0, 0)
    ]
    positions = []
    species = []
    added = set()
    for j, pos in enumerate(doc.positions_abs):
        for d, disp in enumerate(test_displacements):
            image_pos = pos + disp
            if tree.query_ball_point(image_pos, tolerance * 1.2):
                positions.append(image_pos.tolist())
                species.append(doc.atom_types[j])
                added.add((j, d))

    # Can also "self-consistently" check for images close to the just-added images
    # Not sure this is worth enabling
    if n_iter == 2:
        tree = KDTree(positions, compact_nodes=False, copy_data=True)
        for j, pos in enumerate(doc.positions_abs):
            for d, disp in enumerate(test_displacements):
                if (j, d) not in added:
                    image_pos = pos + disp
                    if tree.query_ball_point(image_pos, tolerance):
                        positions.append(image_pos.tolist())
                        species.append(doc.atom_types[j])

    return positions, species


def fresnel_plot(
    structures: Union[Crystal, List[Crystal]],
    figsize: Optional[Tuple[float, float]] = None,
    fig_cols: Optional[int] = None,
    fig_rows: Optional[int] = None,
    labels: Union[bool, List[str]] = True,
    renderer: Optional[Callable] = None,
    camera_patches: Optional[List[Optional[Dict]]] = None,
    **fresnel_view_kwargs,
):
    """Visualize a series of structures as a grid of matplotlib plots.

    Parameters:
        structures: A list or single structure to visualize.
        figsize: The size of the figure (defaults to a 3x4 box for each structure).
        fig_cols: The number of columns in the figure.
        fig_rows: The number of rows in the figure.
        labels: If true, label structures by formula, otherwise label by the passed
            list of strings.
        renderer: A callable to use instead of the default `fresnel.pathtrace`.
            A useful alternative could be `fresnel.preview`. Can also be used to
            parameterise the fresnel rendering with e.g.,
            `partial(fresnel.pathtrace, samples=64)`.
        fresnel_view_kwargs: Any additional kwargs will be passed down to
            `matador.utils.viz_utils.fresnel_view` to control e.g., camera angles.
            The contents of any keywords containing lists with the same length as
            the passed structures will be applied to each structure individually.

    Returns:
        The matplotlib figures, axis grid and a list of fresnel scenes.

    """
    import matplotlib.pyplot as plt

    if isinstance(structures, Crystal):
        structures = [structures]

    if fig_cols is None and fig_rows is None:
        fig_cols = len(structures)
        fig_rows = 1
    elif fig_rows is None:
        fig_rows = -(len(structures) // -fig_cols)
    else:
        fig_cols = -(len(structures) // -fig_rows)

    if figsize is None:
        figsize = (2 * fig_cols, 3 * fig_rows)

    fig, axes = plt.subplots(
        nrows=fig_rows, ncols=fig_cols, figsize=figsize, squeeze=False
    )

    scenes = []

    list_kwargs = []
    for kw in fresnel_view_kwargs:
        if isinstance(fresnel_view_kwargs[kw], list) or isinstance(
            fresnel_view_kwargs[kw], tuple
        ):
            if len(fresnel_view_kwargs[kw]) == len(structures):
                list_kwargs.append(kw)

    for ind, s in enumerate(structures):
        kwargs = {k: fresnel_view_kwargs[k][ind] for k in list_kwargs}
        kwargs.update(
            {
                k: fresnel_view_kwargs[k]
                for k in fresnel_view_kwargs
                if k not in list_kwargs
            }
        )
        scenes.append(fresnel_view(s, **kwargs))

    if camera_patches is not None:
        for ind, scene in enumerate(scenes):
            patch = camera_patches[ind] or {}
            if "view" in patch:
                fit_orthographic_camera(scene, view=patch["view"])
            if "height" in patch:
                scene.camera.height = scene.camera.height + patch["height"]

    if labels is not None:
        if isinstance(labels, bool) and labels:
            labels = [s.formula_tex for s in structures]
        if labels:
            for ind, ax in enumerate(axes.flatten()):
                try:
                    ax.set_title(labels[ind])
                except IndexError:
                    pass

    rerender_scenes_to_axes(scenes, axes, renderer=renderer)

    fig.set_tight_layout(True)

    return fig, axes, scenes


def formula_to_colour(formula: str) -> List[float]:
    """Return an RGBA colour for the given chemical formula, provided
    the formula has 3 or fewer species.

    """
    from matador.utils.chem_utils import get_stoich_from_formula, get_concentration

    stoichiometry = get_stoich_from_formula(formula)
    elements = [d[0] for d in stoichiometry]
    if len(elements) > 3:
        raise RuntimeError(
            f"Cannot mix a colour for more than 3 elements, received: {stoichiometry}"
        )

    concentration = get_concentration(
        stoichiometry, elements=elements, include_end=True
    )

    return colour_from_ternary_concentration(concentration, species=elements)


def colour_from_ternary_concentration(
    conc: Union[Tuple[float, float], Tuple[float, float, float]],
    species: List[str],
    alpha: float = 1,
) -> List[float]:
    """Returns a colour for a ternary composition of the given species, mixing
    the configured VESTA colours.

    Returns:
        RGBA array.

    """
    if len(conc) == 1:
        return get_element_colours()[species[0]] + [alpha]
    elif len(conc) == 2:
        x, y = conc
        z = 1 - x - y
    else:
        x, y, z = conc
    colours = [np.asarray(get_element_colours()[s]) for s in species + ["H"]]
    return ((x * colours[0] + y * colours[1] + z * colours[2])).tolist() + [alpha]
