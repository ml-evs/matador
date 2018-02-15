def construct_network(structure, coordination_cutoff=1.2, num_images=1, debug=False, separate_images=False):
    # iterate over all pairs of atoms in check minimum distance
    import networkx as nx
    crystal_graph = nx.MultiDiGraph()
    atoms = structure.sites
    from itertools import product
    import numpy as np
    from copy import deepcopy
    images = product(range(-num_images, num_images+1), repeat=3)
    image_trans = []
    for image in images:
        image_trans.append(np.zeros((3)))
        for k in range(3):
            image_trans[-1] += image[k] * np.asarray(structure.lattice_cart[k])

    image_number = 0

    for i in range(len(atoms)):
        crystal_graph.add_node(i, species=atoms[i].species, image=False)
        min_dist = 1e20
        for j in range(len(atoms)):
            if i == j:
                continue
            for displacement in image_trans:
                image_atom = deepcopy(atoms[j])
                for k in range(3):
                    image_atom._coords['cartesian'][k] += displacement[k]
                dist = atoms[i].distance_between_sites(image_atom)
                if dist < min_dist:
                    min_dist = dist
        if debug:
            print(min_dist)
        for j in range(len(atoms)):
            if i == j:
                continue
            for displacement in image_trans:
                image_atom = deepcopy(atoms[j])
                for k in range(3):
                    image_atom._coords['cartesian'][k] += displacement[k]
                dist = atoms[i].distance_between_sites(image_atom)
                if dist <= min_dist*coordination_cutoff:
                    if debug:
                        print(i, j, min_dist, displacement)
                    if separate_images and all([val <= 0+1e-8 for val in displacement]):
                        image_number += 1
                        crystal_graph.add_node(j+image_number, species=atoms[j].species, image=True)
                        crystal_graph.add_edge(i, j+image_number, dist=dist)
                    else:
                        crystal_graph.add_edge(i, j, dist=dist)

    return crystal_graph


def draw_network(structure, layout=None):
    import networkx as nx
    from matador.viz import ELEMENT_COLOURS
    import matplotlib.pyplot as plt
    if layout is None:
        pos = nx.spring_layout(structure.network)
    else:
        pos = layout
    try:
        network = structure.network
    except:
        network = structure
    elem_map = ELEMENT_COLOURS
    colours = [elem_map.get(data['species']) for node, data in network.nodes.data()]
    labels = {node: str(data['species']) for node, data in network.nodes.data()}
    nx.draw_networkx_nodes(network, pos, node_color=colours, edgecolors='black', linewidths=2, node_size=1000)
    nx.draw_networkx_edges(network, pos, linewidths=2, node_size=1000)
    nx.draw_networkx_labels(network, pos, labels=labels)
    plt.axis('off')
