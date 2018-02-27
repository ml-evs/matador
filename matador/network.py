""" This file implements turning matador Crystal objects
into CrystalGraph objects.
"""
import networkx as nx


class CrystalGraph(nx.MultiDiGraph):
    def __init__(self, structure=None, graph=None, coordination_cutoff=1.2, num_images=1, debug=False, separate_images=False):
        """ Create networkx.MultiDiGraph object with overloaded functions.

        Args:

            | structure: matador.Crystal, crystal structure to network-ify
            | graph: nx.MultiDiGraph, initialise from graph
            | coordination_cutoff: float, max multiplier of first coordination sphere for edge drawing
            | num_images : int, number of periodic images to include in each direction
            | separate_images: bool, whether or not to include image atoms as new nodes

        """

        super().__init__()

        if graph is None and structure is None:
            raise RuntimeError('No structure or graph to initialise network from.')

        if structure is not None:
            # iterate over all pairs of atoms in check minimum distance
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
                self.add_node(i, species=atoms[i].species, image=False)
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
                                self.add_node(j+image_number, species=atoms[j].species, image=True)
                                self.add_edge(i, j+image_number, dist=dist)
                            else:
                                self.add_edge(i, j, dist=dist)

        elif graph is not None:
            for node, data in graph.nodes.data():
                self.add_node(node, species=data['species'], image=data['image'])
            for node_in, node_out, data in graph.edges.data():
                self.add_edge(node_in, node_out, dist=data['dist'])

    def get_strongly_connected_component_subgraphs(self):
        """ Return generator of strongly-connected subgraphs in CrystalGraph format. """
        return (CrystalGraph(graph=self.subgraph(c)) for c in nx.strongly_connected_components(self))

    def get_bonds_per_atom(self):
        num_bonds = 0
        for node_in in self.nodes():
            for node_out in self.nodes():
                if node_in == node_out:
                    continue
                if self.has_edge(node_in, node_out) and self.has_edge(node_out, node_in):
                    num_bonds += 1
        return num_bonds / self.number_of_nodes()

    @classmethod
    def are_graphs_the_same(graph_1, graph_2, bond_tolerance=1.1):
        raise NotImplementedError


def draw_network(structure, layout=None):
    import networkx as nx
    from matador.viz import ELEMENT_COLOURS
    import matplotlib.pyplot as plt
    try:
        network = structure.network
    except:
        network = structure
    if layout is None:
        pos = nx.spring_layout(network)
    else:
        pos = layout
    fig, ax = plt.subplots()
    elem_map = ELEMENT_COLOURS
    colours = [elem_map.get(data['species']) for node, data in network.nodes.data()]
    labels = {node: str(data['species']) for node, data in network.nodes.data()}
    nx.draw_networkx_nodes(network, pos, node_color=colours, edgecolors='black', linewidths=2, node_size=1000, ax=ax)
    nx.draw_networkx_edges(network, pos, linewidths=2, node_size=1000, ax=ax)
    nx.draw_networkx_labels(network, pos, labels=labels, ax=ax)
    plt.axis('off')
