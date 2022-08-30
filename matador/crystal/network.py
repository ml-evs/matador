""" This file implements turning matador Crystal objects
into CrystalGraph objects.
"""
import networkx as nx
import numpy as np
import itertools

EPS = 1e-12


class CrystalGraph(nx.MultiDiGraph):
    def __init__(
        self,
        structure=None,
        graph=None,
        coordination_cutoff=1.1,
        bond_tolerance=1e20,
        num_images=1,
        debug=False,
        separate_images=False,
        delete_one_way_bonds=False,
        max_bond_length=5,
    ):
        """Create networkx.MultiDiGraph object with extra functionality for atomic networks.

        Keyword Arguments:
            structure (matador.Crystal):  crystal structure to network-ify
            graph (nx.MultiDiGraph): initialise from graph
            coordination_cutoff (float) : max multiplier of first
                coordination sphere for edge drawing
            num_images (int): number of periodic images to include in
                each direction
            separate_images (bool): whether or not to include image
                atoms as new nodes

        """

        super().__init__()

        if graph is None and structure is None:
            raise RuntimeError("No structure or graph to initialise network from.")

        if structure is not None:
            atoms = structure.sites
            num_atoms = len(atoms)
            element_bonds = {}
            images = list(
                itertools.product(range(-num_images, num_images + 1), repeat=3)
            )
            image_number = 0

            # now loop over pairs of atoms and decide whether to draw an edge
            from matador.utils.cell_utils import calc_pairwise_distances_pbc

            distances = calc_pairwise_distances_pbc(
                structure.positions_abs,
                images,
                structure.lattice_cart,
                max_bond_length,
                compress=False,
                debug=True,
            )

            # first over loop all pairs to find the minimum distance between all species pairs
            # and the minimum distance for each atom
            for i, atom in enumerate(atoms):
                self.add_node(i, species=atom.species)

            min_dists = [1e20 for atom in atoms]
            for index in np.where(~distances.mask)[0]:
                image_index = int(index / num_atoms**2)
                i = int((index - image_index * num_atoms**2) / num_atoms)
                j = int((index - image_index * num_atoms**2) % num_atoms)
                atom = atoms[i]
                other_atom = atoms[j]
                if i == j and np.linalg.norm(images[image_index]) <= EPS:
                    continue
                dist = distances[index]
                pair_key = tuple(sorted([atom.species, other_atom.species]))
                if pair_key not in element_bonds or element_bonds[pair_key] > dist:
                    element_bonds[pair_key] = dist
                # find the closest image of an atom to i
                if dist < min_dists[i]:
                    min_dists[i] = dist

            if debug:
                print(min_dists)
                print(element_bonds)

            for index in np.where(~distances.mask)[0]:
                image_index = int(index / num_atoms**2)
                i = int((index - image_index * num_atoms**2) / num_atoms)
                j = int((index - image_index * num_atoms**2) % num_atoms)
                atom = atoms[i]
                other_atom = atoms[j]
                min_dist = min_dists[i]
                if i == j and np.linalg.norm(images[image_index]) <= EPS:
                    continue
                dist = distances[index]
                pair_key = tuple(sorted([atom.species, other_atom.species]))
                if (
                    dist <= min_dist * coordination_cutoff
                    and dist <= element_bonds[pair_key] * bond_tolerance
                ):
                    if separate_images and all(
                        [val <= 0 + 1e-8 for val in images[image_index]]
                    ):
                        image_number += 1
                        self.add_node(j + image_number, species=atoms[j].species)
                        self.add_edge(i, j + image_number, dist=dist)
                        self.add_edge(j + image_number, i, dist=dist)
                    else:
                        is_image = np.linalg.norm(images[image_index]) > EPS
                        self.add_edge(i, j, dist=dist, image=is_image)

        elif graph is not None:
            for node, data in graph.nodes.data():
                self.add_node(
                    node, species=data["species"], image=data.get("image", False)
                )
            for node_in, node_out, data in graph.edges.data():
                if not delete_one_way_bonds or (node_out, node_in) in graph.edges():
                    self.add_edge(
                        node_in,
                        node_out,
                        dist=data.get("dist", 0),
                        image=data.get("image", False),
                    )

    def get_strongly_connected_component_subgraphs(self, delete_one_way_bonds=True):
        """Return generator of strongly-connected subgraphs in CrystalGraph format."""
        return (
            CrystalGraph(
                graph=self.subgraph(c), delete_one_way_bonds=delete_one_way_bonds
            )
            for c in nx.strongly_connected_components(self)
        )

    def get_communities(self, graph=None, **louvain_kwargs):
        """Return list of community subgraphs in CrystalGraph format."""
        import community as louvain

        if graph is None:
            graph = self

        if graph.is_directed():
            undirected_graph = self.remove_directionality(graph=graph)

        partition = louvain.best_partition(undirected_graph, **louvain_kwargs)
        size = len(set(partition.values()))
        subgraphs = [nx.MultiGraph() for i in range(size)]

        for node in partition:
            subgraphs[partition[node]].add_node(
                node,
                species=list(self.nodes(data=True))[list(self.nodes()).index(node)][1][
                    "species"
                ],
            )
        for edge in self.edges():
            if partition[edge[0]] == partition[edge[1]]:
                subgraphs[partition[edge[0]]].add_edge(edge[0], edge[1])

        subgraphs = [CrystalGraph(graph=sg) for sg in subgraphs]

        return subgraphs, partition

    def remove_directionality(self, graph=None):
        if graph is None:
            graph = self
        import networkx as nx

        undirected_graph = nx.MultiGraph()
        for node in graph.nodes(data=True):
            undirected_graph.add_node(node[0], species=node[1]["species"])
        for edge in graph.edges(data=True):
            if (edge[1], edge[0]) not in undirected_graph.edges():
                undirected_graph.add_edge(
                    edge[0],
                    edge[1],
                    dist=edge[2].get("dist", 0),
                    image=edge[2].get("image", False),
                )

        return undirected_graph

    def set_unique_subgraphs(self, method="community"):
        """Filter strongly connected component subgraphs for isomorphism with others inside
        CrystalGraph. Sets self.unique_subgraph to a set of such subgraphs.
        """
        if method == "community":
            self.unique_subgraphs = get_unique_subgraphs(self.get_communities())
        elif method == "strongly_connected":
            self.unique_subgraphs = get_unique_subgraphs(
                self.get_strongly_connected_component_subgraphs()
            )
        elif method == "both":
            strong_subgraphs = self.get_strongly_connected_component_subgraphs()
            community_subgraphs = []
            for sg in strong_subgraphs:
                community_subgraphs.extend(sg.get_communities())

            self.unique_subgraphs = get_unique_subgraphs(community_subgraphs)

    def get_bonds_per_atom(self):
        num_bonds = 0
        for node_in in self.nodes():
            for node_out in self.nodes():
                if node_in == node_out:
                    continue
                if self.has_edge(node_in, node_out) and self.has_edge(
                    node_out, node_in
                ):
                    num_bonds += 1
        return num_bonds / self.number_of_nodes()


def node_match(n1, n2):
    return n1["species"] == n2["species"]


def get_unique_subgraphs(subgraphs):
    """Filter strongly connected component subgraphs for isomorphism with others.

    Input:

        | subgraphs: list(CrystalGraph), list of subgraph objects to filter

    Returns:

        | unique_subgraphs: set(CrystalGraph), set of unique subgraphs

    """
    unique_subgraphs = set()
    for subgraph in subgraphs:
        if not any(
            [
                are_graphs_the_same(subgraph, other_subgraph)
                for other_subgraph in unique_subgraphs
            ]
        ):
            unique_subgraphs.add(subgraph)
    return unique_subgraphs


def are_graphs_the_same(g1, g2, edge_match=None):
    if edge_match is None:

        def edge_match(e1, e2):
            atol = 0.1
            rtol = 0.05
            return (
                abs(e1[0]["dist"] - e2[0]["dist"]) <= atol + rtol * e2[0]["dist"]
                and e1[0]["image"] == e2[0]["image"]
            )

    return nx.is_isomorphic(
        g1,
        g2,
        node_match=lambda n1, n2: n1["species"] == n2["species"],
        edge_match=edge_match,
    )


def draw_network(
    structure,
    layout=None,
    edge_labels=False,
    node_index=False,
    curved_edges=True,
    node_colour="elem",
    partition=None,
    ax=None,
):
    import networkx as nx
    from matador.utils.viz_utils import get_element_colours
    import matplotlib.pyplot as plt

    element_colours = get_element_colours()
    try:
        network = structure.network
    except Exception:
        network = structure
    if layout is None:
        pos = nx.spring_layout(network)
    else:
        pos = layout

    if ax is None:
        fig, ax = plt.subplots()

    if node_colour == "degree":
        coords = list(set(dict(network.degree).values()))
        cmap = plt.cm.get_cmap("Dark2", len(coords)).colors
        colours = [cmap[coords.index(network.degree[node])] for node in network.nodes()]

    elif node_colour == "partition" and partition is not None:
        num_partitions = len(set(partition.values()))
        cmap = plt.cm.get_cmap("Dark2", num_partitions).colors
        colours = [cmap[partition[node]] for node in network.nodes()]

    else:
        elem_map = element_colours
        colours = [elem_map.get(data["species"]) for node, data in network.nodes.data()]

    if node_index:
        labels = {
            node: "{} \\#{}".format(data["species"], node)
            for node, data in network.nodes.data()
        }
    else:
        labels = {node: str(data["species"]) for node, data in network.nodes.data()}

    edge_colours = []
    for edge in network.edges(data=True):
        if edge[2].get("image", True):
            edge_colours.append("grey")
        else:
            edge_colours.append("black")

    nx.draw_networkx_nodes(
        network,
        pos,
        node_color=colours,
        edgecolors="black",
        linewidths=2,
        node_size=1000,
        ax=ax,
    )
    nx.draw_networkx_edges(
        network, pos, edge_color=edge_colours, width=2, node_size=1000, ax=ax
    )
    if edge_labels:
        edge_weight = dict()
        for edge in network.edges(data=True):
            # data = edge[2]
            edge = (edge[0], edge[1])
            if edge not in edge_weight and (edge[1], edge[0]) not in edge_weight:
                edge_weight[edge] = 1
            else:
                if edge in edge_weight:
                    edge_weight[edge] += 1
                else:
                    edge_weight[(edge[1], edge[0])] += 1
        edge_label_dict = edge_weight
        nx.draw_networkx_edge_labels(network, pos, edge_labels=edge_label_dict, ax=ax)
    nx.draw_networkx_labels(network, pos, labels=labels, ax=ax)
    plt.axis("off")
