# coding: utf-8
""" This file implements the Crystal class, a wrapper
to the raw dictionary stored in MongoDB that allows for validation,
manipulation and analysis of the lattice.
"""
from copy import deepcopy
from matador.similarity.pdf_similarity import PDF
import matador.utils.cell_utils as cell_utils
from matador.crystal_site import Site


class Crystal:
    """ Class that wraps the MongoDB document, providing useful
    interfaces for cell manipulation and validation.
    """
    def __init__(self, doc, voronoi=False):

        self._doc = deepcopy(doc)
        self.elems = sorted(list(set(self._doc['atom_types'])))
        self.sites = []
        self.construct_sites(voronoi=voronoi)

        # assume default value for symprec
        self._doc['space_group'] = {0.01: self._doc['space_group']}
        self._space_group = self._doc['space_group']

        for src in self._doc['source']:
            if src.endswith('.castep') or src.endswith('.res'):
                self.root_source = src.split('/')[-1] \
                                      .replace('.res', '').replace('.castep', '')

    def __getitem__(self, key):
        # if array-style access, e.g. crystal[3], return 3rd site object
        if isinstance(key, int):
            return self.sites[key]
        # if dict-style access, try to return key from doc/__dict__
        else:
            return self._doc[key]

    def __setitem__(self, key, item):
        self._doc[key] = item

    def __str__(self):
        return (("{root_source}: {formula}\n"
                 "{num_atoms:<3} atoms\n"
                 "(a, b, c) = {lattice[0][0]:4.4f} Å, {lattice[0][1]:4.4f} Å, {lattice[0][2]:4.4f} Å\n"
                 "(α, β, γ) = {lattice[1][0]:4.4f}° {lattice[1][1]:4.4f}° {lattice[1][2]:4.4f}°\n")
                .format(formula=self.formula, root_source=self.root_source, num_atoms=self.num_atoms, lattice=self.lattice_abc))

    def construct_sites(self, voronoi=False):
        for ind, species in enumerate(self.atom_types):
            position = self._doc['positions_frac'][ind]
            data = {}
            if 'chemical_shifts' in self._doc:
                data['shift'] = self._doc['chemical_shifts'][ind]
            if 'magnetic_shielding_tensor' in self._doc:
                data['shielding'] = self._doc['magnetic_shielding_tensor'][ind]
            if 'atomic_spins' in self._doc:
                data['spin'] = self._doc['atomic_spins'][ind]
            if 'voronoi_substructure' in self._doc:
                data['voronoi'] = self._doc['voronoi_substructure'][ind]
            elif voronoi:
                data['voronoi'] = self.voronoi_substructure[ind]

            self.sites.append(Site(species,
                                   position,
                                   self.lattice_cart,
                                   spin=data.get('spin'),
                                   voronoi_substructure=data.get('voronoi'),
                                   chemical_shift=data.get('shift'),
                                   magnetic_shielding=data.get('shielding')))

    @property
    def atom_types(self):
        return self._doc['atom_types']

    @property
    def num_atoms(self):
        return len(self.sites)

    @property
    def positions_frac(self):
        return self._doc['positions_frac']

    @property
    def stoichiometry(self):
        if '_stoichiometry' not in self.__dict__:
            self._stoichiometry = sorted(self._doc['stoichiometry'])
        return self._stoichiometry

    @property
    def formula(self):
        from matador.utils.chem_utils import get_formula_from_stoich
        return get_formula_from_stoich(self._doc['stoichiometry'], tex=False)

    @property
    def lattice_cart(self):
        if 'lattice_cart' not in self._doc:
            return cell_utils.abc2cart(self._doc['lattice_abc'])
        else:
            return self._doc['lattice_cart']

    @property
    def lattice_abc(self):
        if 'lattice_abc' not in self.__dict__:
            return cell_utils.cart2abc(self['lattice_cart'])
        else:
            return self['lattice_abc']

    @property
    def cell_volume(self):
        if 'cell_volume' not in self.__dict__:
            return cell_utils.cart2volume(self._doc['lattice_cart'])
        else:
            return self._doc['cell_volume']

    @property
    def voronoi_substructure(self):
        from matador.voronoi_interface import get_voronoi_substructure
        if '_voronoi_substructure' not in self.__dict__:
            self._voronoi_substructure = get_voronoi_substructure(self._doc)
            self._doc['voronoi_substructure'] = self._voronoi_substructure
        return self._voronoi_substructure

    @property
    def concentration(self):
        if '_concentration' in self.__dict__:
            return self._concentration
        concentration = []
        for species in self.stoichiometry:
            concentration.append(species[1])
        total = sum(concentration)
        for ind, _ in enumerate(concentration):
            concentration[ind] /= total

        self._concentration = concentration
        return self._concentration

    def get_ordered_concentration(self, elems):
        """ Return a padded and ordered concentration list based on the provided elements.

        Input:

            | elems: list(str), list of element symbols

        Returns:

            | concentration: list(float), concentrations of each element in elems.

        """
        concentration = []
        for species in elems:
            for ind, atom in enumerate(self.stoichiometry):
                if atom[0] == species:
                    concentration.append(atom[1])
                    break
                elif ind == len(self.stoichiometry)-1:
                    concentration.append(0)

        total = sum(concentration)
        for ind, _ in enumerate(concentration):
            concentration[ind] /= total

        return concentration

    def space_group(self, symprec=0.01):
        if symprec not in self._space_group:
            self._doc['space_group'][symprec] = cell_utils.get_spacegroup_spg(self._doc, symprec=symprec)
            self._space_group[symprec] = cell_utils.get_spacegroup_spg(self._doc, symprec=symprec)
        return self._space_group[symprec]

    @property
    def pdf(self, **kwargs):
        self._doc['pdf'] = PDF(self._doc, dr=0.01, gaussian_width=0.01, **kwargs)
        return self._doc['pdf']

    @property
    def coordination_stats(self):
        if '_coordination_stats' in self.__dict__:
            return self._coordination_stats

        coordination_lists = self.coordination_lists

        import numpy as np
        coordination_stats = {}
        for species in self.elems:
            coordination_stats[species] = {}
            for _species in self.elems:
                coordination_stats[species][_species] = {}
                coordination_stats[species][_species]['mean'] = np.mean(coordination_lists[species][_species])
                coordination_stats[species][_species]['median'] = np.median(coordination_lists[species][_species])
                coordination_stats[species][_species]['std'] = np.std(coordination_lists[species][_species])

        self._coordination_stats = coordination_stats
        return coordination_stats

    @property
    def coordination_lists(self):
        if '_coordination_lists' in self.__dict__:
            return self._coordination_lists
        else:
            coordination_lists = {}
            for species in self.elems:
                coordination_lists[species] = {}
                for _species in self.elems:
                    coordination_lists[species][_species] = []
            for site in self:
                for _species in self.elems:
                    if _species in site.coordination:
                        coordination_lists[site.species][_species].append(site.coordination[_species])
                    else:
                        coordination_lists[site.species][_species].append(0)

            self._coordination_lists = coordination_lists

            return coordination_lists

    @property
    def unique_sites(self):
        from matador.similarity.voronoi_similarity import get_unique_sites
        # from collections import defaultdict
        # substructures = []
        # substructure_dict = defaultdict(list)
        # for atom in self.voronoi_substructure:
            # substructure_dict[atom[0]].append(atom[1])
        # print(substructure_dict)
        # array = create_site_array(substructure_dict)
        get_unique_sites(self._doc)
        return self._doc['similar_sites']

    @property
    def network(self):
        from matador.network import construct_network
        if '_network' not in self.__dict__:
            self._network = construct_network(self)
        return self._network

    @property
    def network_stats(self):
        from collections import defaultdict
        coordination = defaultdict(list)
        for node, data in self.network.nodes.data():
            num_edges = len(self.network.edges(node))
            coordination[data['species']].append(num_edges)
        return coordination

    def draw_network(self, layout=None):
        from matador.network import draw_network
        draw_network(self, layout=layout)
