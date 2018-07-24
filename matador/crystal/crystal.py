# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements the Crystal class, a wrapper
to the raw dictionary stored in MongoDB that allows for validation,
manipulation and analysis of the lattice.

"""

from copy import deepcopy
import matador.utils.cell_utils as cell_utils
from matador.similarity.pdf_similarity import PDF
from matador.crystal.site import Site


class Crystal:
    """ Class that wraps the MongoDB document, providing useful
    interfaces for cell manipulation and validation.
    """

    def __init__(self, doc, voronoi=False, network_kwargs=None):
        """ Initialise Crystal object from matador document with Site list
        and any additional abstractions, e.g. voronoi or CrystalGraph.

        Parameters:
            doc (dict): matador document containing structural information

        Keyword Arguments:
           voronoi (bool): whether to compute Voronoi substructure for each site
           network_kwargs (dict): keywords to pass to the CrystalGraph initialiser

        """

        self._doc = deepcopy(doc)
        self.elems = sorted(list(set(self._doc['atom_types'])))
        self.sites = []
        self._construct_sites(voronoi=voronoi)

        if not any(key in self._doc for key in ['lattice_cart', 'lattice_abc']):
            raise RuntimeError('No lattice information found, cannot create Crystal.')

        if network_kwargs is not None:
            self._network_kwargs = network_kwargs
        else:
            self._network_kwargs = {}

        # assign attributes for later
        self._coordination_lists = None
        self._coordination_stats = None
        self._network = None
        self._bond_lengths = None

        # assume default value for symprec
        if 'space_group' in self._doc:
            self._space_group = {0.01: self._doc['space_group']}
        else:
            self._space_group = {}

        # set root source to structure filename
        from matador.utils.chem_utils import get_root_source
        try:
            self.root_source = get_root_source(self._doc['source'])
        except RuntimeError:
            self.root_source = 'xxx'

    def __getitem__(self, key):
        # if array-style access, e.g. crystal[3], return 3rd site object
        if isinstance(key, int):
            return self.sites[key]

        # otherwise, if dict-style access, try to return key from doc/__dict__
        return self._doc[key]

    def __setitem__(self, key, item):
        self._doc[key] = item

    def __contains__(self, key):
        if key in self._doc:
            return True
        return False

    def __str__(self):
        repr_string = "{root_source}: {formula}\n".format(root_source=self.root_source, formula=self.formula)
        repr_string += "{num_atoms:<3} atoms. {space_group:<8}\n".format(num_atoms=self.num_atoms,
                                                                         space_group=self.space_group)

        if 'formation_enthalpy_per_atom' in self._doc:
            repr_string += ("Formation enthalpy = {:6.6f} eV/atom\n".format(self._doc['formation_enthalpy_per_atom']))

        repr_string += (
            "(a, b, c) = {lattice[0][0]:4.4f} Å, {lattice[0][1]:4.4f} Å, {lattice[0][2]:4.4f} Å\n"
            "(α, β, γ) = {lattice[1][0]:4.4f}° {lattice[1][1]:4.4f}° {lattice[1][2]:4.4f}°\n"
            .format(lattice=self.lattice_abc))

        return repr_string

    def get(self, key):
        """ Overload dictionary.get() method.

        Parameters:
            key (str): key to try and obtain.

        Returns:
            doc[key] if it exists, else None.

        """
        return self._doc.get(key)

    def _construct_sites(self, voronoi=False):
        """ Constructs the list of Site objects stored in self.sites.

        Keyword arguments:
            voronoi (bool): whether to calculate the Voronoi substructure
                of each site.

        """
        for ind, species in enumerate(self.atom_types):
            position = self._doc['positions_frac'][ind]
            site_data = {}
            if 'site_occupancy' in self._doc:
                if len(self._doc['site_occupancy']) == len(self._doc['atom_types']):
                    site_data['site_occupancy'] = self._doc['site_occupancy'][ind]
            if 'chemical_shifts' in self._doc:
                if len(self._doc['chemical_shifts']) == len(self._doc['atom_types']):
                    site_data['magres_shift'] = self._doc['chemical_shifts'][ind]
            if 'magnetic_shielding_tensor' in self._doc:
                if len(self._doc['magnetic_shielding_tensor']) == len(self._doc['atom_types']):
                    site_data['magres_shielding'] = self._doc['magnetic_shielding_tensor'][ind]
            if 'chemical_shift_anisos' in self._doc:
                site_data['magres_aniso'] = self._doc['chemical_shift_anisos'][ind]
            if 'chemical_shift_asymmetries' in self._doc:
                site_data['magres_asymm'] = self._doc['chemical_shift_asymmetries'][ind]
            if 'atomic_spins' in self._doc:
                site_data['spin'] = self._doc['atomic_spins'][ind]
            if 'voronoi_substructure' in self._doc:
                site_data['voronoi_substructure'] = self._doc['voronoi_substructure'][ind]
            elif voronoi:
                site_data['voronoi_substructure'] = self.voronoi_substructure[ind]

            self.sites.append(Site(species, position, self.lattice_cart, **site_data))

    @property
    def atom_types(self):
        """ Return list of atom types. """
        return self._doc['atom_types']

    @property
    def num_atoms(self):
        """ Return number of atoms in structure. """
        return len(self.sites)

    @property
    def positions_frac(self):
        """ Return list of fractional positions. """
        return self._doc['positions_frac']

    @property
    def stoichiometry(self):
        """ Return stoichiometry in matador format: a list of
        two-member lists containing element symbol and number
        of atoms per formula unit, sorted in alphabetical order
        by element symbol).

        """
        if 'stoichiometry' not in self._doc:
            from matador.utils.chem_utils import get_stoich
            self._doc['stoichiometry'] = get_stoich(self.atom_types)
        return self._doc['stoichiometry']

    @property
    def formula(self):
        """ Returns chemical formula of structure. """
        from matador.utils.chem_utils import get_formula_from_stoich
        return get_formula_from_stoich(self._doc['stoichiometry'], tex=False)

    @property
    def lattice_cart(self):
        """ Returns Cartesian lattice vectors. """
        if 'lattice_cart' not in self._doc:
            self._doc['lattice_cart'] = cell_utils.abc2cart(self._doc['lattice_abc'])
        return self._doc['lattice_cart']

    @property
    def lattice_abc(self):
        """ Returns lattice parameters. """
        if 'lattice_abc' not in self._doc:
            self._doc['lattice_abc'] = cell_utils.cart2abc(self['lattice_cart'])
        return self._doc['lattice_abc']

    @property
    def cell_volume(self):
        """ Returns cell volume in Å³. """
        if 'cell_volume' not in self.__dict__ or 'cell_volume' not in self._doc:
            self._doc['cell_volume'] = cell_utils.cart2volume(self._doc['lattice_cart'])
        return self._doc['cell_volume']

    @property
    def bond_lengths(self):
        """ Returns a list of ((species_A, species_B), bond_length)),
        sorted by bond length, computed from the network structure of
        the crystal (i.e. first coordination sphere).

        """
        if self._bond_lengths is None:
            self._bond_lengths = []
            for i, j, data in self.network.edges.data():
                self._bond_lengths.append(((self[i].species, self[j].species), data['dist']))
            self._bond_lengths = sorted(self._bond_lengths, key=lambda bond: bond[1])
        return self._bond_lengths

    @property
    def voronoi_substructure(self):
        """ Returns, or calculates  if not present, the Voronoi
        substructure of crystal.

        """
        if 'voronoi_substructure' not in self._doc:
            from matador.plugins.voronoi_interface.voronoi_interface import get_voronoi_substructure
            self._doc['voronoi_substructure'] = get_voronoi_substructure(self._doc)
        return self._doc['voronoi_substructure']

    @property
    def concentration(self):
        """ Returns the an array of floats storing the concentration of
        each element in the structure.

        """
        if 'concentration' in self._doc:
            return self._doc['concentration']

        concentration = []
        for species in self.stoichiometry:
            concentration.append(species[1])
        total = sum(concentration)
        for ind, _ in enumerate(concentration):
            concentration[ind] /= total

        self._doc['concentration'] = concentration
        return self._doc['concentration']

    def get_ordered_concentration(self, elems):
        """ Return a padded and ordered concentration list based on the
        provided elements.

        Parameters:
            elems (list of str): list of element symbols

        Returns:
            concentration: (list of float), concentrations of each
                element in elems.

        """
        concentration = []
        for species in elems:
            for ind, atom in enumerate(self.stoichiometry):
                if atom[0] == species:
                    concentration.append(atom[1])
                    break
                elif ind == len(self.stoichiometry) - 1:
                    concentration.append(0)

        total = sum(concentration)
        for ind, _ in enumerate(concentration):
            concentration[ind] /= total

        return concentration

    @property
    def space_group(self):
        """ Return the space group symbol at the last-used symprec. """
        return self.get_space_group(symprec=self._doc.get('symprec', 0.01))

    def get_space_group(self, symprec=0.01):
        """ Return the space group of the structure at the desired
        symprec. Stores the space group in a dictionary
        `self._space_group` under symprec keys. Updates
        `self._doc['space_group']` and `self._doc['symprec']` with the
        last value calculated.

        Keyword arguments:
            symprec (float): spglib symmetry tolerance.

        """
        if symprec not in self._space_group:
            self._doc['space_group'] = cell_utils.get_spacegroup_spg(self._doc, symprec=symprec)
            self._doc['symprec'] = symprec
            self._space_group[symprec] = self._doc['space_group']
        return self._space_group[symprec]

    @property
    def pdf(self, **kwargs):
        """ Returns a PDF object (pair distribution function) for the
        structure, calculated with default PDF settings, unless kwargs
        are passed.

        """
        self._doc['pdf'] = PDF(self._doc, **kwargs)
        return self._doc['pdf']

    @property
    def coordination_stats(self):
        """ Returns stastistics on the coordination of each site, as
        computed from Voronoi decomposition.

        """
        if self._coordination_stats is None:
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

        return self._coordination_stats

    @property
    def coordination_lists(self):
        """ Returns a dictionary containing pairs of elements and the
        list of coordination numbers per site.

        """
        if self._coordination_lists is None:
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

        return self._coordination_lists

    @property
    def unique_sites(self):
        """ Return unique sites using Voronoi decomposition. """
        from matador.similarity.voronoi_similarity import get_unique_sites
        get_unique_sites(self._doc)
        return self._doc['similar_sites']

    @property
    def network(self):
        """ Returns/constructs a CrystalGraph object of the structure. """
        from matador.crystal.network import CrystalGraph
        if self._network is None:
            self._network = CrystalGraph(self, **self._network_kwargs)
        return self._network

    @property
    def network_stats(self):
        """ Return network-calculated coordnation stats. """
        from collections import defaultdict
        coordination = defaultdict(list)
        for node, data in self.network.nodes.data():
            num_edges = len(self.network.edges(node))
            coordination[data['species']].append(num_edges)
        return coordination

    def draw_network(self, layout=None):
        """ Draw the CrystalGraph network. """
        from matador.network import draw_network
        draw_network(self, layout=layout)
