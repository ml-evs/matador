# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements the Crystal class, a wrapper
to the raw dictionary stored in MongoDB that allows for validation,
manipulation and analysis of the lattice.

"""

from copy import deepcopy
from matador.utils import cell_utils
from matador.orm.orm import DataContainer
from matador.crystal.crystal_site import Site
from matador.utils.chem_utils import get_concentration
from matador.utils.cell_utils import real2recip


class UnitCell:
    """ This class describes a 3D periodic unit cell by its
    Cartesian lattice vectors or lattice parameters, in Å.

    """

    _lattice_abc = None
    _lattice_cart = None
    _volume = None

    def __init__(self, lattice):
        """ Initialise the cell from either Cartesian lattice vectors
        [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]], or lattice
        parameters [[a, b, c], [alpha, beta, gamma]].

        Parameters:
            lattice (:obj:`list` or numpy.ndarray): either Cartesian
                lattice vectors or lattice parameters, stored as lists
                or a numpy arrays.

        """

        self._volume = None
        if len(lattice) == 3:
            if any(len(vec) != 3 for vec in lattice):
                raise RuntimeError('Unable to cast {} into lattice_cart'.format(lattice))
            self.lattice_cart = lattice
        elif len(lattice) == 2:
            if any(len(vec) != 3 for vec in lattice):
                raise RuntimeError('Unable to cast {} into lattice_abc'.format(lattice))
            self.lattice_abc = lattice
        else:
            raise RuntimeError("Unable to create UnitCell from lattice {}".format(lattice))

    @property
    def lattice_cart(self):
        """ The Cartesian lattice vectors as a tuple. """
        return self._lattice_cart

    @lattice_cart.setter
    def lattice_cart(self, new_lattice):
        self._lattice_cart = tuple(tuple(vec) for vec in new_lattice)
        self._lattice_abc = tuple(tuple(elem) for elem in cell_utils.cart2abc(self._lattice_cart))
        self._volume = None

    @property
    def lattice_abc(self):
        """ Lattice parameters as a tuple. """
        return self._lattice_abc

    @lattice_abc.setter
    def lattice_abc(self, new_lattice):
        self._lattice_abc = tuple(tuple(vec) for vec in new_lattice)
        self._lattice_cart = tuple(tuple(vec) for vec in cell_utils.abc2cart(self._lattice_abc))
        self._volume = None

    @property
    def lengths(self):
        """ Lattice vector lengths. """
        return self._lattice_abc[0]

    @lengths.setter
    def lengths(self, new_lengths):
        if len(new_lengths) != 3:
            raise RuntimeError('Expected list of 3 floats for cell vector lengths, received {}'
                               .format(new_lengths))
        self.lattice_abc = [new_lengths, self._lattice_abc[1]]

    @property
    def recip_lattice_cart(self):
        return real2recip(self.lattice_cart)

    @property
    def angles(self):
        """ Lattice vector angles. """
        return self._lattice_abc[1]

    @angles.setter
    def angles(self, new_angles):
        if len(new_angles) != 3:
            raise RuntimeError('Expected list of 3 floats for cell vector angles, received {}'
                               .format(new_angles))
        self.lattice_abc = [self._lattice_abc[0], new_angles]

    @property
    def volume(self):
        """ The cell volume in Å³. """
        if not self._volume:
            self._volume = cell_utils.cart2volume(self._lattice_cart)
        return self._volume


class Crystal(DataContainer):
    """ Class that wraps the MongoDB document, providing useful
    interfaces for cell manipulation and validation.

    Attributes:
        elems (:obj:`list` of :obj:`str`): list of present elements in the crystal,
            sorted in alphabetical order.
        cell (:obj:`UnitCell`): the unit cell constructed from lattice vectors.

    """

    @staticmethod
    def _validate_doc(doc):
        if not any(key in doc for key in ['lattice_cart', 'lattice_abc']):
            raise RuntimeError('No lattice information found, cannot create Crystal.')
        if 'atom_types' not in doc:
            raise RuntimeError('No species information found `"atom_types"`, cannot create Crystal.')
        if not any(key in doc for key in ['positions_frac', 'positions_abs']):
            raise RuntimeError('No position information found `"positions_frac"/"positions_abs"`, cannot create Crystal.')

    def __init__(self, doc, voronoi=False, network_kwargs=None):
        """ Initialise Crystal object from matador document with Site list
        and any additional abstractions, e.g. voronoi or CrystalGraph.

        Parameters:
            doc (dict): document containing structural information, minimal requirement
                is for `atom_types`, one of `lattice_abc` or `lattice_cart`, and one of
                `positions_frac` or `positions_abs` to be present.

        Keyword Arguments:
           voronoi (bool): whether to compute Voronoi substructure for each site
           network_kwargs (dict): keywords to pass to the CrystalGraph initialiser

        """

        self._validate_doc(doc)
        if isinstance(doc, Crystal):
            doc = deepcopy(doc._data)

        super().__init__(doc)

        self.elems = sorted(list(set(self._data['atom_types'])))
        self.sites = []

        # use lattice_cart to construct cell if present, otherwise abc
        self.cell = UnitCell(doc.get('lattice_cart', doc.get('lattice_abc')))

        self._construct_sites(voronoi=voronoi)

        if network_kwargs is not None:
            self._network_kwargs = network_kwargs
        else:
            self._network_kwargs = {}

        # assign attributes for later
        self._coordination_lists = None
        self._coordination_stats = None
        self._network = None
        self._bond_lengths = None
        self._bonding_stats = None

        # assume default value for symprec
        if 'space_group' in self._data:
            self._space_group = {0.01: self._data['space_group']}
        else:
            self._space_group = {}

    def __getitem__(self, key):
        """ If integer key is requested, return index into site array. """
        if isinstance(key, int):
            return self.sites[key]
        return super().__getitem__(key)

    def __str__(self):
        repr_string = "{root_source}: {formula}\n".format(root_source=self.root_source, formula=self.formula)
        repr_string += "{num_atoms:<3} atoms. {space_group:<8}\n".format(num_atoms=self.num_atoms,
                                                                         space_group=self.space_group)

        if 'formation_enthalpy_per_atom' in self._data:
            repr_string += ("Formation enthalpy = {:6.6f} eV/atom\n".format(self._data['formation_enthalpy_per_atom']))

        repr_string += (
            "(a, b, c) = {lattice[0][0]:4.4f} Å, {lattice[0][1]:4.4f} Å, {lattice[0][2]:4.4f} Å\n"
            "(α, β, γ) = {lattice[1][0]:4.4f}° {lattice[1][1]:4.4f}° {lattice[1][2]:4.4f}°\n"
            .format(lattice=self.cell.lattice_abc))

        return repr_string

    def set_positions(self, new_positions, fractional=True):
        if len(new_positions) != self.num_atoms:
            raise RuntimeError('Cannot change size of positions array!')
        if fractional:
            self._data['positions_frac'] = new_positions
            self._data.pop('positions_abs', None)
        else:
            self._data['positions_abs'] = new_positions
            self._data.pop('positions_frac', None)

        self._construct_sites()

    def _construct_sites(self, voronoi=False):
        """ Constructs the list of Site objects stored in self.sites.

        Keyword arguments:
            voronoi (bool): whether to calculate the Voronoi substructure
                of each site.

        """
        self.sites = []
        for ind, species in enumerate(self.atom_types):
            position = self.positions_frac[ind]
            site_data = {}
            if 'site_occupancy' in self._data:
                if len(self._data['site_occupancy']) == len(self._data['atom_types']):
                    site_data['site_occupancy'] = self._data['site_occupancy'][ind]
            if 'chemical_shifts' in self._data:
                if len(self._data['chemical_shifts']) == len(self._data['atom_types']):
                    site_data['magres_shift'] = self._data['chemical_shifts'][ind]
            if 'magnetic_shielding_tensor' in self._data:
                if len(self._data['magnetic_shielding_tensor']) == len(self._data['atom_types']):
                    site_data['magres_shielding'] = self._data['magnetic_shielding_tensor'][ind]
            if 'chemical_shift_anisos' in self._data:
                site_data['magres_aniso'] = self._data['chemical_shift_anisos'][ind]
            if 'chemical_shift_asymmetries' in self._data:
                site_data['magres_asymm'] = self._data['chemical_shift_asymmetries'][ind]
            if 'atomic_spins' in self._data:
                site_data['spin'] = self._data['atomic_spins'][ind]
            if 'voronoi_substructure' in self._data:
                site_data['voronoi_substructure'] = self._data['voronoi_substructure'][ind]
            elif voronoi:
                site_data['voronoi_substructure'] = self.voronoi_substructure[ind]

            self.sites.append(Site(species, position, self.cell.lattice_cart, **site_data))

    @property
    def atom_types(self):
        """ Return list of atom types. """
        return self._data['atom_types']

    @property
    def num_atoms(self):
        """ Return number of atoms in structure. """
        return len(self.sites)

    @property
    def num_elements(self):
        """ Return number of species in the structure. """
        return len(self.elems)

    @property
    def positions_frac(self):
        """ Return list of fractional positions. """
        from matador.utils.cell_utils import cart2frac
        if 'positions_frac' not in self._data:
            self._data['positions_frac'] = cart2frac(self.cell.lattice_cart, self.positions_abs)
        return self._data['positions_frac']

    @property
    def positions_abs(self):
        """ Return list of absolute Cartesian positions. """
        from matador.utils.cell_utils import frac2cart
        if 'positions_abs' not in self._data:
            self._data['positions_abs'] = frac2cart(self.cell.lattice_cart, self.positions_frac)
        return self._data['positions_abs']

    @property
    def site_occupancies(self):
        """ Return a list of site occupancies. """
        if 'site_occupancies' not in self._data:
            self._data['site_occupancies'] = [site.occupancy for site in self]
        return self._data['site_occupancies']

    @property
    def stoichiometry(self):
        """ Return stoichiometry in matador format: a list of
        two-member lists containing element symbol and number
        of atoms per formula unit, sorted in alphabetical order
        by element symbol).

        """
        if 'stoichiometry' not in self._data:
            from matador.utils.chem_utils import get_stoich
            self._data['stoichiometry'] = get_stoich(self.atom_types)
        return self._data['stoichiometry']

    @property
    def concentration(self):
        """ Return concentration of each species in stoichiometry. """
        if 'concentration' not in self._data:
            self._data['concentration'] = get_concentration(
                self.stoichiometry,
                [elem[0] for elem in self.stoichiometry],
                include_end=True
            )
        return self._data['concentration']

    @property
    def formula(self):
        """ Returns chemical formula of structure. """
        from matador.utils.chem_utils import get_formula_from_stoich
        return get_formula_from_stoich(self.stoichiometry, tex=False)

    @property
    def formula_tex(self):
        """ Returns chemical formula of structure in LaTeX format. """
        from matador.utils.chem_utils import get_formula_from_stoich
        return get_formula_from_stoich(self.stoichiometry, tex=True)

    @property
    def cell_volume(self):
        """ Returns cell volume in Å³. """
        return self.cell.volume

    @property
    def lattice_cart(self):
        """ The Cartesian lattice vectors as a tuple. """
        return self.cell.lattice_cart

    @property
    def lattice_abc(self):
        """ Lattice parameters as a tuple. """
        return self.cell.lattice_abc

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
        if 'voronoi_substructure' not in self._data:
            from matador.plugins.voronoi_interface.voronoi_interface import get_voronoi_substructure
            self._data['voronoi_substructure'] = get_voronoi_substructure(self._data)
        return self._data['voronoi_substructure']

    @property
    def space_group(self):
        """ Return the space group symbol at the last-used symprec. """
        return self.get_space_group(symprec=self._data.get('symprec', 0.01))

    def get_space_group(self, symprec=0.01):
        """ Return the space group of the structure at the desired
        symprec. Stores the space group in a dictionary
        `self._space_group` under symprec keys. Updates
        `self._data['space_group']` and `self._data['symprec']` with the
        last value calculated.

        Keyword arguments:
            symprec (float): spglib symmetry tolerance.

        """
        if symprec not in self._space_group:
            self._data['space_group'] = cell_utils.get_spacegroup_spg(self._data, symprec=symprec, check_occ=False)
            self._data['symprec'] = symprec
            self._space_group[symprec] = self._data['space_group']
        return self._space_group[symprec]

    @property
    def pdf(self):
        """ Returns a PDF object (pair distribution function) for the
        structure, calculated with default PDF settings.

        """
        from matador.fingerprints.pdf import PDF
        if 'pdf' not in self._data:
            self._data['pdf'] = PDF(self._data, label=self.formula_tex)
        return self._data['pdf']

    @pdf.setter
    def pdf(self, pdf):
        """ Set the PDF to the given PDF object (or None). """
        from matador.fingerprints.pdf import PDF
        if isinstance(pdf, PDF) or pdf is None:
            self._data['pdf'] = pdf

    def calculate_pdf(self, **kwargs):
        """ Calculate and set the PDF with the passed parameters. """
        from matador.fingerprints.pdf import PDF
        if 'pdf' not in self._data:
            self._data['pdf'] = PDF(self._data, label=self.formula_tex, **kwargs)
        return self._data['pdf']

    @property
    def pxrd(self):
        """ Returns a PXRD object (powder xray diffraction) containing
        the XRD pattern for the structure.

        """
        from matador.fingerprints.pxrd import PXRD
        if 'pxrd' not in self._data:
            self._data['pxrd'] = PXRD(self._data)
        return self._data['pxrd']

    @pxrd.setter
    def pxrd(self, pxrd):
        """ Set the PXRD to the given PXRD object (or None). """
        from matador.fingerprints.pxrd import PXRD
        if isinstance(pxrd, PXRD) or pxrd is None:
            self._data['pxrd'] = pxrd

    def calculate_pxrd(self, **kwargs):
        """ Compute and set PXRD with the passed parameters. """
        from matador.fingerprints.pxrd import PXRD
        if 'pxrd' not in self._data:
            self._data['pxrd'] = PXRD(self._data, **kwargs)
        return self._data['pxrd']

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
    def ase_atoms(self):
        """ Returns an ASE Atoms representation of
        the crystal.

        """
        from matador.utils.ase_utils import doc2ase
        return doc2ase(self, add_keys_to_info=False)

    @property
    def pmg_structure(self):
        """ Returns the pymatgen structure representation of
        the crystal.

        """
        from matador.utils.pmg_utils import doc2pmg
        return doc2pmg(self)

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
        from matador.plugins.voronoi_interface.voronoi_similarity import get_unique_sites
        get_unique_sites(self._data)
        return self._data['similar_sites']

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

    @property
    def bonding_stats(self):
        """ Return network-calculated bonding stats.

        Returns:
            dict: sorted dictionary with root atom as keys and bond
                information as values.

        """
        if self._bonding_stats is None:
            from collections import defaultdict
            bonding_dict = defaultdict(dict)
            network = self.network
            for node in network.nodes:
                bonding_dict[node] = {'species': self.sites[node].species, 'position': self.sites[node].coords, 'bonds': []}
            bonds = set()
            for data in network.edges.data():
                atom_1 = data[0]
                atom_2 = data[1]
                pair = tuple(sorted([atom_1, atom_2]))
                if pair in bonds:
                    continue
                else:
                    bonds.add(pair)

                site_1 = self.sites[atom_1]
                site_2 = self.sites[atom_2]

                bond_length = data[2]['dist']
                is_image = bool(data[2]['image'])
                bonding_dict[atom_1]['bonds'].append(
                    {'species': site_2.species, 'index': atom_2, 'length': bond_length,
                     'is_image': is_image, 'position': site_2.coords}
                )
                bonding_dict[atom_2]['bonds'].append(
                    {'species': site_1.species, 'index': atom_1, 'length': bond_length,
                     'is_image': is_image, 'position': site_1.coords}
                )

            for key in bonding_dict:
                bonding_dict[key]['bonds'] = sorted(bonding_dict[key]['bonds'], key=lambda x: x['index'])

            self._bonding_stats = {key: bonding_dict[key] for key in sorted(bonding_dict)}

        return self._bonding_stats

    def draw_network(self, layout=None):
        """ Draw the CrystalGraph network. """
        from matador.crystal.network import draw_network
        draw_network(self, layout=layout)
