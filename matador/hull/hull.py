# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements convex hull functionality from database
queries.

"""


from copy import deepcopy
from typing import List
from collections import defaultdict
import re
import warnings

from bson.son import SON
import pymongo as pm
import numpy as np

from matador.utils.print_utils import print_notify, print_warning
from matador.utils.chem_utils import parse_element_string, get_padded_composition, get_num_intercalated
from matador.utils.chem_utils import get_generic_grav_capacity, get_formula_from_stoich, get_stoich_from_formula
from matador.utils.chem_utils import get_formation_energy
from matador.utils.cursor_utils import set_cursor_from_array, get_array_from_cursor
from matador.utils.cursor_utils import display_results, recursive_get
from matador.utils.cursor_utils import filter_cursor_by_chempots
from matador.battery import Electrode, VoltageProfile
from matador.hull.phase_diagram import PhaseDiagram

# general small number used when comparing energies to zero
EPS = 1e-8
# small number used for checking boundaries
BOUNDARY_EPS = 1e-12


class QueryConvexHull:
    """ Construct a binary or ternary phase diagram from a
    matador.query.DBQuery object, or a list of structures.

    Attributes:
        cursor (list): list of all structures used to create phase diagram.
        hull_cursor (list): list of all documents within hull_cutoff.
        chempot_cursor (list): list of chemical potential documents.
        structures (numpy.ndarray): all structures used to create hull.
        hull_dist (np.ndarray): array of distances from hull for each structure.
        species (list): list of chemical potential symbols.
        num_elements (int): number of elements present in the chemical potentials.
        elements (list): the elements present in the convex hull.
        voltage_data (list): if voltage_curve() has been called, then
            this is a lsit of VoltageProfile objects containing Q, x, V and reaction pathways.
        volume_data (dict): if volume_curve() has been called, then
            this is a dictionary containing x and volumes per B (in AxBy).

    """

    def __init__(self, query=None, cursor=None, elements=None, species=None, voltage=False, volume=False, subcmd=None,
                 plot_kwargs=None, lazy=False, energy_key='enthalpy_per_atom', client=None, collections=None, db=None,
                 **kwargs):
        """ Initialise the class from either a DBQuery or a cursor (list
        of matador dicts) and construct the appropriate phase diagram.

        Keyword arguments:
            query (matador.query.DBQuery): object containing structures,
            cursor (list(dict)): alternatively specify list of matador documents.
            species (list(str)): list of elements/chempots to use, used to provide a useful order,
            voltage (bool): whether or nto to compute voltages relative for insertion of first entry in species,
            volume (bool): whether or not to compute volume expansion relative to first entry in species,
            energy_key (str): key under which the desired energy *per atom* is stored.
            lazy (bool): if True, do not create hull until `self.create_hull()` is called
            chempots (list(float)): list of chemical potential values to use.
            elements (list(str)): deprecated form `species`.
            kwargs (dict): mostly CLI arguments, see matador hull --help for full options.
            plot_kwargs (dict): arguments to pass to plot_hull function
            client (pymongo.MongoClient): optional client to pass to DBQuery.
            collections (dict of pymongo.collections.Collection): optional dict of collections to pass to DBQuery.
            db (str): db name to connect to in DBQuery.

        """
        self.args = dict()
        if query is not None:
            self.args.update(query.args)
        self.args.update(kwargs)

        if subcmd is not None:
            warnings.warn("subcmd will soon be deprecated, please pass the equivalent flag as a kwarg (e.g. voltage=True)")
            if subcmd == 'voltage':
                voltage = True
        self.args['subcmd'] = subcmd

        if plot_kwargs is None:
            plot_kwargs = {}
        self.args['plot_kwargs'] = plot_kwargs

        self.from_cursor = False
        self.plot_params = False

        self.compute_voltages = voltage
        self.compute_volumes = volume

        if query is None and cursor is None:
            # if no query or cursor passed, push all kwargs to new query
            from matador.query import DBQuery
            kwargs.pop('intersection', None)
            query = DBQuery(
                subcmd='hull',
                intersection=True,
                client=client,
                collections=collections,
                **kwargs
            )

        self._query = query

        if self._query is not None:
            # this isn't strictly necessary but it maintains the sanctity of the query results
            self.cursor = list(deepcopy(query.cursor))
            self.args['use_source'] = False
            if self._query.args['subcmd'] not in ['hull', 'voltage', 'hulldiff']:
                print_warning('Query was not prepared with subcmd=hull, so cannot guarantee consistent formation energies.')
        else:
            self.cursor = list(cursor)
            self.from_cursor = True
            self.args['use_source'] = True

        # set up attributes for later
        self.structures = None
        self.convex_hull = None
        self.chempot_cursor = None
        self.hull_cursor = None
        self.phase_diagram = None
        self.hull_dist = None
        self.species = None
        self.voltage_data: List[VoltageProfile] = []
        self.volume_data = defaultdict(list)
        self.elements = []
        self.num_elements = 0
        self.species = self._get_species(species, elements)
        self._dimension = len(self.species)

        # tracker for whether per_b fields have been setup
        self._per_b_done = False

        self.energy_key = energy_key
        if not self.energy_key.endswith('_per_atom'):
            warnings.warn('Appending per_atom to energy_key {}'.format(self.energy_key))
            self.energy_key += '_per_atom'

        self._extensive_energy_key = self.energy_key.split('_per_atom')[0]

        if self.args.get('hull_cutoff') is not None:
            self.hull_cutoff = float(self.args['hull_cutoff'])
        else:
            self.hull_cutoff = 0.0

        if self.cursor is None:
            raise RuntimeError('Failed to find structures to create hull!')

        self._non_elemental = False
        assert isinstance(self.species, list)
        for _species in self.species:
            if len(parse_element_string(_species, stoich=True)) > 1:
                self._non_elemental = True

        if not lazy:
            self.create_hull()

    def create_hull(self):
        """ Begin the hull creation routines and perform the
        post-processing specified by initial arguments.

        """
        if self.args.get('uniq'):
            from matador.utils.cursor_utils import filter_unique_structures
            if self.args.get('uniq') is True:
                sim_tol = 0.1
            else:
                sim_tol = self.args.get('uniq')
            print_notify('Filtering for unique structures...')
            self.cursor = filter_unique_structures(
                self.cursor,
                args=self.args,
                quiet=True,
                sim_tol=sim_tol,
                hull=True,
                energy_key=self.energy_key
            )

        self.construct_phase_diagram()

        if not self.hull_cursor:
            print_warning('No structures on hull with chosen chemical potentials.')
        else:
            print_notify('{} structures found within {} eV of the hull, including chemical potentials.'
                         .format(len(self.hull_cursor), self.hull_cutoff))

        display_results(self.hull_cursor, hull=True, energy_key=self.energy_key, **self.args)

        if self.compute_voltages:
            print("Constructing electrode system with active ion: {}".format(self.species[0]))
            self.voltage_curve([doc for doc in self.hull_cursor if doc['hull_distance'] <= 1e-9])

        if self.compute_volumes:
            self.volume_curve()

        if not self.args.get('no_plot'):
            if self.compute_voltages and self.voltage_data:
                print('plotting voltage')
                self.plot_voltage_curve(show=False)
            if self.compute_volumes and self.volume_data:
                self.plot_volume_curve(show=False)

            self.plot_hull(**self.args['plot_kwargs'], debug=self.args.get('debug'), show=True)

    def __repr__(self):
        return display_results(self.hull_cursor, args=self.args, hull=True, energy_key=self.energy_key, return_str=True)

    @property
    def savefig(self):
        """ True if any figure type argument was passed. """
        return any([self.args.get('pdf'), self.args.get('png'), self.args.get('svg')])

    @property
    def savefig_ext(self):
        """ Returns the figure extension to save with. """
        for ext in ['pdf', 'png', 'svg']:
            if self.args.get(ext):
                return ext

    def plot_hull(self, **kwargs):
        """ Hull plot helper function. """
        from matador import plotting
        if self._dimension == 3:
            ax = plotting.plot_ternary_hull(self, **kwargs)
        elif self._dimension == 2:
            ax = plotting.plot_2d_hull(self, **kwargs)
        else:
            print_notify('Unable to plot phase diagram of dimension {}.'.format(self._dimension))
        return ax

    def plot_voltage_curve(self, **kwargs):
        """ Voltage plot helper function. """
        from matador import plotting
        savefig = None
        if self.savefig:
            savefig = '{}_voltage.{}'.format(''.join(self.elements), self.savefig_ext)

        plotting.plot_voltage_curve(
            self.voltage_data,
            expt=self.args.get('expt'),
            savefig=savefig,
            labels=self.args.get('labels'),
            expt_label=self.args.get('expt_label'),
            **kwargs
        )

    def plot_volume_curve(self, **kwargs):
        """ Volume plot helper function. """
        from matador import plotting
        plotting.plot_volume_curve(
            self,
            **self.args['plot_kwargs'],
            show=False
        )

    def fake_chempots(self, custom_elem=None):
        """ Spoof documents for command-line chemical potentials.

        Keyword arguments:
            custom_elem (list(str)): list of element symbols to generate chempots for.

        """
        self.chempot_cursor = []
        print('Generating fake chempots...')

        if custom_elem is None:
            custom_elem = self.species

        if len(custom_elem) != len(self.species):
            raise RuntimeError('Wrong number of compounds/chemical potentials specified: {} vs {}'
                               .format(custom_elem, self.args.get('chempots')))
        for i, _ in enumerate(self.args.get('chempots')):
            self.chempot_cursor.append(dict())
            self.chempot_cursor[i]['stoichiometry'] = get_stoich_from_formula(custom_elem[i])
            self.chempot_cursor[i][self.energy_key] = -1*abs(self.args.get('chempots')[i])
            self.chempot_cursor[i][self._extensive_energy_key] = self.chempot_cursor[i][self.energy_key] * \
                sum(elem[1] for elem in self.chempot_cursor[i]['stoichiometry'])
            self.chempot_cursor[i]['num_fu'] = 1
            self.chempot_cursor[i]['num_atoms'] = 1
            self.chempot_cursor[i]['text_id'] = ['command', 'line']
            self.chempot_cursor[i]['_id'] = None
            self.chempot_cursor[i]['source'] = ['command_line']
            self.chempot_cursor[i]['space_group'] = 'xxx'
            self.chempot_cursor[i][self._extensive_energy_key + '_per_b'] = self.chempot_cursor[i][self.energy_key]
            self.chempot_cursor[i]['num_a'] = 0
            self.chempot_cursor[i]['cell_volume'] = 1
            self.chempot_cursor[i]['concentration'] = [1 if i == ind else 0 for ind in range(self._dimension-1)]
        self.chempot_cursor[0]['num_a'] = float('inf')
        notify = 'Custom chempots:'
        for chempot in self.chempot_cursor:
            notify += '{:3} = {} eV/fu, '.format(get_formula_from_stoich(chempot['stoichiometry'], sort=False),
                                                 chempot[self._extensive_energy_key])

        if self.args.get('debug'):
            for match in self.chempot_cursor:
                print(match)
        print(len(notify) * '─')
        print(notify)
        print(len(notify) * '─')

    def construct_phase_diagram(self):
        """ Create a phase diagram with arbitrary chemical potentials.

        Expects self.cursor to be populated with structures and chemical potential
        labels to be set under self.species.

        """
        self.set_chempots()
        self.cursor = filter_cursor_by_chempots(self.species, self.cursor)

        formation_key = 'formation_{}'.format(self.energy_key)
        extensive_formation_key = 'formation_{}'.format(self._extensive_energy_key)
        for ind, doc in enumerate(self.cursor):
            self.cursor[ind][formation_key] = get_formation_energy(self.chempot_cursor, doc,
                                                                   energy_key=self.energy_key)
            self.cursor[ind][extensive_formation_key] = doc[formation_key] * doc['num_atoms']

        if self._non_elemental and self.args.get('subcmd') in ['voltage', 'volume']:
            raise NotImplementedError('Pseudo-binary/pseudo-ternary voltages not yet implemented.')

        self.phase_diagram = PhaseDiagram(self.cursor, formation_key, self._dimension)
        # aliases for data stored in phase diagram
        self.structures = self.phase_diagram.structures
        self.hull_dist = self.phase_diagram.hull_dist
        self.convex_hull = self.phase_diagram.convex_hull

        # ensure hull cursor is sorted by enthalpy_per_atom,
        # then by concentration, as it will be by default if from database
        hull_cursor = [self.cursor[idx] for idx in np.where(self.hull_dist <= self.hull_cutoff + EPS)[0]]
        # TODO: check why this fails when the opposite way around
        hull_cursor = sorted(hull_cursor, key=lambda doc: (doc['concentration'], recursive_get(doc, self.energy_key)))

        # by default hull cursor includes all structures within hull_cutoff
        # if summary requested and we're in hulldiff mode, filter hull_cursor for lowest per stoich
        if self.args.get('summary') and self.args['subcmd'] == 'hulldiff':
            tmp_hull_cursor = []
            compositions = set()
            for ind, member in enumerate(hull_cursor):
                formula = get_formula_from_stoich(member['stoichiometry'], sort=True)
                if formula not in compositions:
                    compositions.add(formula)
                    tmp_hull_cursor.append(member)

            hull_cursor = tmp_hull_cursor

        self.hull_cursor = hull_cursor

    def set_chempots(self, energy_key=None):
        """ Search for chemical potentials that match the structures in
        the query cursor and add them to the cursor. Also set the concentration
        of chemical potentials in :attr:`cursor`, if not already set.

        """
        if energy_key is None:
            energy_key = self.energy_key
        query = self._query
        query_dict = dict()
        species_stoich = [sorted(get_stoich_from_formula(spec, sort=False)) for spec in self.species]
        self.chempot_cursor = []

        if self.args.get('chempots') is not None:
            self.fake_chempots(custom_elem=self.species)

        elif self.from_cursor:
            chempot_cursor = sorted([doc for doc in self.cursor if doc['stoichiometry'] in species_stoich],
                                    key=lambda doc: recursive_get(doc, energy_key))

            for species in species_stoich:
                for doc in chempot_cursor:
                    if doc['stoichiometry'] == species:
                        self.chempot_cursor.append(doc)
                        break

            if len(self.chempot_cursor) != len(self.species):
                raise RuntimeError('Found {} of {} required chemical potentials'
                                   .format(len(self.chempot_cursor), len(self.species)))

            if self.args.get('debug'):
                print([mu['stoichiometry'] for mu in self.chempot_cursor])

        else:
            print(60 * '─')
            self.chempot_cursor = len(self.species) * [None]
            # scan for suitable chem pots in database
            for ind, elem in enumerate(self.species):

                print('Scanning for suitable', elem, 'chemical potential...')
                query_dict['$and'] = deepcopy(list(query.calc_dict['$and']))

                if not self.args.get('ignore_warnings'):
                    query_dict['$and'].append(query.query_quality())

                if len(species_stoich[ind]) == 1:
                    query_dict['$and'].append(query.query_composition(custom_elem=[elem]))
                else:
                    query_dict['$and'].append(query.query_stoichiometry(custom_stoich=[elem]))

                # if oqmd, only query composition, not parameters
                if query.args.get('tags') is not None:
                    query_dict['$and'].append(query.query_tags())

                mu_cursor = query.repo.find(SON(query_dict)).sort(energy_key, pm.ASCENDING)
                if mu_cursor.count() == 0:
                    print_notify('Failed... searching without spin polarization field...')
                    scanned = False
                    while not scanned:
                        for idx, dicts in enumerate(query_dict['$and']):
                            for key in dicts:
                                if key == 'spin_polarized':
                                    del query_dict['$and'][idx][key]
                                    break
                            if idx == len(query_dict['$and']) - 1:
                                scanned = True
                    mu_cursor = query.repo.find(SON(query_dict)).sort(energy_key, pm.ASCENDING)

                if mu_cursor.count() == 0:
                    raise RuntimeError('No chemical potentials found for {}...'.format(elem))

                self.chempot_cursor[ind] = mu_cursor[0]
                if self.chempot_cursor[ind] is not None:
                    print('Using', ''.join([self.chempot_cursor[ind]['text_id'][0], ' ',
                                            self.chempot_cursor[ind]['text_id'][1]]), 'as chem pot for', elem)
                    print(60 * '─')
                else:
                    raise RuntimeError('No possible chem pots available for {}.'.format(elem))

            for i, mu in enumerate(self.chempot_cursor):
                self.chempot_cursor[i][self._extensive_energy_key + '_per_b'] = mu[energy_key]
                self.chempot_cursor[i]['num_a'] = 0

            self.chempot_cursor[0]['num_a'] = float('inf')

        # don't check for IDs if we're loading from cursor
        if not self.from_cursor:
            ids = [doc['_id'] for doc in self.cursor]
            if self.chempot_cursor[0]['_id'] is None or self.chempot_cursor[0]['_id'] not in ids:
                self.cursor.insert(0, self.chempot_cursor[0])
            for match in self.chempot_cursor[1:]:
                if match['_id'] is None or match['_id'] not in ids:
                    self.cursor.append(match)

        # add faked chempots to overall cursor
        elif self.args.get('chempots') is not None:
            self.cursor.insert(0, self.chempot_cursor[0])
            self.cursor.extend(self.chempot_cursor[1:])

        # find all elements present in the chemical potentials
        elements = []
        for mu in self.chempot_cursor:
            for elem, _ in mu['stoichiometry']:
                if elem not in elements:
                    elements.append(elem)
        self.elements = elements
        self.num_elements = len(elements)

    @staticmethod
    def filter_cursor_by_chempots(species, cursor):
        """ For the desired chemical potentials, remove any incompatible structures
        from cursor.

        Parameters:
            species (list): list of chemical potential formulae.
            cursor (list): list of matador documents to filter.

        Returns:
            list: the filtered cursor.
        """

        return filter_cursor_by_chempots(species, cursor)

    def _setup_per_b_fields(self):
        """ Calculate the enthalpy and volume per "B" in A_x B. """
        if self._per_b_done:
            return
        for ind, doc in enumerate(self.cursor):
            nums_b = len(self.species[1:]) * [0]
            for elem in doc['stoichiometry']:
                for chem_pot_ind, chem_pot in enumerate(self.species[1:]):
                    if elem[0] == chem_pot:
                        nums_b[chem_pot_ind] += elem[1]

            num_b = sum(nums_b)
            num_fu = doc['num_fu']
            if num_b == 0:
                self.cursor[ind][self._extensive_energy_key + '_per_b'] = 12345e5
                self.cursor[ind]['cell_volume_per_b'] = 12345e5
            else:
                self.cursor[ind][self._extensive_energy_key + '_per_b'] = doc[self._extensive_energy_key] / (num_b * num_fu)
                if 'cell_volume' in doc:
                    self.cursor[ind]['cell_volume_per_b'] = doc['cell_volume'] / (num_b * num_fu)

        capacities = np.zeros((len(self.cursor)))
        for i, doc in enumerate(self.cursor):
            concs = self.cursor[i]['concentration']
            concs = get_padded_composition(doc['stoichiometry'], self.elements)
            capacities[i] = get_generic_grav_capacity(concs, self.elements)
        set_cursor_from_array(self.cursor, capacities, 'gravimetric_capacity')
        self._per_b_done = True

    def voltage_curve(self, hull_cursor):
        """ Take a computed convex hull and calculate voltages for either binary or ternary
        systems. Sets the self.voltage_data attribute with various fields.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        """
        if not self._non_elemental:
            self._setup_per_b_fields()
        if self._dimension == 2:
            self._calculate_binary_voltage_curve(hull_cursor)
        elif self._dimension == 3:
            self._calculate_ternary_voltage_curve(hull_cursor)
        else:
            raise RuntimeError('Unable to calculate voltage curve for hull of dimension {}'.format(self._dimension))

        self._scale_voltages(self.species[0])

        for profile in self.voltage_data:
            print(profile.voltage_summary(csv=self.args.get('csv', False)))

    def volume_curve(self):
        """ Take stable compositions and volume and calculate
        volume expansion per "B" in AB binary.

        """

        if not self._non_elemental and not self._per_b_done:
            self._setup_per_b_fields()
        if self._dimension == 2:
            self._calculate_binary_volume_curve()
        elif self._dimension == 3:
            # dimension 3 is implemented inside the voltage curve call.
            pass
        else:
            raise NotImplementedError(
                "Volume curves have only been implemented for binary and ternary phase diagrams."
            )

        self.print_volume_summary()

    def print_volume_summary(self):
        """ Prints a volume data summary.

        If self.args['csv'] is True, save the summary to a file.

        """
        for reaction_idx, _ in enumerate(self.volume_data['Q']):
            data_str = '# Reaction {} \n'.format(reaction_idx+1)
            data_str += '# ' + ''.join(get_formula_from_stoich(self.volume_data['endstoichs'][reaction_idx])) + '\n'
            data_str += '# {:>10}\t{:>14} \t{:>14}\n'.format('Q (mAh/g)', 'Volume (A^3)', 'Volume ratio with bulk')
            for idx, _ in enumerate(self.volume_data['electrode_volume'][reaction_idx]):
                data_str += ('{:>10.2f} \t{:14.2f} \t{:14.2f}'
                             .format(self.volume_data['Q'][reaction_idx][idx],
                                     self.volume_data['electrode_volume'][reaction_idx][idx],
                                     self.volume_data['volume_ratio_with_bulk'][reaction_idx][idx]))
                if idx != len(self.volume_data['Q'][reaction_idx]) - 1:
                    data_str += '\n'
            if self.args.get('csv'):
                with open(''.join(self.species) + '_volume.csv', 'w') as f:
                    f.write(data_str)
            print('\nVolume data:')
            print('\n' + data_str)

    def get_hull_distances(self, *args, **kwargs):
        """ Wrapper to PhaseDiagram.get_hull_distances. """
        return self.phase_diagram.get_hull_distances(*args, **kwargs)

    def _calculate_binary_voltage_curve(self, hull_cursor):
        """ Generate binary voltage curve, setting the self.voltage_data
        dictionary.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        """
        mu_enthalpy = get_array_from_cursor(self.chempot_cursor, self.energy_key)
        # take a copy of the cursor so it can be reordered
        _hull_cursor = deepcopy(hull_cursor)
        x = get_num_intercalated(_hull_cursor)
        _x_order = np.argsort(x)
        capacities = np.asarray([
            get_generic_grav_capacity(
                get_padded_composition(doc['stoichiometry'], self.elements),
                self.elements)
            for doc in _hull_cursor
        ])
        set_cursor_from_array(_hull_cursor, x, 'conc_of_active_ion')
        set_cursor_from_array(_hull_cursor, capacities, 'gravimetric_capacity')
        _hull_cursor = sorted(_hull_cursor, key=lambda doc: doc['conc_of_active_ion'])
        capacities = capacities[_x_order]
        x = np.sort(x)
        stable_enthalpy_per_b = get_array_from_cursor(_hull_cursor,
                                                      self._extensive_energy_key + '_per_b')

        x, uniq_idxs = np.unique(x, return_index=True)
        stable_enthalpy_per_b = stable_enthalpy_per_b[uniq_idxs]
        capacities = capacities[uniq_idxs]

        V = np.zeros_like(x)
        V[1:] = - np.diff(stable_enthalpy_per_b) / np.diff(x) + mu_enthalpy[0]
        V[0] = V[1]

        reactions = [((None, get_formula_from_stoich(doc['stoichiometry'])),) for doc in _hull_cursor[1:]]

        # make V, Q and x available for plotting, stripping NaNs to re-add later
        # in the edge case of duplicate chemical potentials
        capacities, unique_caps = np.unique(capacities, return_index=True)
        non_nan = np.argwhere(np.isfinite(capacities))
        V = (np.asarray(V)[unique_caps])[non_nan].flatten().tolist()
        x = (np.asarray(x)[unique_caps])[non_nan].flatten().tolist()
        capacities = capacities[non_nan].flatten().tolist()

        x.append(np.nan)
        capacities.append(np.nan)
        V.append(0.0)

        average_voltage = Electrode.calculate_average_voltage(capacities, V)
        profile = VoltageProfile(
            voltages=V,
            capacities=capacities,
            average_voltage=average_voltage,
            starting_stoichiometry=[[self.species[1], 1.0]],
            reactions=reactions,
            active_ion=self.species[0]
        )
        self.voltage_data.append(profile)

    def _calculate_ternary_voltage_curve(self, hull_cursor):
        """ Calculate tenary voltage curve, setting self.voltage_data.
        First pass written by James Darby, jpd47@cam.ac.uk.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        """

        # construct working array of concentrations and energies
        stoichs = get_array_from_cursor(hull_cursor, 'stoichiometry')

        # do another convex hull on just the known hull points, to allow access to useful indices
        import scipy.spatial

        # construct working array of concentrations and energies
        points = np.hstack((
            get_array_from_cursor(hull_cursor, 'concentration'),
            get_array_from_cursor(hull_cursor, self.energy_key).reshape(len(hull_cursor), 1)
        ))
        stoichs = get_array_from_cursor(hull_cursor, 'stoichiometry')

        # do another convex hull on just the known hull points, to allow access to useful indices
        convex_hull = scipy.spatial.ConvexHull(points)

        endpoints, endstoichs = Electrode._find_starting_materials(convex_hull.points, stoichs)
        print('{} starting point(s) found.'.format(len(endstoichs)))
        for endstoich in endstoichs:
            print(get_formula_from_stoich(endstoich), end=' ')
        print('\n')

        # iterate over possible delithiated phases
        for reaction_ind, endpoint in enumerate(endpoints):
            print(30 * '-')
            print('Reaction {}, {}:'.format(reaction_ind+1, get_formula_from_stoich(endstoichs[reaction_ind])))
            reactions, capacities, voltages, average_voltage, volumes = self._construct_electrode(
                convex_hull, endpoint, endstoichs[reaction_ind], hull_cursor
            )

            profile = VoltageProfile(
                starting_stoichiometry=endstoichs[reaction_ind],
                reactions=reactions,
                capacities=capacities,
                voltages=voltages,
                average_voltage=average_voltage,
                active_ion=self.species[0]
            )
            self.voltage_data.append(profile)

            if not self._non_elemental:
                self.volume_data['x'].append(capacities)
                self.volume_data['Q'].append(capacities)
                self.volume_data['electrode_volume'].append(volumes)
                self.volume_data['endstoichs'].append(endstoichs[reaction_ind])
                self.volume_data['volume_ratio_with_bulk'].append(np.asarray(volumes) / volumes[0])
                self.volume_data['volume_expansion_percentage'].append(((np.asarray(volumes) / volumes[0]) - 1) * 100)
                self.volume_data['hull_distances'].append(np.zeros_like(capacities))

    def _calculate_binary_volume_curve(self):
        """ Take stable compositions and volume and calculate volume
        expansion per "B" in AB binary.

        """
        stable_comp = get_array_from_cursor(self.hull_cursor, 'concentration')
        stable_comp, unique_comp_inds = np.unique(stable_comp, return_index=True)
        for doc in self.hull_cursor:
            if 'cell_volume_per_b' not in doc:
                raise RuntimeError("Document missing key `cell_volume_per_b`: {}".format(doc))
        stable_stoichs = get_array_from_cursor(self.hull_cursor, 'stoichiometry')[unique_comp_inds]
        stable_vol = get_array_from_cursor(self.hull_cursor, 'cell_volume_per_b')[unique_comp_inds]
        stable_cap = get_array_from_cursor(self.hull_cursor, 'gravimetric_capacity')[unique_comp_inds]
        hull_distances = get_array_from_cursor(self.hull_cursor, 'hull_distance')[unique_comp_inds]
        with np.errstate(divide='ignore'):
            stable_x = stable_comp / (1 - stable_comp)
        non_nans = np.argwhere(np.isfinite(stable_x))
        self.volume_data['x'].append(stable_x[non_nans].flatten())
        self.volume_data['Q'].append(stable_cap[non_nans].flatten())
        self.volume_data['electrode_volume'].append(stable_vol[non_nans].flatten())
        self.volume_data['volume_ratio_with_bulk'].append((stable_vol[non_nans] / stable_vol[0]).flatten())
        self.volume_data['volume_expansion_percentage'].append(((stable_vol[non_nans] / stable_vol[0]) - 1) * 100)
        self.volume_data['hull_distances'].append(hull_distances[non_nans].flatten())
        self.volume_data['endstoichs'].append(stable_stoichs[non_nans].flatten()[0])

    def _get_species(self, species=None, elements=None):
        """ Try to determine species for phase diagram from arguments or
        passed data. Sets the `self.species` attribute.

        Keyword arguments:
            species (str/list): the species kwarg passed to __init__.
            elements (str/list): the elements kwarg passed to __init__.

        """
        # try to determine species for phase diagram, from arguments and data
        species = species or elements
        if species is None:
            # handles element list passed by query: should only ever be a list of one entry
            if self.args.get('composition') is not None:
                if isinstance(self.args.get('composition'), list):
                    species = self.args.get('composition')[0]
                else:
                    species = self.args.get('composition')
                if ':' in species:
                    species = species.split(':')
        if species is None:
            if self.from_cursor:
                species = list({atom for doc in self.cursor for atom in doc['atom_types']})
            else:
                raise RuntimeError("Unable to determine species in hull, please specify with "
                                   "e.g. `elements=['Li', 'P']")

        # handles when species is e.g. ['LiCo2:Sn2S']
        if isinstance(species, str):
            if ':' in species:
                species = species.split(':')
            else:
                species = [spec for spec in re.split(r'([A-Z][a-z]*)', species) if spec]

        # edge case where user passes e.g. ['KP'], when they mean ['K', 'P'].
        if isinstance(species, list) and len(species) == 1:
            tmp_species = []
            for spec in species:
                if ':' in spec:
                    tmp_species.append(spec.split(':'))
                else:
                    tmp_species.extend([item for item in re.split(r'([A-Z][a-z]*)', spec) if item])
            species = tmp_species

        return species

    def _construct_electrode(self, hull, endpoint, endstoich, hull_cursor):

        intersections, crossover = self._find_hull_pathway_intersections(
            endpoint, endstoich, hull
        )
        return self._compute_voltages_from_intersections(
            intersections, hull, crossover, endstoich, hull_cursor
        )

    def _compute_voltages_from_intersections(
            self, intersections, hull, crossover, endstoich, hull_cursor
    ):
        """ Traverse the composition pathway and its intersections with hull facets
        and compute the voltage and volume drops along the way, for a ternary system
        with N 3-phase regions.

        Parameters:
            intersections (np.ndarray): Nx3 array containing the list of face indices
                crossed by the path.
            hull (scipy.spatial.ConvexHull): the temporary hull object that is referred
                to by any indices.
            crossover (np.ndarray): (N+1)x3 array containing the coordinates at which
                the composition pathway crosses the hull faces.
            endstoich (list): a matador stoichiometry that defines the end of the pathway.
            hull_cursor (list): the actual structures used to create the hull.

        """
        if len(crossover) != len(intersections) + 1:
            raise RuntimeError("Incompatible number of crossovers ({}) and intersections ({})."
                               .format(np.shape(crossover), np.shape(intersections)))

        # set up output arrays
        voltages = np.empty(len(intersections) + 1)
        volumes = np.empty(len(intersections))
        reactions = []

        # set up some useful arrays for later
        points = hull.points
        mu_enthalpy = get_array_from_cursor(self.chempot_cursor, self.energy_key)
        if not self._non_elemental:
            input_volumes = get_array_from_cursor(hull_cursor, 'cell_volume') / get_array_from_cursor(hull_cursor, 'num_atoms')
        capacities = [get_generic_grav_capacity(point, self.species) for point in crossover]

        # initial composition of one formula unit of the starting material
        initial_comp = get_padded_composition(endstoich, self.species)

        # loop over intersected faces and compute the voltage and volume at the
        # crossover with that face
        for ind, face in enumerate(intersections):

            simplex_index = int(face[0])

            final_stoichs = [hull_cursor[idx]['stoichiometry'] for idx in hull.simplices[simplex_index]]
            atoms_per_fu = [sum([elem[1] for elem in stoich]) for stoich in final_stoichs]

            energy_vec = points[hull.simplices[simplex_index], 2]
            comp = points[hull.simplices[simplex_index], :]
            comp[:, 2] = 1 - comp[:, 0] - comp[:, 1]

            # normalize the crossover composition to one formula unit of the starting electrode
            norm = np.asarray(crossover[ind][1:]) / np.asarray(initial_comp[1:])
            ratios_of_phases = np.linalg.solve(comp.T, crossover[ind] / norm[0])

            # remove small numerical noise
            ratios_of_phases[np.where(ratios_of_phases < EPS)] = 0

            # create a list containing the sections of the balanced reaction for printing
            balanced_reaction = []
            for i, ratio in enumerate(ratios_of_phases):
                if ratio < EPS:
                    continue
                else:
                    formula = get_formula_from_stoich(final_stoichs[i])
                    balanced_reaction.append((ratio / atoms_per_fu[i], formula))
            reactions.append(balanced_reaction)

            # compute the voltage as the difference between the projection of the gradient down to pure active ion,
            # and the reference chemical potential
            comp_inv = np.linalg.inv(comp.T)
            V = -(comp_inv.dot([1, 0, 0])).dot(energy_vec)
            V = V + mu_enthalpy[0]

            # compute the volume of the final electrode
            if not self._non_elemental:
                if ind <= len(intersections) - 1:
                    volume_vec = input_volumes[hull.simplices[simplex_index]]
                    volumes[ind] = np.dot(ratios_of_phases, volume_vec)

            # double up on first voltage
            if ind == 0:
                voltages[0] = V
            voltages[ind+1] = V

        average_voltage = Electrode.calculate_average_voltage(capacities, voltages)

        # print the reaction over a few lines, remove 1s and rounded ratios
        print(
            ' ---> '.join(
                [' + '.join(
                    ["{}{}".format(
                        str(round(chem[0], 3)) + ' ' if abs(chem[0] - 1) > EPS else '',
                        chem[1])
                     for chem in region])
                 for region in reactions])
        )

        return reactions, capacities, voltages, average_voltage, volumes

    def _find_hull_pathway_intersections(self, endpoint, endstoich, hull):
        """ This function traverses a ternary phase diagram on a path from `endpoint`
        towards [1, 0, 0], i.e. a pure phase of the active ion. The aim is to find
        all of the faces and intersections along this pathway, making sure to only
        add unique intersections, and making appropriate (symmetry-breaking) choices
        of which face to use. We proceed as follows:

            1. Filter out any vertices that contain only binary phases, or that contain
                only the chemical potentials.
            2. Loop over all faces of the convex hull and test if the pathway touches/intersects
                that face. If it does, ascertain where the intersection occurs, and whether
                it is at one, two or none of the vertices of the face. Create an array of unique
                intersections for this face and save them for later.
            3. Loop over the zeros found between for each face. If the pathway only grazes a single
                point on this face, ignore it. If one of the edges is parallel to the pathway,
                make sure to include this face and the faces either side of it.

        Parameters:
            endpoint (np.ndarray): array containing composition (2D) and energy of
                start point.

        """
        import scipy.spatial.distance
        import itertools

        endpoint = endpoint[:-1]
        compositions = np.zeros_like(hull.points)
        compositions[:, 0] = hull.points[:, 0]
        compositions[:, 1] = hull.points[:, 1]
        compositions[:, 2] = 1 - hull.points[:, 0] - hull.points[:, 1]

        # define the composition pathway through ternary space
        gradient = endpoint[1] / (endpoint[0] - 1)
        # has to intersect [1, 0] so c = y0 - m*x0 = -m
        y0 = -gradient

        # filter the vertices we're going to consider
        skip_inds = set()
        for simp_ind, simplex in enumerate(hull.simplices):
            vertices = np.asarray([compositions[i] for i in simplex])
            # skip all simplices that contain only binary phases
            for i in range(3):
                if np.max(vertices[:, i]) < BOUNDARY_EPS:
                    skip_inds.add(simp_ind)
                    continue
            # skip the outer triangle formed by the chemical potentials
            if all(np.max(vertices, axis=-1) > 1 - BOUNDARY_EPS):
                skip_inds.add(simp_ind)
                continue

        two_phase_crossover = []
        three_phase_crossover = []

        for simp_ind, simplex in enumerate(hull.simplices):
            if simp_ind in skip_inds:
                continue
            vertices = np.asarray([compositions[i] for i in simplex])

            # put each vertex of triangle into line equation and test their signs
            test = vertices[:, 0] * gradient + y0 - vertices[:, 1]
            test[np.abs(test) < BOUNDARY_EPS] = 0.0
            test = np.sign(test)

            # if there are two different signs in the test array, then the line intersects/grazes this face triangle
            if len(np.unique(test)) > 1:
                # now find intersection points and split them into one, two and three-phase crossover regions
                zeros = [val for zero in np.where(test == 0) for val in zero.tolist()]
                num_zeros = len(zeros)
                zero_points = []
                # skip_vertex = None
                for zero in zeros:
                    zero_pos = compositions[simplex[zero]]
                    if num_zeros == 2:
                        zero_points.append(zero_pos)
                        two_phase_crossover.append((zero_pos, simp_ind))

                # if we have already found both crossovers, skip to next face
                if num_zeros == 2:
                    continue

                # now find the non-trivial intersections
                for i, j in itertools.combinations(simplex, r=2):
                    # if skip_vertex in [i, j]:
                    #     continue

                    A = compositions[i]
                    B = compositions[j]
                    C = endpoint
                    D = [1, 0]

                    def coeffs(X, Y):
                        # find line equation for edge AB
                        # a x + b y = c
                        return (Y[1] - X[1], X[0] - Y[0])

                    a1, b1 = coeffs(A, B)
                    c1 = a1 * A[0] + b1 * A[1]

                    a2, b2 = coeffs(C, D)
                    c2 = a2 * C[0] + b2 * C[1]
                    det = a1 * b2 - a2 * b1
                    if det == 0:
                        continue

                    x = (b2 * c1 - b1 * c2) / det
                    y = (a1 * c2 - a2 * c1) / det
                    intersection_point = np.asarray([x, y, 1-x-y])

                    # now have to test whether that intersection point is inside the triangle
                    # which we can do by testing that is inside the rectangle spanned by AB
                    x_bound = sorted([A[0], B[0]])
                    y_bound = sorted([A[1], B[1]])
                    if (
                        x_bound[0] - BOUNDARY_EPS/2 <= intersection_point[0] <= x_bound[1] + BOUNDARY_EPS/2 and
                        y_bound[0] - BOUNDARY_EPS/2 <= intersection_point[1] <= y_bound[1] + BOUNDARY_EPS/2
                    ):
                        # three_phase_crossover.append((intersection_point, simp_ind))
                        zero_points.append(intersection_point)

                unique_zeros = []
                for zero in zero_points:
                    if (len(unique_zeros) >= 1 and
                            np.any(scipy.spatial.distance.cdist(unique_zeros, zero.reshape(1, 3)) < BOUNDARY_EPS)):
                        pass
                    else:
                        unique_zeros.append(zero)

                num_zeros = len(unique_zeros)
                if num_zeros == 1:
                    # we never need to use faces which just touch the triangle
                    continue
                elif num_zeros == 2:
                    three_phase_crossover.append((unique_zeros[0], simp_ind))
                    three_phase_crossover.append((unique_zeros[1], simp_ind))

        # loop over points and only add the unique ones to the crossover array
        # if a point exists as a two-phase region, then always add that simplex
        crossover = [[1, 0, 0]]
        simplices = {}
        zeros = sorted(two_phase_crossover + three_phase_crossover, key=lambda x: x[0][0])
        # for multiplicity, zeros in enumerate([one_phase_crossover, two_phase_crossover, three_phase_crossover]):
        for zero, simp_ind in zeros:
            if len(crossover) >= 1 and np.any(scipy.spatial.distance.cdist(crossover, zero.reshape(1, 3)) < BOUNDARY_EPS):
                pass
            else:
                if simp_ind not in simplices:
                    crossover.append(zero)
                    # if multiplicity > 0:
                    simplices[simp_ind] = zero[0]

        intersections = sorted(simplices.items(), key=lambda x: x[1])

        if len(intersections) == 0:
            raise RuntimeError(
                'No intermediate structures found for starting point {}.'
                .format(get_formula_from_stoich(endstoich))
            )

        return intersections, sorted(crossover, key=lambda x: x[0])

    def _scale_voltages(self, species):
        """ Rescale voltages to account for the valence of the active ion,
        as stored in the Electrode class.

        """
        valence_factor = Electrode.valence_data.get(species, None)

        if valence_factor is None:
            warnings.warn("Unable to find valence of species {}, treating it as 1.".format(species))
            return

        if valence_factor != 1:
            for ind, start_point in enumerate(self.voltage_data):
                for j in range(len(self.voltage_data[ind].voltages)):
                    self.voltage_data[ind].voltages[j] *= valence_factor
            self.voltage_data[ind].average_voltage *= valence_factor
