# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements convex hull functionality from database
queries.

"""


from copy import deepcopy
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
from matador.export import generate_hash
from matador.utils.cursor_utils import filter_cursor_by_chempots
from matador.hull.phase_diagram import PhaseDiagram

EPS = 1e-12


class QueryConvexHull:
    """ Construct a binary or ternary phase diagram from a
    matador.query.DBQuery object, or a list of structures.

    Attributes:
        cursor (list): list of all structures used to create phase diagram.
        hull_cursor (list): list of all documents within hull_cutoff.
        chempot_cursor (list): list of chemical potential documents.
        structures (numpy.ndarray): all structures used to create hull.
        structure_slice (numpy.ndarray): array of concentrations and
            formation energies for each structure with E_F <= 0. The indices
            in hull correspond to entries of this array.
        hull_dist (np.ndarray): array of distances from hull for each structure.
        species (list): list of chemical potential symbols.
        num_elements (int): number of elements present in the chemical potentials.
        elements (list): the elements present in the convex hull.
        voltage_data (dict): if voltage_curve() has been called, then
            this is a dictionary containing Q, x, V and reaction pathways.
        volume_data (dict): if volume_curve() has been called, then
            this is a dictionary containing x and volumes per B (in AxBy).

    """

    def __init__(self, query=None, cursor=None, elements=None, species=None, subcmd='hull',
                 plot_kwargs=None, lazy=False, energy_key='enthalpy_per_atom', **kwargs):
        """ Initialise the class from either a DBQuery or a cursor (list
        of matador dicts) and construct the appropriate phase diagram.

        Keyword arguments:
            query (matador.query.DBQuery): object containing structures,
            cursor (list(dict)): alternatively specify list of matador recursive_get(doc, self.energy_key).
            species (list(str)): list of elements/chempots to use, used to provide a useful order,
            elements (list(str)): deprecated form of the above.
            subcmd (str): either 'hull' or 'voltage',
            lazy (bool): if True, do not create hull until self.create_hull() is called
            energy_key (str): key under which the desired energy_per_atom is stored.
            kwargs (dict): mostly CLI arguments, see matador hull --help for full options.
            plot_kwargs (dict): arguments to pass to plot_hull function

        """
        self.args = dict()
        if query is not None:
            self.args.update(query.args)
        self.args.update(kwargs)

        self.devel = True
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd
        if plot_kwargs is None:
            plot_kwargs = {}
        if plot_kwargs.get('show') is None:
            plot_kwargs['show'] = True
        self.args['plot_kwargs'] = plot_kwargs

        self._query = query
        self.from_cursor = False
        self.plot_params = False

        if self._query is not None:
            self.cursor = list(deepcopy(query.cursor))
            self.args['use_source'] = False
            if self._query.args['subcmd'] not in ['hull', 'voltage', 'hulldiff']:
                print_warning('Query was not prepared with subcmd=hull, so cannot guarantee consistent formation energies.')
        else:
            self.cursor = list(deepcopy(cursor))
            self.from_cursor = True
            self.args['use_source'] = True

        # set up attributes for later
        self.chempot_cursor = None
        self.hull_cursor = None
        self.structure_slice = None
        self.phase_diagram = None
        self.hull_dist = None
        self.species = None
        self.voltage_data = {}
        self.volume_data = {}
        self.elements = []
        self.num_elements = 0

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

        self._set_species(species, elements)

        self._non_elemental = False
        assert isinstance(self.species, list)
        for _species in self.species:
            if len(parse_element_string(_species, stoich=True)) > 1:
                self._non_elemental = True

        self._dimension = len(self.species)
        if self._dimension > 2 and self._query is not None and not self._query.args.get('intersection'):
            raise SystemExit('Please query with -int/--intersection when creating ternary+ hulls.')

        if not lazy:
            self.create_hull()

    def create_hull(self):
        """ Begin the hull creation routines and perform the
        post-processing specified by initial arguments.

        """
        self.construct_phase_diagram()

        if not self.hull_cursor:
            print_warning('No structures on hull with chosen chemical potentials.')
        else:
            print_notify('{} structures found within {} eV of the hull, including chemical potentials.'
                         .format(len(self.hull_cursor), self.hull_cutoff))

        display_results(self.hull_cursor, args=self.args, hull=True, energy_key=self.energy_key)

        if not self.args.get('no_plot'):
            from matador import plotting

        if self.args['subcmd'] == 'voltage':
            self.voltage_curve([doc for doc in self.hull_cursor if doc['hull_distance'] <= 1e-9])
            if not self.args.get('no_plot'):
                plotting.plot_voltage_curve(self)
                self.plot_hull(**self.args['plot_kwargs'], debug=self.args.get('debug'))

        if self.args.get('volume'):
            self.volume_curve()
            if not self.args.get('no_plot'):
                plotting.plot_volume_curve(self)

        if self.args['subcmd'] == 'hull' and not self.args.get('no_plot'):
            self.plot_hull(**self.args['plot_kwargs'], debug=self.args.get('debug'))

        if self.args.get('uniq') and self.args['subcmd'] != 'swaps':
            from matador.similarity.similarity import get_uniq_cursor
            if self.args.get('uniq') is True:
                sim_tol = 0.1
            else:
                sim_tol = self.args.get('uniq')
            print_notify('Filtering for unique structures...')
            unique_set, _, _, _ = get_uniq_cursor(self.hull_cursor, debug=self.args.get('debug'), sim_tol=sim_tol)
            old_cursor_len = len(self.hull_cursor)
            self.hull_cursor = [self.hull_cursor[ind] for ind in unique_set]
            display_results(self.hull_cursor, hull=True, args=self.args, energy_key=self.energy_key)
            print('Filtered {} down to {}'.format(old_cursor_len, len(unique_set)))

    @property
    def savefig(self):
        """ True if any figure type argument was passed. """
        return any([self.args.get('pdf'), self.args.get('png'), self.args.get('svg')])

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

    def fake_chempots(self, custom_elem=None):
        """ Spoof documents for command-line chemical potentials.

        Keyword arguments:
            custom_elem (list(str)): list of element symbols to generate chempots for.

        """
        self.chempot_cursor = []

        if custom_elem is None:
            custom_elem = self.species

        if len(custom_elem) != len(self.species):
            raise RuntimeError('Wrong number of compounds/chemical potentials specified: {} vs {}'.format(custom_elem, self.args.get('chempots')))
        for i, _ in enumerate(self.args.get('chempots')):
            self.chempot_cursor.append(dict())
            self.chempot_cursor[i]['stoichiometry'] = get_stoich_from_formula(custom_elem[i])
            self.chempot_cursor[i][self.energy_key] = -1*abs(self.args.get('chempots')[i])
            self.chempot_cursor[i][self._extensive_energy_key] = self.chempot_cursor[i][self.energy_key] * \
                sum(elem[1] for elem in self.chempot_cursor[i]['stoichiometry'])
            self.chempot_cursor[i]['num_fu'] = 1
            self.chempot_cursor[i]['num_atoms'] = 1
            self.chempot_cursor[i]['text_id'] = ['command', 'line']
            self.chempot_cursor[i]['_id'] = generate_hash(hash_len=10)
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
        for ind, doc in enumerate(self.cursor):
            self.cursor[ind][formation_key] = get_formation_energy(self.chempot_cursor, doc,
                                                                   energy_key=self.energy_key)

        if not self._non_elemental:
            self._setup_per_b_fields()
        elif self.args.get('subcmd') in ['voltage', 'volume']:
            raise NotImplementedError('Pseudo-binary/ternary voltages not yet implemented.')

        self.phase_diagram = PhaseDiagram(self.cursor, formation_key, self._dimension)
        # aliases for data stored in phase diagram
        self.structures = self.phase_diagram.structures
        self.structure_slice = self.phase_diagram.structure_slice
        self.hull_dist = self.phase_diagram.hull_dist
        set_cursor_from_array(self.cursor, self.hull_dist, 'hull_distance')
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

    def set_chempots(self):
        """ Search for chemical potentials that match the structures in
        the query cursor and add them to the cursor. Also set the concentration
        of chemical potentials in cursor, if not already set.

        """
        query = self._query
        query_dict = dict()
        species_stoich = [sorted(get_stoich_from_formula(spec, sort=False)) for spec in self.species]
        self.chempot_cursor = []

        if self.args.get('chempots') is not None:
            self.fake_chempots(custom_elem=self.species)

        elif self.from_cursor:
            chempot_cursor = sorted([doc for doc in self.cursor if doc['stoichiometry'] in species_stoich],
                                    key=lambda doc: recursive_get(doc, self.energy_key))

            for species in species_stoich:
                for doc in chempot_cursor:
                    if doc['stoichiometry'] == species:
                        self.chempot_cursor.append(doc)
                        break

            if len(self.chempot_cursor) != len(self.species):
                raise RuntimeError('Found {} of {} required chemical potentials'.format(len(self.chempot_cursor), len(self.species)))

            # for ind, doc in enumerate(self.chempot_cursor):
                # self.chempot_cursor[ind][self._extensive_energy_key + '_per_b'] = doc[self.energy_key]
                # self.chempot_cursor[ind]['cell_volume_per_b'] = float('inf') if ind == 0 else 0.0

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

                mu_cursor = query.repo.find(SON(query_dict)).sort(self.energy_key, pm.ASCENDING)
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
                    mu_cursor = query.repo.find(SON(query_dict)).sort(self.energy_key, pm.ASCENDING)

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
                self.chempot_cursor[i][self._extensive_energy_key + 'per_b'] = mu[self.energy_key]
                self.chempot_cursor[i]['num_a'] = 0

            self.chempot_cursor[0]['num_a'] = float('inf')

        # don't check for IDs if we're loading from cursor
        if not self.from_cursor:
            ids = [doc['_id'] for doc in self.cursor]
            if self.chempot_cursor[0]['_id'] not in ids:
                self.cursor.insert(0, self.chempot_cursor[0])
            for match in self.chempot_cursor[1:]:
                if match['_id'] not in ids:
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
                self.cursor[ind]['cell_volume_per_b'] = doc['cell_volume'] / (num_b * num_fu)

        capacities = np.zeros((len(self.cursor)))
        for i, doc in enumerate(self.cursor):
            concs = self.cursor[i]['concentration']
            concs = get_padded_composition(doc['stoichiometry'], self.elements)
            capacities[i] = get_generic_grav_capacity(concs, self.elements)
        set_cursor_from_array(self.cursor, capacities, 'gravimetric_capacity')

    def voltage_curve(self, hull_cursor):
        """ Take a computed convex hull and calculate voltages for either binary or ternary
        systems. Sets the self.voltage_data attribute with various fields.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        """
        if self._dimension == 2:
            self._calculate_binary_voltage_curve(hull_cursor)
        elif self._dimension == 3:
            self._calculate_ternary_voltage_curve(hull_cursor)
        else:
            raise RuntimeError('Unable to calculate voltage curve for hull of dimension {}'.format(self._dimension))

        data_str = ''
        for ind, path in enumerate(self.voltage_data['Q']):
            if ind != 0:
                data_str += '\n'
            if self._dimension == 3:
                data_str += '# ' + get_formula_from_stoich(self.voltage_data['endstoichs'][ind]) + '\n'
            else:
                data_str += '# ' + ''.join(self.species) + '\n'
            data_str += '# {:^10} \t{:^10} \t{:^10}\n'.format('x', 'Q (mAh/g)', 'Voltage (V)')
            for idx, _ in enumerate(path):
                data_str += '{:>10.2f} \t{:>10.2f} \t{:>10.4f}'.format(self.voltage_data['x'][ind][idx],
                                                                       self.voltage_data['Q'][ind][idx],
                                                                       self.voltage_data['voltages'][ind][idx])
                if idx != len(path) - 1:
                    data_str += '\n'
        if self.args.get('csv'):
            with open(''.join(self.species) + '_voltage.csv', 'w') as f:
                f.write(data_str)

        print('\nVoltage data:')
        print('\n' + data_str)

    def volume_curve(self):
        """ Take stable compositions and volume and calculate
        volume expansion per "B" in AB binary.

        """

        if self._dimension == 2:
            self._calculate_binary_volume_curve()
        else:
            raise NotImplementedError('Volume curves have only been implemented for binary phase diagrams.')

        data_str = ''
        data_str += '# ' + ''.join(self.species) + '\n'
        data_str += '# {:>10},\t{:>10}\n'.format('x in {d[0]}x{d[1]}'.format(d=self.species),
                                                 'Volume per {} (Ang^3)'.format(self.species[1]))
        for idx, _ in enumerate(self.volume_data['x']):
            data_str += '{:>10.2f},\t{:>10.4f}'.format(self.volume_data['x'][idx],
                                                       self.volume_data['vol_per_y'][idx])
            if idx != len(self.volume_data['x']) - 1:
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
        x = get_num_intercalated(hull_cursor)
        # sort for voltage calculation
        capacities = get_array_from_cursor(hull_cursor, 'gravimetric_capacity')
        capacities = capacities[np.argsort(x)]
        stable_enthalpy_per_b = get_array_from_cursor(hull_cursor,
                                                      self._extensive_energy_key + '_per_b')[np.argsort(x)]

        x = np.sort(x)
        x, uniq_idxs = np.unique(x, return_index=True)
        stable_enthalpy_per_b = stable_enthalpy_per_b[uniq_idxs]
        capacities = capacities[uniq_idxs]
        V = []
        for i, _ in enumerate(x):
            V.append(
                -(stable_enthalpy_per_b[i] - stable_enthalpy_per_b[i-1]) / (x[i] - x[i-1]) + mu_enthalpy[0])
        V[0] = V[1]
        V[-1] = 0
        # make V, Q and x available for plotting
        self.voltage_data['voltages'] = []
        self.voltage_data['voltages'].append(V)
        self.voltage_data['Q'] = []
        self.voltage_data['Q'].append(capacities)
        self.voltage_data['x'] = []
        self.voltage_data['x'].append(x)

    def _calculate_ternary_voltage_curve(self, hull_cursor):
        """ Calculate tenary voltage curve, setting self.voltage_data.
        First pass written by James Darby, jpd47@cam.ac.uk.

        Parameters:
            hull_cursor (list(dict)): list of structures to include in the voltage curve.

        """
        points = np.hstack((
            get_array_from_cursor(hull_cursor, 'concentration'),
            get_array_from_cursor(hull_cursor, self.energy_key).reshape(len(hull_cursor), 1)))
        stoichs = get_array_from_cursor(hull_cursor, 'stoichiometry')
        mu_enthalpy = get_array_from_cursor(self.chempot_cursor, self.energy_key)
        enthalpy_active_ion = mu_enthalpy[0]
        # do another convex hull on just the known hull points, to allow access to useful indices
        import scipy.spatial
        hull = scipy.spatial.ConvexHull(points)

        endpoints = []
        endstoichs = []
        for ind, point in enumerate(points):
            if point[0] == 0 and point[1] != 0 and point[1] != 1:
                if not any([point.tolist() == test_point.tolist() for test_point in endpoints]):
                    endpoints.append(point)
                    endstoichs.append(stoichs[ind])
        print('{} starting point(s) found.'.format(len(endstoichs)))
        for endstoich in endstoichs:
            print(get_formula_from_stoich(endstoich), end=' ')
        print('\n')

        # iterate over possible endpoints of delithiation
        _reactions = []
        _voltages = []
        _capacities = []
        for reaction_ind, endpoint in enumerate(endpoints):
            ratio = endpoint[1] / (1 - endpoint[0] - endpoint[1])

            print(30 * '-')
            print('Reaction {}, {}:'.format(reaction_ind, get_formula_from_stoich(endstoichs[reaction_ind])))

            y0 = endpoint[1] / (1 - endpoint[0])
            simp_in = 0
            intersections = []
            crossover = []
            # put starting point into crossover just in case it is not detected
            sum_conc = sum([int(species[1]) for species in endstoichs[reaction_ind]])
            conc = [int(species[1]) / sum_conc for species in endstoichs[reaction_ind]]
            if len(conc) == 2:
                conc.insert(0, 0.0)
            crossover.append(conc)
            for simplex in hull.simplices:
                tints = []
                for i in range(3):
                    j = (i + 1) % 3
                    e = points[simplex[i], 0]
                    f = points[simplex[i], 1]
                    g = points[simplex[j], 0] - points[simplex[i], 0]
                    h = points[simplex[j], 1] - points[simplex[i], 1]

                    x1 = e
                    y1 = f
                    z1 = 1 - x1 - y1
                    x2 = points[simplex[j], 0]
                    y2 = points[simplex[j], 1]
                    z2 = 1 - x2 - y2

                    if np.abs(h + g * y0) > EPS:
                        tin = (e * h + g * y0 - f * g) / (h + g * y0)
                        s2 = (y0 - e * y0 - f) / (h + g * y0)
                        if tin >= 0 and tin <= 1 and s2 >= 0 and s2 <= 1:
                            tints = np.append(tints, tin)
                            a = 1
                            # x1-x2 never == 0 on points we care about
                            if np.abs(x1 - x2) > EPS:
                                b = (y1 - y2) / (x1 - x2)
                                c = (z1 - z2) / (x1 - x2)
                                x_cross = tin
                                y_cross = b * (tin - x1) / a + y1
                                z_cross = c * (tin - x1) / a + z1
                                # only append unique points
                                if (len(crossover) == 0 or not np.any([np.isclose([x_cross, y_cross, z_cross], val)
                                                                       for val in crossover])):
                                    if y1 != 0 and y2 != 0 and round(float(z1 / y1), 5) == round(float(
                                            z2 / y2), 5) and round(float(z1 / y1), 5) == round(ratio, 5):
                                        pass
                                    else:
                                        crossover.append([x_cross, y_cross, z_cross])
                if len(tints) != 0:
                    temp = [simp_in, np.amin(tints), np.amax(tints)]
                    # condition removes the big triangle and the points which only graze the line of interest
                    if all([temp[2] > EPS, temp[1] < 1, temp[2] - temp[1] > EPS, temp[2] - temp[1] < 1,
                            temp[1] != temp[2]]):
                        intersections = np.append(intersections, temp)
                simp_in += 1

            # if tie line runs from fully de-lithiated to pure lithium (for example), then print and skip
            if len(intersections) == 0:
                print_notify('No intermediate structures found for starting point {}.'.format(get_formula_from_stoich(endstoichs[reaction_ind])))
                continue

            intersections = np.asarray(intersections)
            intersections = intersections.reshape(-1, 3)
            intersections = intersections[intersections[:, 1].argsort()]
            ends_of_rows = []
            min_values = []
            rows_to_keep = []
            # remove row corresponding to largest triangle, i.e. chempots only, and near duplicates (i.e. points with multiple tie-lines)
            for ind, row in enumerate(intersections):
                if not (row[1:].tolist() == [0, 1] or row[1:].tolist() in ends_of_rows or
                        np.any(np.isclose(row.tolist()[1], [val for val in min_values]))):
                    rows_to_keep.append(ind)
                    ends_of_rows.append(row[1:].tolist())
                    min_values.append(row.tolist()[1])
            intersections = intersections[rows_to_keep]

            voltages = []
            crossover = sorted(crossover)
            capacities = sorted([get_generic_grav_capacity(point, self.species) for point in crossover])
            reactions = []
            reaction = [get_formula_from_stoich(endstoichs[reaction_ind])]
            reactions.append(reaction)
            for ind, face in enumerate(intersections):
                simplex_index = int(face[0])
                reaction = []
                reaction = [get_formula_from_stoich(hull_cursor[idx]['stoichiometry'])
                            for idx in hull.simplices[simplex_index]
                            if get_formula_from_stoich(hull_cursor[idx]['stoichiometry']) not in reaction]
                reactions.append(reaction)
                print('{d[0]} + {d[1]} + {d[2]}'.format(d=reaction))
                energy_vec = points[hull.simplices[simplex_index], 2]
                comp = points[hull.simplices[simplex_index], :]
                comp[:, 2] = 1 - comp[:, 0] - comp[:, 1]

                comp = comp.T
                comp_inv = np.linalg.inv(comp)

                V = -(comp_inv.dot([1, 0, 0])).dot(energy_vec)
                V = V + enthalpy_active_ion
                # double up on first voltage
                if ind == 0:
                    voltages.append(V)
                if ind != len(intersections) - 1:
                    print(5 * (ind + 1) * ' ' + ' ---> ', end='')
                voltages.append(V)

            _reactions.append(reactions)
            _capacities.append(capacities)
            _voltages.append(voltages)
            print('\n')
        assert len(_capacities) == len(_voltages)

        self.voltage_data['x'] = _capacities
        self.voltage_data['Q'] = _capacities
        self.voltage_data['voltages'] = _voltages
        self.voltage_data['reactions'] = _reactions
        self.voltage_data['endstoichs'] = endstoichs

    def _calculate_binary_volume_curve(self):
        """ Take stable compositions and volume and calculate volume
        expansion per "B" in AB binary.

        """
        stable_comp = get_array_from_cursor(self.hull_cursor, 'concentration')
        for doc in self.hull_cursor:
            if 'cell_volume_per_b' not in doc:
                print(doc)
        stable_vol = get_array_from_cursor(self.hull_cursor, 'cell_volume_per_b')
        self.volume_data['x'] = np.asarray([comp / (1 - comp) for comp in stable_comp[:-1]]).flatten()
        self.volume_data['vol_per_y'] = np.asarray([vol for vol in stable_vol[:-1]])
        self.volume_data['bulk_volume'] = self.volume_data['vol_per_y'][0]
        self.volume_data['bulk_species'] = self.species[-1]
        self.volume_data['volume_ratio_with_bulk'] = self.volume_data['vol_per_y'] / self.volume_data['bulk_volume']

    def _set_species(self, species=None, elements=None):
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

        self.species = species
