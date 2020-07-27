# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements the base class for parameterised phase diagrams. """

import tqdm
from matador.hull import PhaseDiagram, QueryConvexHull
from matador.utils.cursor_utils import filter_cursor_by_chempots, recursive_get, recursive_set, set_cursor_from_array
from matador.utils.chem_utils import get_formation_energy, get_root_source, get_formula_from_stoich


class EnsembleHull:
    """ Class to create and store an ensemble of composition vs energy convex
    hulls from cursor data. The variable energies must be stored under a given
    key, e.g. `doc['_beef'][energy_key][beef_index]`, as specified by init.

    If `calculate_formation_energy` is False, this class will assume that the cursor
    already contains the appropriate formation energies.

    Data must be stored in the following way under each document in cursor::

        {...,
         data_key: {
            parameter_key: "<list of parameter values>",
            energy_key: "<list of energies at parameter values>",
         },
        }

    Hull data will be stored as arrays per document under `doc[data_key]['hull_distance']`
    and ``doc[data_key]['formation_' + energy_key]``.

    Attributes:
        phase_diagrams (list of :obj:`matador.hull.PhaseDiagram`): list of phase diagram
            objects for each parameter value.

    """
    def __init__(self, cursor=None, data_key=None, hull=None, energy_key='enthalpy_per_atom', chempot_energy_key=None,
                 calculate_formation_energy=True, num_samples=None, parameter_key=None, species=None,
                 voltage=False, verbosity=None, **kwargs):
        """ Initialise EnsembleHull from a cursor, with other keywords
        following QueryConvexHull.

        Keyword arguments:
            cursor (list[dict]): list of matador documents containing
                variable parameter data for energies.
            data_key (str): the key under which all parameter data is
                stored to the variable parameter, e.g. `_beef` or `_temperature`.
            hull (QueryConvexHull): hull object to use instead of cursor.
            energy_key (str): the key under `parameter_key` to use to create
                the hulls.
            chempot_energy_key (str): the key used to create the first convex hull.
            parameter_key (str): the key pertaining to the variable parameter
                itself, e.g. `temperature` or `thetas`.
            calculate_formation_energy (bool): if False, it will be assumed that
                formation energies can already be found under the `data_key`.
            num_samples (int): use up to this many samples in creating the hull.
            species (list[str]): list of elements/chempots to use, in
                the desired order.
            voltage (bool): whether or not to compute voltage curves.
            plot_kwargs (dict): arguments to pass to plot_hull function.
            kwargs (dict): other arguments to pass to QueryConvexHull.

        """
        # sometimes the first hull needs to be made with a different key
        if chempot_energy_key is not None:
            self.chempot_energy_key = chempot_energy_key
        else:
            self.chempot_energy_key = energy_key

        if cursor is not None:
            self.cursor = cursor

        if hull is None:
            print("Creating hull")
            self.hull = QueryConvexHull(
                cursor=self.cursor,
                energy_key=self.chempot_energy_key,
                species=species,
                voltage=voltage,
                no_plot=True,
                lazy=False,
                **kwargs
            )
        else:
            self.hull = hull
            if cursor is None:
                self.cursor = self.hull.cursor

        self.energy_key = energy_key

        self.species = self.hull.species
        self.dimension = self.hull._dimension

        self.verbosity = verbosity
        # set up relative keys
        self.formation_key = 'formation_' + self.energy_key
        self.data_key = data_key
        self.parameter_key = parameter_key

        if self.parameter_key is None:
            self._parameter_keys = None
        else:
            self._parameter_keys = [self.data_key] + [parameter_key]
        self._formation_keys = [self.data_key] + [self.formation_key]
        self._hulldist_keys = [self.data_key] + ['hull_distance']
        self._energy_keys = [self.data_key] + [self.energy_key]

        self.phase_diagrams = []

        if self.parameter_key is None:
            parameter_iterable = recursive_get(self.cursor[0], self._energy_keys)
        else:
            parameter_iterable = recursive_get(self.cursor[0], self._parameter_keys)

        print('Found {} entries under data key: {}.'.format(len(parameter_iterable), self.data_key))

        if calculate_formation_energy:
            self.cursor = filter_cursor_by_chempots(self.species, self.cursor)
            self.cursor = sorted(
                self.cursor,
                key=lambda doc: (recursive_get(doc, self.chempot_energy_key), doc['concentration'])
            )

            # allocate formation energy array
            for ind, doc in enumerate(self.cursor):
                recursive_set(doc, self._formation_keys, [None] * len(parameter_iterable))

        # allocate hull distance array
        for ind, doc in enumerate(self.cursor):
            recursive_set(doc, self._hulldist_keys, [None] * len(parameter_iterable))

        n_hulls = len(parameter_iterable)
        if num_samples is not None:
            parameter_iterable = parameter_iterable[:num_samples]
            print('Using {} out of {} possible phase diagrams.'.format(num_samples, n_hulls))
        else:
            num_samples = n_hulls

        for param_ind, parameter in enumerate(tqdm.tqdm(parameter_iterable)):

            if calculate_formation_energy:
                # if formation energies aren't already present, compute them here
                for ind, doc in enumerate(self.cursor):
                    if self.parameter_key is not None:
                        assert recursive_get(doc, self._parameter_keys + [param_ind]) == parameter

                    formation_energy = get_formation_energy(self.hull.chempot_cursor, doc,
                                                            energy_key=self._energy_keys + [param_ind])
                    recursive_set(self.cursor[ind], self._formation_keys + [param_ind], formation_energy)

            # compute phase diagram based on formation energies at each parameter set
            self.phase_diagrams.append(PhaseDiagram(self.cursor,
                                                    self._formation_keys + [param_ind],
                                                    self.dimension))

            # set hull distances with values from each phase diagram
            set_cursor_from_array(self.cursor,
                                  self.phase_diagrams[-1].hull_dist,
                                  self._hulldist_keys + [param_ind])

        self.stability_histogram = self.generate_stability_statistics()

    def generate_stability_statistics(self, group_by='structure'):
        """ Creates a histogram that counts how many times each structure
        is found to be stable in the ensemble.

        Keyword arguments:
            group_by (str): either 'structure' or 'formula' for bar groupings.

        """
        from collections import defaultdict
        histogram = defaultdict(int)
        for pd in self.phase_diagrams:
            for doc in pd.stable_structures:
                if group_by == 'formula':
                    histogram[get_formula_from_stoich(doc['stoichiometry'])] += 1
                else:
                    histogram[get_root_source(doc)] += 1
        return histogram

    def generate_composition_stability(self):
        # Here we are trying to make a bar chart of x-axis chempot where a structure 
        # is stable and y axis bars of that structure
        from collections import defaultdict
        histogram = defaultdict(list)
        chempot = ''
        space = []
        
        # First identify the chemical potential that has been altered, and get the 'space'
        # i.e. range along which it has been altered
        for mu in self.chempot_cursor:
            energies = mu[self.data_key][self.energy_key]
            if (energies[1]-energies[-1])/len(energies) != 0:
                chempot = get_formula_from_stoich(mu.get('stoichiometry'))
                space = energies
                break

        # now create a histogram of all the structures and uner which chempot theyre stable
        for hull in self.phase_diagrams:
            for doc in hull.stable_structures:
                #if the structure is not yet been seen add it to the hist        
                form = get_formula_from_stoich(doc.get('stoichiometry'))
                if doc not in self.chempot_cursor:
                    histogram['%s %s'%(form,doc.get('space_group'))].append(space[hull.formation_key[2]])

        return chempot, histogram
    
    def plot_composition_stability(self, **kwargs):
        """Plots stable structures at each chempot in bar graph form"""
        from matador.plotting.hull_plotting import plot_ensemble_stability
        chempot, hist = self.generate_composition_stability()
        plot_ensemble_stability(self, self.data_key, chempot, hist, **kwargs)

    def plot_hull(self, **kwargs):
        """ Hull plot helper function. """
        from matador.plotting.hull_plotting import plot_ensemble_hull
        return plot_ensemble_hull(self, formation_energy_key=self.formation_key, **kwargs)
