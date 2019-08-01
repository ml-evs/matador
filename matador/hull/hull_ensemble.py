# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements parameterised phase diagrams. """


from matador.hull import PhaseDiagram, QueryConvexHull
from matador.utils.cursor_utils import filter_cursor_by_chempots, recursive_get, recursive_set, set_cursor_from_array
from matador.utils.chem_utils import get_formation_energy, get_root_source


class EnsembleHull(QueryConvexHull):
    """ Class to create and store an ensemble of composition vs energy convex
    hulls from cursor data. The variable energies must be stored under a given
    key, e.g. `doc['_beef'][energy_key][beef_index]`, as specified by init.

    Data must be stored in the following way under each document in cursor:

        ```
        {...,
         data_key: {
                    parameter_key: <list of parameter values>,
                    energy_key: <list of energies at parameter values>,
                   },
         ...
        }
        ```
    Hull data will be stored as arrays per document under `doc[data_key]['hull_distance']`
    and ``doc[data_key]['formation_' + energy_key]``.

    Inherits the attributes of matador.hull.QueryConvexHull, with many set to
    None.

    Attributes:
        phase_diagrams (list of :obj:`matador.hull.PhaseDiagram`): list of phase diagram
            objects for each parameter value.

    """
    def __init__(self, cursor, data_key, energy_key='enthalpy_per_atom', num_samples=None,
                 parameter_key=None, species=None, subcmd='hull', verbosity=None, **kwargs):
        """ Initialise EnsembleHull from a cursor, with other keywords
        following QueryConvexHull.

        Parameters:
            cursor (list[dict]): list of matador documents containing
                variable parameter data for energies.
            data_key (str): the key under which all parameter data is
                stored to the variable parameter, e.g. `_beef` or `_temperature`.

        Keyword arguments:
            energy_key (str): the key under `parameter_key` to use to create
                the hulls.
            parameter_key (str): the key pertaining to the variable parameter
                itself, e.g. `temperature` or `thetas`.
            num_samples (int): use up to this many samples in creating the hull.
            species (list[str]): list of elements/chempots to use, in
                the desired order.
            subcmd (str): either 'hull' or 'voltage'.
            plot_kwargs (dict): arguments to pass to plot_hull function.
            kwargs (dict): other arguments to pass to QueryConvexHull.

        """
        super().__init__(cursor=cursor,
                         energy_key=energy_key,
                         species=species,
                         subcmd=subcmd,
                         no_plot=True,
                         lazy=False,
                         **kwargs)

        if self.phase_diagram is None:
            del self.phase_diagram
        if self.structure_slice is None:
            del self.structure_slice
        if self.hull_dist is None:
            del self.hull_dist

        self.from_cursor = True
        self.verbosity = verbosity

        # set up relative keys
        self.formation_key = 'formation_' + self.energy_key
        if isinstance(data_key, str):
            self.data_key = [data_key]
        else:
            self.data_key = data_key
        if parameter_key is None:
            self.parameter_key = None
        else:
            self.parameter_key = self.data_key + [parameter_key]
        self.formation_key = self.data_key + [self.formation_key]
        self.hulldist_key = self.data_key + ['hull_distance']
        self.energy_key = self.data_key + [self.energy_key]

        self.phase_diagrams = []

        self.set_chempots()
        self.cursor = sorted(filter_cursor_by_chempots(self.species, self.cursor),
                             key=lambda doc: (recursive_get(doc, self.energy_key), doc['concentration']))

        if parameter_key is None:
            parameter_iterable = recursive_get(self.chempot_cursor[0], self.energy_key)
        else:
            parameter_iterable = recursive_get(self.chempot_cursor[0], self.parameter_key)

        print('Found {} entries under data key: {}.'.format(len(parameter_iterable), self.data_key))

        # allocate formation energy and hull distance arrays
        for ind, doc in enumerate(self.cursor):
            recursive_set(doc, self.formation_key, [None] * len(recursive_get(doc, self.energy_key)))
            recursive_set(doc, self.hulldist_key, [None] * len(recursive_get(doc, self.energy_key)))

        n_hulls = len(parameter_iterable)
        if num_samples is not None:
            parameter_iterable = parameter_iterable[:num_samples]
            print('Using {} out of {} possible phase diagrams.'.format(num_samples, n_hulls))
        else:
            num_samples = n_hulls

        for param_ind, parameter in enumerate(parameter_iterable):
            for ind, doc in enumerate(self.cursor):
                if parameter_key is not None:
                    assert recursive_get(doc, self.parameter_key + [param_ind]) == parameter

                formation_energy = get_formation_energy(self.chempot_cursor, doc,
                                                        energy_key=self.energy_key + [param_ind])
                recursive_set(self.cursor[ind], self.formation_key + [param_ind], formation_energy)
            self.phase_diagrams.append(PhaseDiagram(self.cursor,
                                                    self.formation_key + [param_ind],
                                                    self._dimension))
            set_cursor_from_array(self.cursor,
                                  self.phase_diagrams[-1].hull_dist,
                                  self.hulldist_key + [param_ind])

        self.stability_histogram = self.generate_stability_statistics()

    def generate_stability_statistics(self):
        """ Creates a histogram that counts how many times each structure
        is found to be stable in the ensemble.

        """
        from collections import defaultdict
        histogram = defaultdict(int)
        for pd in self.phase_diagrams:
            for doc in pd.stable_structures:
                histogram[get_root_source(doc)] += 1
        return histogram

    def plot_hull(self, **kwargs):
        """ Hull plot helper function. """
        from matador.plotting.hull_plotting import plot_ensemble_hull
        plot_ensemble_hull(self, '_beef', **kwargs)
