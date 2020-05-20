# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the `VariableMuHull` class for assessing
phase stability with variable chemical potentials.

"""

import copy
import itertools
from collections import defaultdict

import numpy as np

from matador.hull import QueryConvexHull
from matador.hull.hull_ensemble import EnsembleHull
from matador.utils.chem_utils import get_formation_energy, get_formula_from_stoich


class VariableMuHull(EnsembleHull):
    """ Leverages `EnsembleHull` to construct variable chemical potential
    hulls from total energies. The procedure is as follows:

    Create a normal QueryConvexHull from the passed cursor, then construct
    the parametererised the formation energies adjusted by the tweaked chemical
    potentials. The parent `EnsembleHull` then recomputes the phase diagrams
    with these new formation energies.

    """

    data_key = "variable_mu"

    def __init__(
        self, cursor=None, hull=None, energy_key='enthalpy_per_atom', variable_chempot_formulae=None,
        chempot_ranges=(-1, 1), chempot_samples=10, **kwargs
    ):

        if hull is None and cursor is not None:
            standard_hull = QueryConvexHull(cursor=cursor, **kwargs)
        elif hull is not None and cursor is None:
            standard_hull = hull
        else:
            raise RuntimeError("Please specify either a hull or a cursor.")

        chempot_values = {
            get_formula_from_stoich(doc['stoichiometry']): doc[energy_key]
            for doc in standard_hull.chempot_cursor
        }

        if variable_chempot_formulae is None:
            variable_chempot_formulae = list(chempot_values.keys())

        self.parameter_key = "mu_offset"
        self.formation_energy_key = "formation_" + energy_key
        self.energy_key = energy_key

        offset_range = np.linspace(*chempot_ranges, chempot_samples)
        chempot_cursor = copy.deepcopy(standard_hull.chempot_cursor)
        _cursor = copy.deepcopy(standard_hull.cursor)

        all_offsets = itertools.product(*len(variable_chempot_formulae) * [offset_range])
        for doc_ind, doc in enumerate(_cursor):
            _cursor[doc_ind][self.data_key] = defaultdict(list)

        for param_ind, offsets in enumerate(all_offsets):
            for mu_ind, bare_mu in enumerate(chempot_cursor):
                for ind, mu in enumerate(variable_chempot_formulae):
                    if get_formula_from_stoich(bare_mu['stoichiometry']) == mu:
                        chempot_cursor[mu_ind][energy_key] = chempot_values[mu] + offsets[ind]
            for doc_ind, doc in enumerate(_cursor):
                for ind, mu in enumerate(variable_chempot_formulae):
                    if get_formula_from_stoich(doc['stoichiometry']) == mu:
                        _cursor[doc_ind][energy_key] = chempot_values[mu] + offsets[ind]
                eform = get_formation_energy(chempot_cursor, doc, energy_key=energy_key)
                _cursor[doc_ind][self.data_key][self.formation_energy_key].append(eform)
                _cursor[doc_ind][self.data_key][self.parameter_key].append(offsets)

        super().__init__(
            hull=hull,
            cursor=_cursor,
            data_key=self.data_key,
            energy_key=self.energy_key,
            calculate_formation_energy=False,
            chempot_energy_key=energy_key,
            parameter_key=self.parameter_key,
            **kwargs
        )

    def plot_hull(self, **kwargs):
        """ Hull plot helper function. """
        from matador.plotting.hull_plotting import plot_ensemble_hull
        return plot_ensemble_hull(self, formation_energy_key=self.formation_key, **kwargs)
