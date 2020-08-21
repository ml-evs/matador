# coding: utf-8
# Distributed under the terms of the MIT license.

from typing import Dict, List, Optional
import numpy as np
import warnings

from matador.plotting.plotting import plotting_function
from matador.crystal import Crystal

__all__ = ("MagresReferencer", )


class MagresReferencer:
    """ Class for referencing NMR predictions with experimental data. """

    def __init__(
        self,
        structures_exp: Dict[str, Crystal],
        shifts_exp: List[Dict[str, List[float]]],
        species: str,
        structures: Optional[List[Crystal]] = None,
        warn_tolerance: float = 0.1,
    ):
        self.structures_exp = structures_exp
        self.shifts_exp = shifts_exp
        self.species = species
        self.structures = structures
        self.warn_tolerance = warn_tolerance

        self._calc_shifts = []
        self._expt_shifts = []
        self._fit_weights = []
        self._fit_structures = []
        self._fitted = False

        for formula in self.shifts_exp:
            if self.species not in self.shifts_exp[formula]:
                continue

            if formula not in self.structures_exp:
                warnings.warn(f"Missing {formula} in reference structures.")

            self.match_exp_structure_shifts(
                self.structures_exp[formula], self.shifts_exp[formula][self.species]
            )

        self.fit()
        self.print_fit_summary()

        if self.structures is not None:
            self.set_shifts_from_fit(self.structures)

    def match_exp_structure_shifts(self, structure, shifts):
        relevant_sites = [site for site in structure if site.species == self.species]
        calc_shifts = sorted(
            [site["chemical_shielding_iso"] for site in relevant_sites]
        )
        _shifts = shifts
        if (
            len(_shifts) <= len(relevant_sites)
            and len(relevant_sites) % len(_shifts) == 0
        ):
            multiplier = len(relevant_sites) // len(_shifts)
            _shifts = [shift for cell in [_shifts] * multiplier for shift in cell]
            _weights = [1 / multiplier for shift in _shifts]
        else:
            raise RuntimeError(
                f"Incompatible shift sizes: {len(relevant_sites)} (theor.) vs {len(_shifts)} (expt.), "
                "please provide commensurate cells."
            )

        _shifts = sorted(_shifts, reverse=True)

        self._calc_shifts.extend(calc_shifts)
        self._expt_shifts.extend(_shifts)
        self._fit_weights.extend(_weights)
        self._fit_structures.extend(len(_shifts) * [structure.formula_tex])

        return _shifts, _weights, calc_shifts

    def set_shifts_from_fit(self, structures):
        for ind, struc in enumerate(structures):
            for jnd, site in enumerate(struc):
                if site.species == self.species:
                    site["chemical_shift_iso"] = self.predict(site["chemical_shielding_iso"])[0]

    def fit(self):
        import statsmodels.api as sm

        _fit_shifts = sm.add_constant(self._calc_shifts)
        self.fit_model = sm.regression.linear_model.WLS(self._expt_shifts, _fit_shifts, weights=self._fit_weights)
        self.fit_results = self.fit_model.fit()
        self.fit_intercept = self.fit_results.params[0]
        self.fit_gradient = self.fit_results.params[1]
        self.fit_rsquared = self.fit_results.rsquared
        self._fitted = True

        if abs(self.fit_gradient + 1.0) > self.warn_tolerance:
            warnings.warn(
                f"{self.__class__.__name__} fit gradient was {self.fit_gradient:.2f}, "
                f"outside of tolerated range -1.0 ± {self.warn_tolerance:.2f}"
            )

    def predict(self, shifts):
        _shifts = np.asarray(shifts)
        return self.fit_gradient * _shifts + self.fit_intercept, self.fit_results.bse[1] * _shifts + self.fit_results.bse[0]

    def print_fit_summary(self):
        if self._fitted:
            print("Performed WLS fit for: δ_expt = m * σ_calc + c")
            print(f"m = {self.fit_gradient:3.3f} ± {self.fit_results.bse[1]:3.3f}")
            print(f"c = {self.fit_intercept:3.3f} ± {self.fit_results.bse[0]:3.3f} ppm")
            print(f"R² = {self.fit_rsquared:3.3f}.")
        else:
            raise RuntimeError("Fit has not yet been performed.")

    @plotting_function
    def plot_fit(self, ax=None, padding=100, figsize=None):
        import matplotlib.pyplot as plt
        import seaborn as sns

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        ax.grid(False)
        ax.set_xlim(np.min(self._calc_shifts) - padding, np.max(self._calc_shifts) + padding)
        ax.set_ylim(np.min(self._expt_shifts) - padding, np.max(self._expt_shifts) + padding)
        ax = sns.regplot(
            y=self._expt_shifts,
            x=self._calc_shifts,
            scatter=False,
            ax=ax,
            line_kws={"linewidth": 0},
            color="grey",
            truncate=False,
        )

        sns.scatterplot(
            y=self._expt_shifts,
            x=self._calc_shifts,
            hue=self._fit_structures,
            palette="Dark2",
            ax=ax,
            zorder=100,
        )

        ax.plot(
            np.asarray(ax.get_xlim()),
            self.fit_gradient * np.asarray(ax.get_xlim()) + self.fit_intercept,
            label=f"$m = {self.fit_gradient:3.3f}; c = {self.fit_intercept:3.0f} ppm; R^2 = {self.fit_rsquared:3.3f}$",
            zorder=10,
            lw=1.5,
            c="grey",
        )
        ax.legend(ncol=len(set(self._fit_structures)) // 5)
        ax.set_ylabel("$\\delta_\\mathrm{expt}$ (ppm)")
        ax.set_xlabel("$\\sigma_\\mathrm{calc}$ (ppm)")

        return ax

    @plotting_function
    def plot_fit_and_predictions(self, ax=None, padding=100):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        self.plot_fit(ax=ax, padding=padding)

        for doc in self.structures:
            _calc_shifts = [site["chemical_shielding_iso"] for site in doc if site.species == self.species]
            _predicted_shifts, _predicted_errs = self.predict(_calc_shifts)
            ax.scatter(_predicted_shifts, _calc_shifts, s=5, c='k')
            ax.errorbar(_predicted_shifts, _calc_shifts, fmt='None', xerr=_predicted_errs, lw=0.5, c='k')

        return ax
