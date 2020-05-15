# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functionality to fit E(V) equations
of state to ab initio data, primarily to calculate bulk moduli.

"""

import re

import numpy as np
from matador.utils.chem_utils import eV_PER_ANGSTROM_CUBED_TO_GPa


def get_equation_of_state(seed, plot=False):
    """ Extract E(V) data from CASTEP files and perform
    fits for the equation of state and bulk moduli.

    Parameters:
        seed (str/dict): filename or scraped dictionary to fit.

    Keyword arguments:
        plot (bool): plot the fitted EoS.

    """
    if isinstance(seed, dict):
        results = seed
    else:
        from matador.scrapers import castep2dict
        results, success = castep2dict(seed, intermediates=True, db=False)
        if not success:
            raise RuntimeError(results)

    volumes = []
    energies = []
    for snapshot in results['intermediates']:
        volumes.append(snapshot['cell_volume'])
        energies.append(snapshot['total_energy'])
    volumes.append(results['cell_volume'])
    energies.append(results['total_energy'])
    volumes = np.asarray(volumes)
    energies = np.asarray(energies)

    if len(volumes) < 3:
        raise RuntimeError('Seed {} does not contain enough energies vs volumes to fit a bulk modulus')

    energies = energies[np.argsort(volumes)]
    volumes = np.sort(volumes)

    E_0 = np.min(energies)
    V_0 = volumes[np.argmin(energies)]

    types_of_fit = EquationOfState.__subclasses__()

    results['eos'] = []
    results['summary'] = []
    probes_all = []
    curves = []
    labels = []

    for eos_type in types_of_fit:
        eos = eos_type(E_0, V_0)
        eos.fit(volumes, energies)

        results['summary'].append(eos_type.__name__ + '\n')
        results['summary'].append('fitting parameters = {d[0]:.6f}, {d[1]:.6f}\n'.format(d=eos.fit_parameters))
        results['summary'].append('rrrmsd = {:10.10f} %\n'.format(eos.rrmsd))
        results['summary'].append('bulk modulus = {d[0]:.6f} +/- {d[1]:.6f} GPa\n'
                                  .format(d=(eos.bulk_modulus, eos.bulk_modulus_err)))
        results['summary'].append(80*'-' + '\n')
        probes = np.linspace(min(volumes), max(volumes), num=100)
        fitted_curve = eos.evaluate(probes, *eos.popt)
        results['eos'].append(eos)
        probes_all.append(probes)
        curves.append(fitted_curve)
        name = eos_type.__name__.replace('EOS', '')
        label = '-'.join(re.findall('[A-Z][^A-Z]*', name))
        labels.append(label)

    if plot:
        import matplotlib.pyplot as plt
        _, ax = plt.subplots(1, figsize=(8, 6))
        ax.set_ylabel('Total energy (eV)')
        ax.set_xlabel('Cell volume ($\\AA^3$)')
        ls = ['--', '-.', ':']
        for ind, probe in enumerate(probes_all):
            ax.plot(probes_all[ind], curves[ind], label=labels[ind], ls=ls[ind])
        ax.plot(volumes, energies, marker='o', label='DFT data')

        ax.legend(loc=0)
        plt.savefig(seed + '.png')

    return results


class EquationOfState:
    """ Abstract class for E(V) isothermal equations of state. Used to
    perform least squares fitting to arbitrary functional form.

    Attributes:
        E_0 (float): equilbirium lattice energy.
        V_0 (float): equilibrium lattice volume.
        popt (:obj:`list` of :obj:`float`): fitting parameters for
            particular EOS.
        p0 (:obj:`list` of :obj:`float`): initial guess at fitting parameters.
        rrmsd (float): relative root mean square deviation of fit, as
            defined by [2].

    """

    def __init__(self, E_0, V_0):
        """ Set up EOS ready for fit.

        Parameters:
            E_0 (float): equilibrium energy.
            V_0 (float): equilibrium volume.

        """
        self.E_0 = E_0
        self.V_0 = V_0
        self.popt = None
        self.rrmsd = None
        self.p0 = [0.01, 1]

    def fit(self, volumes, energies):
        """ Perform a least squares fit on the volumes and energies.

        Parameters:
            volumes (:obj:`list` of :obj:`float`): calculated volumes.
            energies (:obj:`list` of :obj:`float`): calculated energies.

        """
        from scipy import optimize
        self.popt, _ = optimize.leastsq(self.residual, self.p0, args=(volumes, energies))
        self.rrmsd = np.sqrt(np.sum(((energies - self.evaluate(volumes, *self.popt)) / energies)**2) / (len(energies) - 1))

    def evaluate(self, V, B, C):
        """ Evaluate the EOS.

        Parameters:
            V (:obj:`list` of :obj:`float`): volumes at which to test.
            B (float): the first fitting parameter defined in [2].
            C (float): the second fitting parameter defined in [2].

        """
        raise NotImplementedError

    def residual(self, guess, V, E):
        """ Calculate the resdiual of the current fit.

        Parameters:
            guess (:obj:`list` of :obj:`float`): interim fitting parameters.
            V (:obj:`list` of :obj:`float`): volumes at which to test.
            E (:obj:`list` of :obj:`float`): energies at the volumes V.

        """
        return E - self.evaluate(V, *guess)

    @property
    def bulk_modulus(self):
        """ Returns the bulk modulus predicted by the fit. """
        bulk_modulus = self.get_bulk_modulus()
        return bulk_modulus

    @property
    def bulk_modulus_err(self):
        """ Returns the estimated error on the bulk modulus. """
        return self.bulk_modulus * self.rrmsd

    @property
    def fit_parameters(self):
        """ Return list of final fitting parameters. """
        return self.popt


class BirchMurnaghanEulerianEOS(EquationOfState):
    """ Implements the 3rd order Birch-Murnaghan EOS [1] for given data, as
    provided by Ref. [2] in the Eulerian frame.

    [1]. Francis Birch, Phys. Rev. 71 809 (1947)
         DOI: 10.1103/PhysRev.71.809.

    [2]. K. Latimer, S. Dwaraknath, K. Mathew, D. Winston, K. A. Persson,
         npj Comput. Mater. 2018 41 2018, 4, 40.
         DOI: 10.1038/s41524-018-0091-x.

    """

    def evaluate(self, V, B, C):
        nu = V / self.V_0
        nu = (nu**(-2.0/3.0) - 1)
        return self.E_0 + B * self.V_0 * (nu**2 + 0.5 * C * nu**3)

    def get_bulk_modulus(self):
        """ Return the bulk modulus of this particular fit. """
        if self.popt is None:
            raise RuntimeError('No fit performed.')
        return self.popt[0] * 8.0/9.0 * eV_PER_ANGSTROM_CUBED_TO_GPa


class PoirerTarantolaEOS(EquationOfState):
    """ Implements the logarithmic Poirer-Tarantola [3] for given data, as
    provided by Ref. [2].

    [2]. K. Latimer, S. Dwaraknath, K. Mathew, D. Winston, K. A. Persson,
         npj Comput. Mater. 2018 41 2018, 4, 40.
         DOI: 10.1038/s41524-018-0091-x.

    [3]. J. Poirer, A. Tarantola, Phys. Earth Planet. Inter. 109 1-8 (1998).

    """

    def evaluate(self, V, B, C):
        nu = V / self.V_0
        return self.E_0 + B * self.V_0 * np.log(nu)**2 * (3 - C*np.log(nu))

    def get_bulk_modulus(self):
        """ Return the bulk modulus of this particular fit. """
        if self.popt is None:
            raise RuntimeError('No fit performed.')
        return self.popt[0] * 6 * eV_PER_ANGSTROM_CUBED_TO_GPa


class TaitEOS(EquationOfState):
    """ Implements the exponential Tait EOS [4] for given data, as
    provided by Ref. [2].

    [2]. K. Latimer, S. Dwaraknath, K. Mathew, D. Winston, K. A. Persson,
         npj Comput. Mater. 2018 41 2018, 4, 40.
         DOI: 10.1038/s41524-018-0091-x.

    [4]. J. H. Dymond, R. Malhotra, Int. J. Thermophys. 9 941-951 (1988).

    """

    def evaluate(self, V, B, C):
        nu = V / self.V_0
        return self.E_0 + (B*self.V_0 / C) * (nu - 1 + (1/C) * (np.exp(C*(1-nu)) - 1))

    def get_bulk_modulus(self):
        """ Return the bulk modulus of this particular fit. """
        if self.popt is None:
            raise RuntimeError('No fit performed.')
        return self.popt[0] * eV_PER_ANGSTROM_CUBED_TO_GPa
