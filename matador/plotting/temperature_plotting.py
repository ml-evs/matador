# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule contains functions to plot temperature-dependent
quantities for particular structures, e.g. free energies and heat
capacities.

"""

from matador.plotting.plotting import plotting_function


@plotting_function
def plot_thermo_curves(seed, show=False, **kwargs):
    """ Plot free energies and heat capacity as a function of temperature, as
    calculated from a CASTEP thermodynamics run in <seed>.castep.

    Parameters:
        seed (str): filename of thermodynamics <seed>.castep file

    Keyword Arguments:
        plot_energy (bool): whether to plot T vs Energy
        plot_heat_cap (bool): whether to plot T vs Cv

    """
    import matplotlib.pyplot as plt
    from matador.scrapers.castep_scrapers import castep2dict
    from matador.utils.chem_utils import AVOGADROS_NUMBER, ELECTRON_CHARGE

    # set defaults
    prop_defaults = {'plot_energy': True, 'plot_heat_cap': True}
    prop_defaults.update(kwargs)
    kwargs = prop_defaults

    if kwargs['plot_energy'] and not kwargs['plot_heat_cap']:
        _, ax_energy = plt.subplots(figsize=(5, 5))
    elif kwargs['plot_energy'] and kwargs['plot_heat_cap']:
        _, ax_grid = plt.subplots(1, 2, figsize=(10, 5), sharey=False,
                                  gridspec_kw={'width_ratios': [1, 1],
                                               'wspace': 0.15,
                                               'left': 0.15})
        ax_energy = ax_grid[0]
        ax_Cv = ax_grid[1]
    elif not kwargs['plot_energy'] and kwargs['plot_heat_cap']:
        _, ax_Cv = plt.subplots(figsize=(5, 5))

    data, success = castep2dict(seed, db=False)
    if not success:
        raise RuntimeError('Failed to read {}: {}'.format(seed, data))

    temps = data['thermo_temps']
    Cv = data['thermo_heat_cap_Cv']
    free_energy = data['thermo_free_energy_F']
    entropy = data['thermo_entropy_S']
    enthalpy = data['thermo_enthalpy_E']
    # E0 = data['zero_point_E']

    entropyT = []
    for i, val in enumerate(entropy):
        entropyT.append(temps[i] * val * (1. / ELECTRON_CHARGE) * (1. / AVOGADROS_NUMBER))  # J/mol/K --> eV

    if kwargs['plot_energy']:
        ax_energy.plot(temps, free_energy, marker='o', lw=1, ls='-', alpha=1, label='Free Energy')
        ax_energy.plot(temps, entropyT, marker='o', lw=1, ls='-', alpha=1, label='Entropy*T')
        ax_energy.plot(temps, enthalpy, marker='o', lw=1, ls='-', alpha=1, label='Enthalpy')
        ax_energy.set_xlabel('Temperature (K)')
        ax_energy.set_title('Energy (eV)')
        ax_energy.legend()

    if kwargs['plot_heat_cap']:
        ax_Cv.plot(temps, Cv, marker='o', lw=1, ls='-', alpha=1, label='Heat Capacity (J/mol/K)')
        ax_Cv.set_title('Heat Capacity (J/mol/K)')
        ax_Cv.set_xlabel('Temperature (K)')

    if any([kwargs.get('pdf'), kwargs.get('svg'), kwargs.get('png')]):
        if kwargs.get('pdf'):
            plt.savefig(seed.replace('.castep', '') + '_thermoplots.pdf', dpi=300, transparent=True)
        if kwargs.get('svg'):
            plt.savefig(seed.replace('.castep', '') + '_thermoplots.svg', dpi=300, transparent=True)
        if kwargs.get('png'):
            plt.savefig(seed.replace('.castep', '') + '_thermoplots.png', dpi=300, transparent=True)

    elif show:
        plt.show()
