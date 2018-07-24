# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the "TDHull" class
for assessing phase stability from finite temperature free energies.

Written by Angela F. Harper afh41@cam.ac.uk

"""

import numpy as np


from matador.hull.hull import QueryConvexHull
from matador.plotting.plotting import plotting_function
from matador.utils.chem_utils import get_formula_from_stoich


class TDHull(object):
    """ Use QueryConvexHull to construct several hulls at different temperatures,
    based on the free_energy in those hulls rather than the enthalpy.
    This is implemented only for a list of structures, i.e. a folder
    of castep files converted to a cursor using something akin to:

    `cursor, success = castep2dict(glob('/path/to/folder/*.castep'), db=False)`.

    """

    def __init__(self, elements=None, cursor=None, subcmd='hull',
                 temperature_list=None, energy_key='free_energy', **kwargs):
        """Initialize a class from a cursor (list of matador dicts and construct
        a temperature dependent phase diagram or voltage profile. Only 2D.

        Keyword Arguments:
            cursor (list(dict)): specify list of matador docs
            subcmd (str): either 'hull' or 'voltage' FIXME voltage is not implemented
            temperatures (list): list of temperatures to calculate at
            energy_key (str): type of energy to use (enthalpy or free_energy)
            kwargs (dict): arguments from traditional matador options

        """

        # add in kwargs no_temp_plot to override no_plot since it is referenced below
        if 'no_plot' in kwargs:
            kwargs['no_temp_plot'] = kwargs['no_plot']
            del kwargs['no_plot']
        self.args = kwargs

        self.cursor = cursor
        if self.args.get('subcmd') is None:
            self.args['subcmd'] = subcmd

        self.temperature_list = temperature_list
        self.energy_key = energy_key
        self.subcmd = subcmd
        self.elements = elements
        self.thermo_energy_key = 'thermo_' + self.energy_key
        self.thermo_energy_key_pa = 'thermo_'+self.energy_key + '_per_atom'
        if self.cursor is None:
            raise RuntimeError('No structures in cursor!')

        if self.subcmd == 'hull':
            if self.temperature_list is None:
                # plot one hull at 0K first and here energy_key will become energy_key_per_atom
                print('Plotting hull using %s' % self.energy_key)
                QueryConvexHull(cursor=self.cursor, energy_key=self.energy_key, **kwargs)
            else:
                print('Plotting temperature dependent hull...')
                self.plot_td_hull(**kwargs)
        if self.subcmd == 'voltage':
            if self.temperature_list is None:
                print('Plotting voltage at 0K using %s' % self.energy_key)
                QueryConvexHull(cursor=self.cursor, energy_key=self.energy_key, **kwargs)
            else:
                print('Plotting temperature dependent voltage curve...')
                self.plot_td_voltage(**kwargs)

    @plotting_function
    def plot_td_hull(self, **kwargs):
        """ Plot 2D Temperature Dependent hull at temperatures
        in temperature_list
        """
        from matador.plotting import plot_2d_hull
        import matplotlib.pyplot as plt

        _, ax = plt.subplots(figsize=(8, 6))
        temperature_colors = plt.cm.get_cmap('brg')(np.linspace(0, 1, (len(self.temperature_list)+1)*2))

        hull = QueryConvexHull(no_plot=True,
                               cursor=self.cursor,
                               energy_key=self.energy_key,
                               elements=self.elements,
                               **kwargs)
        current_color = temperature_colors[0]

        # plot static hull first (enthalpy)
        ax = plot_2d_hull(hull, ax=ax, show=False, plot_points=False)
        self.plot_hull_line(hull=hull, ax=ax, label='Static Lattice', color='k')

        # get zero point and temperature dependent energies
        self.set_zp_td_energy()

        # plot 0K hull
        hull = QueryConvexHull(no_plot=True, cursor=self.cursor,
                               energy_key=self.energy_key, elements=self.elements, **kwargs)

        current_color = temperature_colors[0]
        self.plot_hull_line(hull=hull, ax=ax, label='0 K', color=current_color)

        # plot temperature dependent hulls
        for temp_ind, temperature1 in enumerate(sorted(self.temperature_list)):
            if self.args.get('quiet') is not True:
                print('Temperature is {}'.format(temperature1))
            current_color = temperature_colors[temp_ind+1]
            hull = QueryConvexHull(no_plot=True, cursor=self.cursor,
                                   energy_key=self.thermo_energy_key_pa,
                                   temperature=temperature1, elements=self.elements,
                                   **kwargs)
            self.plot_hull_line(hull=hull, ax=ax, label='{} K'.format(temperature1), color=current_color)

        eform_limits = (np.min(hull.structures[:, 1]), np.max(hull.structures[:, 1]))
        lims = (-0.1 if eform_limits[0] >= 0 else eform_limits[0] - 0.15,
                eform_limits[1] if eform_limits[0] >= 0 else 0.1)
        ax.set_ylim(lims)

        ax.set_yticklabels(['{:.2f}'.format(val) for val in ax.get_yticks()])
        ax.legend(loc=9)

        if self.args.get('savefig'):
            plt.savefig('{d[0]}{d[1]}_tdhull.png'.format(d=self.elements), bbox_inches='tight')

        if not self.args.get('no_temp_plot'):
            plt.show()

    @plotting_function
    def plot_td_voltage(self, **kwargs):
        """ Plot 2D Temperature Dependent hull at temperatures
        in temperature_list
        """
        import matplotlib.pyplot as plt

        _, axQ = plt.subplots(figsize=(8, 6))
        temperature_colors = plt.cm.get_cmap('brg')(np.linspace(0, 1, (len(self.temperature_list)+1)*2))

        # first plot static line
        hull = QueryConvexHull(no_plot=True,
                               cursor=self.cursor,
                               energy_key=self.energy_key,
                               elements=self.elements,
                               **kwargs)

        self.plot_voltage_line(hull, axQ=axQ, label='Static Lattice', color='k')

        # get zero point corrected energy
        self.set_zp_td_energy()

        # plot 0K hull next with zero point correction
        hull = QueryConvexHull(no_plot=True,
                               cursor=self.cursor,
                               energy_key=self.energy_key,
                               elements=self.elements,
                               **kwargs)

        current_color = temperature_colors[0]
        self.plot_voltage_line(hull, axQ=axQ, label='0 K', color=current_color)

        # plot temperature dependent curves
        for temp_ind, temperature1 in enumerate(sorted(self.temperature_list)):
            if not self.args.get('quiet'):
                print('Temperature is %f' % temperature1)
            current_color = temperature_colors[temp_ind+1]
            hull = QueryConvexHull(no_plot=True,
                                   cursor=self.cursor,
                                   energy_key=self.thermo_energy_key_pa,
                                   temperature=temperature1,
                                   elements=self.elements,
                                   **kwargs)

            self.plot_voltage_line(hull, axQ=axQ, label='{} K'.format(temperature1), color=current_color)

        # set axes labels
        axQ.set_ylabel('Voltage (V) vs {}$^+$/{}'.format(hull.elements[0], hull.elements[0]))
        axQ.set_xlabel('Gravimetric cap. (mAh/g)')
        _, end = axQ.get_ylim()
        axQ.set_ylim(0, 1.1 * end)
        _, end = axQ.get_xlim()
        axQ.set_xlim(0, 1.1 * end)
        axQ.grid('off')
        plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)

        axQ.legend(loc='upper right')

        if self.args.get('savefig'):
            plt.savefig('{d[0]}{d[1]}_tdvoltage.png'.format(d=self.elements),
                        bbox_inches='tight', dpi=300)

        if not self.args.get('no_temp_plot'):
            plt.show()

    @staticmethod
    def plot_hull_line(hull, ax, label='', color=None, plot_points=True):
        """ Add one convex hull tie-line to a given set of axes.

        Parameters:
            hull (matador.QueryConvexHull): the hull to add to the plot.
            ax (matplotlib.Axes): the axes object add the line to.

        Keyword arguments:
            label (str): the label to attach to the hull.
            color (str): the color to use for the hull.
            plot_points (bool): whether to plot all points on the hull,
                or just those on the tie-line.

        """

        tie_line = hull.structure_slice[hull.hull.vertices]
        if plot_points:
            ax.scatter(hull.structures[np.argsort(hull.hull_dist), 0][::-1],
                       hull.structures[np.argsort(hull.hull_dist), -1][::-1],
                       lw=1, alpha=1, edgecolor='k', zorder=10000,
                       c=color)
        ax.plot(np.sort(tie_line[:, 0]), tie_line[np.argsort(tie_line[:, 0]), 1],
                c=color, lw=2, alpha=1, marker='o', markeredgecolor='k', label=label,
                zorder=99999)

        return ax

    def set_zp_td_energy(self):
        """ Set the zero-point energy and temperature-dependent
        energy fields for each document in the cursor.

        """
        for doc in self.cursor:
            doc[self.energy_key] += doc['zero_point_E']
            doc[self.energy_key + '_per_atom'] += (doc['zero_point_E'] / doc['num_atoms'])
        # calculate thermo_free_energy_per_atom for each doc in cursor
        for doc in self.cursor:
            doc[self.thermo_energy_key_pa] = {}
            for temp in doc['thermo_temps']:
                # calculate the thermo_free_energy_per_atom
                doc[self.thermo_energy_key_pa][temp] = (doc[self.energy_key + '_per_atom'] +
                                                        doc[self.thermo_energy_key][temp] / doc['num_atoms'])
                doc[self.thermo_energy_key][temp] += doc[self.energy_key]

    def plot_voltage_line(self, hull, axQ, label='', color='k'):
        """ Add one voltage curve to a given set of axes.

        Parameters:
            hull (matador.QueryConvexHull): the hull to add to the plot.
            axQ (matplotlib.Axes): the axes object add the line to.

        Keyword arguments:
            label (str): the label to attach to the voltage curve.
            color (str): the color to use for the curve.

        """

        for ind, voltage in enumerate(hull.voltage_data['voltages']):
            for i in range(1, len(voltage) - 1):
                if i == 1 and hull.args.get('expt'):
                    axQ.plot([hull.voltage_data['Q'][ind][i-1], hull.voltage_data['Q'][ind][i]],
                             [voltage[i], voltage[i]],
                             marker='*',
                             lw=2,
                             c=color,
                             label='DFT (this work)')
                elif i == 1 and len(hull.voltage_data['voltages']) != 1:
                    axQ.plot([hull.voltage_data['Q'][ind][i-1], hull.voltage_data['Q'][ind][i]],
                             [voltage[i], voltage[i]],
                             marker='o',
                             markersize=5,
                             lw=2,
                             c=color,
                             label=get_formula_from_stoich(hull.endstoichs[ind], tex=True))
                elif i == 1:
                    axQ.plot([hull.voltage_data['Q'][ind][i-1], hull.voltage_data['Q'][ind][i]],
                             [voltage[i], voltage[i]],
                             marker='o',
                             markersize=5,
                             lw=2,
                             c=color, label=label)
                else:
                    axQ.plot([hull.voltage_data['Q'][ind][i-1], hull.voltage_data['Q'][ind][i]],
                             [voltage[i], voltage[i]],
                             marker='o',
                             markersize=5,
                             lw=2,
                             c=color)
                if i != len(voltage) - 2:
                    axQ.plot([hull.voltage_data['Q'][ind][i], hull.voltage_data['Q'][ind][i]],
                             [voltage[i], voltage[i + 1]],
                             marker='o',
                             markersize=5,
                             lw=2,
                             c=color)

            if hull.args.get('labels'):
                eps = 1e-9
                hull.label_cursor = [doc for doc in hull.hull_cursor if doc['hull_distance'] <= 0 + eps]
                hull.label_cursor = hull.label_cursor[1:-1]
                for i in range(len(hull.label_cursor)):
                    axQ.annotate(get_formula_from_stoich(hull.label_cursor[i]['stoichiometry'],
                                                         elements=hull.elements, tex=True),
                                 xy=(hull.voltage_data['Q'][0][i+1]+0.02*max(hull.voltage_data['Q'][0]),
                                     hull.voltage_data['voltages'][0][i+1]+0.02*max(hull.voltage_data['voltages'][0])),
                                 textcoords='data',
                                 ha='center',
                                 zorder=9999)
            if hull.args.get('expt') or len(hull.voltage_data['voltages']) != 1:
                axQ.legend(loc=1)
