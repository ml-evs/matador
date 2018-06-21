# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the "TDHull" class 
for assessing phase stability from finite temperature free energies.

Written by Angela F. Harper afh41@cam.ac.uk

"""

from matador.hull.hull import QueryConvexHull

class TDHull(object):
    """Use QueryConvexHull to construct several tempearture dependent hulls
    based on the free_energy in those hulls rather than the enthalpy
    This is implemented only for a list of structures, i.e. a folder
    of castep files converted to a cursor using something akin to:
    cursor = [castep2dict(castep,db=False)[0] for castep_file in glob('path/to/folder')]
    """

    def __init__(self, elements=None,cursor=None, subcmd='hull',temperature_list=None, energy_key='free_energy',**kwargs):
        """Initialize a class from a cursor (list of matador dicts and construct 
        a temperature dependent phase diagram or voltage profile. Only 2D.

        Keyword Arguments:
            cursor (list(dict)): specify list of matador docs
            subcmd (str): either 'hull' or 'voltage' FIXME voltage is not implemented
            temperatures (list): list of temperatures to calculate at 
            energy_key (str): type of energy to use (enthalpy or free_energy)
            kwargs (dict): arguments from traditional matador options
        """
        
        #add in kwargs no_temp_plot to override no_plot since it is referenced below
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
            raise RuntimeError('Failed to find structures to create hull!')
        if self.subcmd == 'hull':
            if self.temperature_list is None:
            #plot one hull at 0K first and here energy_key will become energy_key_per_atom
                print('Plotting hull using %s'%self.energy_key)
                QueryConvexHull(cursor=self.cursor,energy_key=self.energy_key,**kwargs)
            else:
                print('Plotting temperature dependent hull...')
                self.plot_td_hull(**kwargs)
        if self.subcmd == 'voltage':
            if self.temperature_list is None:
                print('Plotting voltage at 0K using %s'%self.energy_key)
                QueryConvexHull(cursor=self.cursor,energy_key=self.energy_key,**kwargs) 
            else:
                print('Plotting temperature dependent voltage curve...')
                self.plot_td_voltage(**kwargs)

    def plot_td_hull(self, **kwargs):
        """ Plot 2D Temperature Dependent hull at temperatures
        in temperature_list
        """
        import matplotlib
        import matplotlib.pyplot as plt
        import seaborn as sns
        from matador.plotting import plot_2d_hull
        import numpy as np
        #initialize plotting params
        matplotlib.rcParams['figure.figsize'] = (8,6)
        fig, ax = plt.subplots(figsize=(8,6))
        temperature_colours = plt.cm.get_cmap('brg')(np.linspace(0, 1, (len(self.temperature_list)+1)*2))
        hull = QueryConvexHull(no_plot=True,cursor=self.cursor,
            energy_key=self.energy_key,elements=self.elements,**kwargs)
        hull.set_plot_param() #this is a necessary line to initialize plotting
        current_color=temperature_colours[0]


        # plot static hull first (enthalpy)
        ax = plot_2d_hull(hull,ax=ax,show=False, plot_points=False,)
        self.plot_hull_line(hull=hull,ax=ax,label='Static Lattice',color='k')

        # get zero point and temperature dependent energies 
        self.get_zp_td_energy() 

        #plot 0K hull 
        hull = QueryConvexHull(no_plot=True,cursor=self.cursor,
            energy_key=self.energy_key,elements=self.elements,**kwargs)
        current_color=temperature_colours[0]
        self.plot_hull_line(hull=hull,ax=ax,label='0 K',color=current_color)

        #plot temperature dependent hulls
        for temp_ind, temperature1 in enumerate(sorted(self.temperature_list)):
            if not self.args.get('quiet') == True:
                print('Temperature is %f'%temperature1)
            current_color = temperature_colours[temp_ind+1]
            hull = QueryConvexHull(no_plot=True,cursor=self.cursor,
                energy_key=self.thermo_energy_key_pa,
                temperature=temperature1,elements=self.elements,**kwargs)
            self.plot_hull_line(hull=hull,ax=ax,label='{} K'.format(temperature1),color=current_color)

        #add nice plotting features at end
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)
        ax.set_ylim([-0.18,0.05])
        ax.set_yticks(np.arange(-0.18,0.05,0.05))
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)
        #ax.set_title('A' + r'$_\mathrm{x}$' + 'B' + r'$_\mathrm{1-x}$')
        #ax.set_xlabel(r'x in {}$_\mathrm{{x}}${}$_\mathrm{{1-x}}$'.format('A', 'B'))
        ax.set_yticklabels(['{:.2f}'.format(val) for val in ax.get_yticks()]) 
        ax.legend(loc='lower right')
        if self.args.get('savefig') == True:
            plt.savefig('%s%s_tdhull.png'%(self.elements[0],self.elements[1]), bbox_inches='tight',dpi=300)
        if not self.args.get('no_temp_plot'):
            plt.show()      
        return 
            
             
    def plot_td_voltage(self, **kwargs):
        """ Plot 2D Temperature Dependent hull at temperatures
        in temperature_list
        """
        import matplotlib
        import matplotlib.pyplot as plt
        import seaborn as sns
        import numpy as np
        
        #setup figure for plotting
        matplotlib.rcParams['figure.figsize'] = (8,6)
        fig, axQ = plt.subplots(figsize=(8,6))
        temperature_colours = plt.cm.get_cmap('brg')(np.linspace(0, 1, (len(self.temperature_list)+1)*2))
        
        #first plot static line
        hull = QueryConvexHull(no_plot=True,cursor=self.cursor,energy_key=self.energy_key,elements=self.elements,**kwargs)
        self.plot_voltage_line(hull,axQ=axQ,label='Static Lattice',color='k')

        #get zero point corrected energy 
        self.get_zp_td_energy()
 
        #plot 0K hull next with zero point correction
        hull = QueryConvexHull(no_plot=True,cursor=self.cursor,
            energy_key=self.energy_key,elements=self.elements,**kwargs)
        hull.set_plot_param() #this is a necessary line to initialize plotting
        current_color=temperature_colours[0]
        self.plot_voltage_line(hull,axQ=axQ,label='0 K',color=current_color)

        # plot temperature dependent curves
        for temp_ind, temperature1 in enumerate(sorted(self.temperature_list)):
            if not self.args.get('quiet'):
                print('Temperature is %f'%temperature1)
            current_color = temperature_colours[temp_ind+1]
            hull = QueryConvexHull(no_plot=True,cursor=self.cursor,
                energy_key=self.thermo_energy_key_pa,
                temperature=temperature1,elements=self.elements,**kwargs)
            self.plot_voltage_line(hull,axQ=axQ,label='{} K'.format(temperature1),color=current_color)
    
        #set axes labels
        axQ.set_ylabel('Voltage (V) vs {}$^+$/{}'.format(hull.elements[0], hull.elements[0]))
        axQ.set_xlabel('Gravimetric cap. (mAh/g)')
        _, end = axQ.get_ylim()
        axQ.set_ylim(0, 1.1 * end)
        _, end = axQ.get_xlim()
        axQ.set_xlim(0, 1.1 * end)
        axQ.grid('off')
        plt.tight_layout(pad=0.0, h_pad=1.0, w_pad=0.2)
        try:
            import seaborn as sns
            sns.despine()
            dark_grey = '#262626'
            for spine in ['left', 'bottom']:
                axQ.spines[spine].set_linewidth(0.5)
                axQ.spines[spine].set_color(dark_grey)
        except ImportError:
            pass
        for item in ([axQ.title, axQ.xaxis.label, axQ.yaxis.label] +
            axQ.get_xticklabels() + axQ.get_yticklabels()):
            item.set_fontsize(20)
        axQ.legend(loc='upper right')
        if self.args.get('savefig') == True:
            plt.savefig('%s%s_tdvoltage.png'%(self.elements[0],self.elements[1]), bbox_inches='tight',dpi=300)
        if not self.args.get('no_temp_plot'):
            plt.show()      
        return 

    def plot_hull_line(self,hull,ax,label='',color='k'):
        '''plot one hull curve on a given set of axes
        '''
        import numpy as np

        tie_line = hull.structure_slice[hull.hull.vertices]
        ax.scatter(tie_line[:,0],tie_line[:,1],
            c=color,
            marker='o', label=label,
            zorder=99999,
            edgecolor='k',
            lw=1,
            alpha=1)
        ax.scatter(hull.structures[np.argsort(hull.hull_dist), 0][::-1],
            hull.structures[np.argsort(hull.hull_dist), -1][::-1],
            lw=1, alpha=1,edgecolor='k',zorder=10000,
            c=color)
        ax.plot(np.sort(tie_line[:,0]), tie_line[np.argsort(tie_line[:,0]), 1],
            c=color, lw=2, alpha=1, zorder=1000)
        return

    def get_zp_td_energy(self): 
        '''get zero point energy and temperature dependent energy
        '''
        for doc in self.cursor:
            doc[self.energy_key] += doc['zero_point_E']
            doc[self.energy_key + '_per_atom'] += (doc['zero_point_E']  / doc['num_atoms'])
        #calculate thermo_free_energy_per_atom for each doc in cursor
        for doc in self.cursor:
            doc[self.thermo_energy_key_pa] = {}
            for temp in doc['thermo_temps']:
                #calculate the thermo_free_energy_per_atom
                doc[self.thermo_energy_key_pa][temp] = doc[self.energy_key + '_per_atom'] + (doc[self.thermo_energy_key][temp]) / doc['num_atoms']
                doc[self.thermo_energy_key][temp] += doc[self.energy_key]
        return

    def plot_voltage_line(self,hull,axQ,label='',color='k'):
        '''plot one voltage curve on a given set of axes'''
        for ind, voltage in enumerate(hull.voltages):
            for i in range(1, len(voltage) - 1):
                if i == 1 and hull.args.get('expt'):
                    axQ.plot([hull.Q[ind][i - 1], hull.Q[ind][i]], [voltage[i], voltage[i]],
                         marker='*',
                         lw=2,
                         c=color, 
                         label='DFT (this work)')
                elif i == 1 and len(hull.voltages) != 1:
                    axQ.plot([hull.Q[ind][i - 1], hull.Q[ind][i]], [voltage[i], voltage[i]],
                         marker='o',
                         markersize=5,
                         lw=2,
                         c=color,
                         label=get_formula_from_stoich(hull.endstoichs[ind], tex=True))
                elif i == 1:
                    axQ.plot([hull.Q[ind][i - 1], hull.Q[ind][i]], [voltage[i], voltage[i]],
                         marker='o',
                         markersize=5,
                         lw=2,
                         c=color, label=label)
                else:
                    axQ.plot([hull.Q[ind][i - 1], hull.Q[ind][i]], [voltage[i], voltage[i]],
                         marker='o',
                         markersize=5,
                         lw=2,
                         c=color)
                if i != len(voltage) - 2:
                    axQ.plot([hull.Q[ind][i], hull.Q[ind][i]], [voltage[i], voltage[i + 1]],
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
                             xy=(hull.Q[0][i+1]+0.02*max(hull.Q[0]),
                                 hull.voltages[0][i+1]+0.02*max(hull.voltages[0])),
                             textcoords='data',
                             ha='center',
                            zorder=9999)
            if hull.args.get('expt') or len(hull.voltages) != 1:
                axQ.legend(loc=1)
        return
