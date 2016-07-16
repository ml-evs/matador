#!/usr/bin/python
# coding: utf-8
""" This file implements the creation of
new input files from a query and a desired
level of accuracy and atomic swaps.
"""
from __future__ import print_function
from scrapers.castep_scrapers import cell2dict
from scrapers.castep_scrapers import param2dict
from print_utils import print_notify, print_success, print_warning
from export import query2files
from traceback import print_exc
import re


class Polisher:
    """ This class handles the creation of
    input files from database queries that have
    a new level of accuracy.
    """
    def __init__(self, cursor, *args):
        """ Initialise class with query cursor
        and arguments.
        """
        self.args = args[0]
        # define some swap macros
        self.periodic_table = dict()
        self.periodic_table['I'] = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
        self.periodic_table['II'] = ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']
        self.periodic_table['III'] = ['B', 'Al', 'Ga', 'In', 'Tl']
        self.periodic_table['IV'] = ['C', 'Si', 'Ge', 'Sn', 'Pb']
        self.periodic_table['V'] = ['N', 'P', 'As', 'Sb', 'Bi']
        self.periodic_table['VI'] = ['O', 'S', 'Se', 'Te', 'Po']
        self.periodic_table['VII'] = ['F', 'Cl', 'Br', 'I', 'At']
        self.periodic_table['Tran'] = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                                       'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                                       'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
        self.periodic_table['Lan'] = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                                      'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
        self.periodic_table['Act'] = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
                                      'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
        self.template_structure = None
        try:
            self.cursor = cursor.clone()
        except:
            self.cursor = list(cursor)
        # parse new parameters
        self.cell_dict, self.param_dict = self.get_accuracy()
        if self.args['subcmd'] == 'swaps':
            self.swap_counter = 0
            self.parse_swaps()
            swap_cursor = []
            for doc in self.cursor[:]:
                doc, counter = self.atomic_swaps(doc)
                self.swap_counter += counter
                if counter == 1:
                    swap_cursor.append(doc)
            self.cursor = swap_cursor
            if self.swap_counter > 0:
                print_success('Performed swaps on ' + str(self.swap_counter) + ' structures.')
            else:
                print_warning('No swaps performed.')
        polish_cursor = []
        for doc in self.cursor[:]:
            polish_cursor.append(self.change_accuracy(doc))
        self.cursor = polish_cursor
        self.args['cell'] = False
        self.args['param'] = False
        self.args['res'] = True
        query2files(self.cursor, self.args)

    def get_accuracy(self):
        """ Read the correct key-value pairs
        either from the template file or template
        structure, and set some useful defaults.
        """
        final_param = dict()
        final_cell = dict()
        # set some decent defaults
        default_parameters = dict()
        default_parameters['task'] = 'GeometryOptimization'
        # to maintain kpoints_mp_spacing
        default_parameters['geom_max_iter'] = 3
        default_parameters['finite_basis_corr'] = 0
        default_parameters['write_bib'] = False
        default_parameters['geom_method'] = 'lbfgs'
        default_parameters['page_wvfns'] = 0
        default_parameters['write_cell_structure'] = True
        default_parameters['calculate_stress'] = True
        default_parameters['calculate_stress'] = True
        if self.args.get('with') is not None:
            self.template_seedname = self.args.get('with')
            try:
                cell_dict, success = cell2dict(self.template_seedname)
                if not success:
                    raise RuntimeError('Failed to read cell file.')
                param_dict, success = param2dict(self.template_seedname)
                if not success:
                    raise RuntimeError('Failed to read param file.')
            except Exception:
                print_exc()
                exit()
            default_parameters.update(param_dict)
            final_cell = cell_dict
            final_param = default_parameters
        elif self.args.get('to') is not None:
            doc = self.template_structure
            param_list = ['cut_off_energy', 'xc_functional',
                          'spin_polarized']
            cell_list = ['kpoints_mp_spacing', 'species_pot',
                         'external_pressure']
            for param in [param for param in param_list if param in doc]:
                param_dict[param] = doc[param]
            for cell in [cell for cell in cell_list if cell in doc]:
                cell_dict[cell] = doc[cell]
            default_parameters.update(param_dict)
            final_cell = cell_dict
            final_param = default_parameters
        else:
            final_param = default_parameters
            final_cell = dict()
        # scrub sources to prevent overwriting
        try:
            del final_cell['source']
            del final_param['source']
        except:
            pass
        return final_cell, final_param

    def change_accuracy(self, doc):
        """ Augment a document to have the desired
        parameters for polishing.
        """
        doc.update(self.cell_dict)
        doc.update(self.param_dict)
        return doc

    def parse_swaps(self):
        """ Parse command line options into valid
        atomic species swaps.
        e.g. --swap Li,As |--> ['Li', 'As']
        """
        self.swap_pairs = []
        swap_list = self.args.get('swap')[0].split(',')
        for swap in swap_list:
            if len(swap) <= 1:
                exit('Not enough arguments for swap!')
            tmp_list = re.split(r'([A-Z][a-z]*)', swap)
            for ind, tmp in enumerate(tmp_list):
                if tmp == '[':
                    while tmp_list[ind+1] != ']':
                        tmp_list[ind] += tmp_list[ind+1]
                        del tmp_list[ind+1]
                    tmp_list[ind] += ']'
            while ']' in tmp_list:
                tmp_list.remove(']')
            while '' in tmp_list:
                tmp_list.remove('')
            swap_test = tmp_list
            for ind, atom in enumerate(tmp_list):
                if '[' in atom:
                    group = atom.strip(']').strip('[')
                    atoms = self.periodic_table[group]
                    swap_test[ind] = atoms
                else:
                    swap_test[ind] = [atom]
            for old_atom in swap_test[0]:
                for new_atom in swap_test[1]:
                    if old_atom != new_atom:
                        self.swap_pairs.append([old_atom, new_atom])
        for pair in self.swap_pairs:
            print_notify('Swapping all ' + pair[0] + ' for ' + pair[1])

    def atomic_swaps(self, doc):
        """ Swap atomic species according to parsed
        options.
        """
        swapped = False
        for ind, atom in enumerate(doc['atom_types']):
            for swap_pair in self.swap_pairs:
                if atom == swap_pair[0]:
                    doc['atom_types'][ind] = swap_pair[1]
                    swapped = True
        if swapped:
            return doc, 1
        else:
            return doc, 0
