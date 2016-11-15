#!/usr/bin/python
# coding: utf-8
""" This file implements the creation of
new input files from a query and a desired
level of accuracy and atomic swaps.
"""

from __future__ import print_function

# matador functionality
from scrapers.castep_scrapers import cell2dict
from scrapers.castep_scrapers import param2dict
from print_utils import print_notify, print_success, print_warning, print_failure
from chem_utils import get_periodic_table
from export import query2files

# standard library
from traceback import print_exc
from copy import deepcopy
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
        self.periodic_table = get_periodic_table()
        del self.periodic_table['X']
        self.template_structure = None
        try:
            self.cursor = list(cursor)
            if self.args.get('top') is not None:
                self.cursor = self.cursor[:self.args.get('top')]
        except:
            print_exc()
            print_failure('Something went wrong!')
            exit()

        # parse new parameters
        self.cell_dict, self.param_dict = self.get_accuracy()
        if self.args['subcmd'] == 'swaps':
            self.swap_counter = 0
            self.parse_swaps()
            swap_cursor = []
            for doc in self.cursor[:]:
                docs, counter = self.atomic_swaps(doc)
                self.swap_counter += counter
                if counter > 0:
                    swap_cursor.extend(docs)
            self.cursor = swap_cursor
            if self.swap_counter > 0:
                print_success('Performed ' + str(self.swap_counter) + ' swaps.')
            else:
                print_warning('No swaps performed.')

        polish_cursor = []
        for doc in self.cursor[:]:
            polish_cursor.append(self.change_accuracy(doc))
        self.cursor = polish_cursor

        # query2files(self.cursor, self.args)

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

        # split by comma to get pairs of swaps
        if self.args.get('debug'):
            print(self.args.get('swap'))
        if len(self.args.get('swap')) > 1:
            print_failure('Detected whitespace in your input, ' +
                          'clear it and try again.')
            exit()
        swap_list = self.args.get('swap')[0].split(':')
        if self.args.get('debug'):
            print(swap_list)
        for swap in swap_list:
            if len(swap) <= 1:
                exit('Not enough arguments for swap!')
            tmp_list = re.split(r'([A-Z][a-z]*)', swap)
            if self.args.get('debug'):
                print(tmp_list)
            # scrub square brackets
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
            if self.args.get('debug'):
                print(tmp_list)

            # parse list of elements or group
            for ind, atom in enumerate(tmp_list):
                if '[' in atom:
                    group = atom.strip(']').strip('[')
                    if group == 'X':
                        print_failure('Cannot swap from macro [X], please reconsider...')
                        exit()
                    if group in self.periodic_table:
                        atoms = self.periodic_table[group]
                    else:
                        atoms = group.split(',')
                    swap_test[ind] = atoms
                else:
                    swap_test[ind] = [atom]
            self.swap_pairs.append(swap_test)
            if self.args.get('debug'):
                print(self.swap_pairs)
        for pair in self.swap_pairs:
            print_notify('Swapping all ' + str(pair[0]) + ' for ' + str(pair[1]))

    def atomic_swaps(self, source_doc):
        """ Swap atomic species according to parsed
        options.
        """
        doc = source_doc
        new_doc = deepcopy(doc)
        swapped_docs = []
        swapped = False
        # iterate over sets of swaps
        for swap_pair in self.swap_pairs:
            if self.args.get('debug'):
                print(swap_pair)
            # for each atom to be swapped
            for swap_atom in swap_pair[0]:
                # if structure contains an atom to be swapped
                if swap_atom in new_doc['atom_types']:
                    # iterate over new atoms and swap
                    for new_atom in swap_pair[1]:
                        # don't swap structure to itself
                        if new_atom != swap_atom:
                            # iterate over structure and swap all to new atom
                            for ind, source_atom in enumerate(new_doc['atom_types']):
                                if source_atom == swap_atom:
                                    new_doc['atom_types'][ind] = new_atom
                                    swapped = True
        if swapped:
            # add to list of all structures post-swapping
            swapped_docs.append(new_doc)
        return swapped_docs, len(swapped_docs)
