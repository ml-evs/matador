#!/usr/bin/python
# coding: utf-8
""" This file implements the creation of
new input files from a query and a desired
level of accuracy and atomic swaps.
"""

from __future__ import print_function
# matador modules
from .scrapers.castep_scrapers import cell2dict
from .scrapers.castep_scrapers import param2dict
from .utils.print_utils import print_success, print_warning, print_failure
from .utils.chem_utils import get_periodic_table
# standard library
from traceback import print_exc
from copy import deepcopy
import re


class Polisher(object):
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

    def parse_swaps(self, swap_args=None):
        """ Parse command line options into valid
        atomic species swaps.

        e.g. --swap LiP:NaAs

            ==> [[['Li'], ['P']], [['Na'], ['P']].

        Handles multiple many-to-many swaps, macros for
        groups of the periodic table, and wildcards.
        """

        self.swap_pairs = []
        if swap_args is None:
            swap_args = self.args.get('swap')
        if len(swap_args) > 1:
            print_failure('Detected whitespace in your input, ' +
                          'clear it and try again.')
            exit()
        swap_list = swap_args[0].split(':')
        for swap in swap_list:
            if len(swap) <= 1:
                exit('Not enough arguments for swap!')
            # check is both options are groups
            if '][' in swap:
                tmp_list = [x for x in swap.split('][') if x is not '']
            # check if only first option is group
            elif swap[0] is '[':
                tmp_list = [x for x in swap.split(']') if x is not '']
            # check if only last option is group
            elif swap[-1] is ']':
                tmp_list = [x for x in swap.split('[') if x is not '']
            # check if no groups
            else:
                tmp_list = [x for x in re.split(r'([A-Z][a-z]*)', swap) if x is not '']
            for ind, tmp in enumerate(tmp_list):
                tmp_list[ind] = self._atoms_to_list(tmp)
            assert(len(tmp_list) == 2)
            self.swap_pairs.append(tmp_list)
            self.construct_swap_options()

    def _atoms_to_list(self, atom_string):
        """ For a given set of atoms in a string,
        parse any macros and return a list of options.

        e.g. '[V' -> [<all group V atoms>],
        and 'V' -> ['V'].
        """
        if '[' in atom_string or ']' in atom_string:
            group = atom_string.replace('[', '')
            group = group.replace(']', '')
            if group in self.periodic_table:
                atom_list = self.periodic_table[group]
            else:
                atom_list = group.split(',')
        else:
            return [atom_string]
        return [x.strip() for x in atom_list]

    def construct_swap_options(self):
        """ Iterate over possible combinations of multiple
        many-to-many swaps and create a dict for each swap.
        """
        from itertools import product
        self.swap_dict_list = []
        for branch in product(*([pair[1] for pair in self.swap_pairs])):
            self.swap_dict_list.append(dict())
            for ind, pair in enumerate(self.swap_pairs):
                for swap_from in pair[0]:
                    if swap_from != branch[ind]:
                        self.swap_dict_list[-1][swap_from] = branch[ind]

    def atomic_swaps(self, source_doc):
        """ Swap atomic species according to parsed
        options.
        """
        doc = source_doc
        new_doc = deepcopy(doc)
        swapped_docs = []
        swapped = False
        for swap in self.swap_dict_list:
            for ind, source_atom in enumerate(doc['atom_types']):
                if source_atom in swap:
                    new_doc['atom_types'][ind] = swap[source_atom]
                    swapped = True
            if swapped:
                swapped_doc = deepcopy(new_doc)
                swapped_docs.append(swapped_doc)
        return swapped_docs, len(swapped_docs)
