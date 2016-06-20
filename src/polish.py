#!/usr/bin/python
# coding: utf-8
""" This file implements the creation of
new input files from a query and a desired
level of accuracy.
"""
from __future__ import print_function
from scrapers.castep_scrapers import cell2dict
from scrapers.castep_scrapers import param2dict
from export import query2files
from traceback import print_exc


class Polisher:
    """ This class handles the creation of
    input files from database queries that have
    a new level of accuracy.
    """
    def __init__(self, cursor, template_structure=None, *args):
        """ Initialise class with query cursor
        and arguments.
        """
        self.args = args[0]
        self.template_structure = template_structure
        self.cursor = cursor
        # parse new parameters
        self.cell, self.param = self.get_accuracy()
        if self.args['subcmd'] == 'swaps':
            # to-do parse swap command line
            self.parse_swaps()
            swap_cursor = []
            for doc in self.cursor:
                swap_cursor.append(self.atomic_swaps(doc))
            self.cursor = swap_cursor
        polish_cursor = []
        for doc in self.cursor:
            polish_cursor.append(self.change_accuracy(doc))
        self.cursor = polish_cursor
        query2files(self.cursor)

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
        if self.args.get('with') is not None:
            self.template_seedname = self.args.get('with')
            try:
                cell_dict, cell_success = cell2dict(self.template_seedname)
                param_dict, param_success = param2dict(self.template_seedname)
                if not cell_success:
                    if not param_success:
                        raise RuntimeError('Failed to read cell and param file.')
                    raise RuntimeError('Failed to read cell file.')
                if not param_success:
                    raise RuntimeError('Failed to read param file.')
            except Exception:
                print_exc()
                exit()
        elif self.args.get('to') is not None:
            doc = self.template_structure
            self.param_dict = dict()
            param_list = ['cut_off_energy', 'xc_functional',
                          'spin_polarized']
            cell_list = ['kpoints_mp_spacing', 'species_pot',
                         'external_pressure']
            for param in [param for param in param_list if param in doc]:
                param_dict[param] = doc[param]
            for cell in [cell for cell in cell_list if cell in doc]:
                cell_dict[cell] = doc[cell]
        final_param = param_dict.update(default_parameters)
        final_cell = final_cell.update(cell_dict)
        return final_cell, final_param

    def change_accuracy(self, doc):
        """ Augment a document to have the desired
        parameters for polishing.
        """
        doc = doc.update(self.cell)
        doc = doc.update(self.param)
        return doc

    def parse_swaps(self):
        """ Parse command line options into valid
        atomic species swaps.
        e.g. --swap Li,As |--> ['Li', 'As']
        """
        self.swap_pairs = []
        # read in pairs of atoms, check they are valid

    def atomic_swaps(self, doc):
        """ Swap atomic species according to parsed
        options.
        """
        for ind, atom in enumerate(doc['atom_types']):
            for swap_pair in self.swap_pairs:
                if atom == swap_pair[0]:
                    doc['atom_types'][ind] = swap_pair[1]
        return doc
