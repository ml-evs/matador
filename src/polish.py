#!/usr/bin/python
# coding: utf-8
""" This file implements the creation of 
new input files from a query and a desired
level of accuracy.
"""
from __future__ import print_function
from scrapers.castep_scrapers import cell2dict
from scrapers.castep_scrapers import param2dict
from query import query2files
from traceback import print_exc


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
        self.cursor = cursor
        # parse new parameters
        if self.args.get('to') is not None:
            self.template_structure = query.template_structure
        if self.args.get('with') is not None:
            self.template_seedname = self.args.get('with')
            try:
               cell_dict, cell_success = cell2dict(self.template_seedname) 
               param_dict, param_success = param2dict(self.template_seedname) 
               if not cell_sucess:
                   if not param_success:
                       raise RuntimeError('Failed to read cell and param file.')
                    raise RuntimeError('Failed to read cell file.')
               if not param_sucess:
                   raise RuntimeError('Failed to read param file.')
            except Exception:
                print_exc()
                exit()
        if self.args['subcmd'] == 'swaps':
            # to-do parse swap command line
            self.parse_swaps()
            swap_cursor = []
            for doc in self.cursor:
                swap_cursor.append(self.atomic_swaps(doc, self.args))
            self.cursor = swap_cursor
        polish_cursor = []
        for doc in self.cursor:
            polish_cursor.append(self.change_accuracy(doc, self.args))
        self.cursor = polish_cursor
        query2files(self.cursor)

    def change_accuracy(self):
        """ Augment a document to have the desired
        parameters for polishing.
        """
    
    def parse_swaps(self):
        """ Parse command line options into valid
        atomic species swaps.
        """

        
    def atomic_swaps(self):
        """ Swap atomic species according to parsed
        options.
        """
