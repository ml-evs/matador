# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements atomic swaps through the `AtomicSwapper` class. """


import re
from copy import deepcopy
from matador.utils.print_utils import print_success, print_warning, print_failure
from matador.utils.chem_utils import get_periodic_table, get_stoich, get_formula_from_stoich, get_root_source


class AtomicSwapper:
    """ This class handles the creation of input files from database
    queries that have swapped atoms.

    """
    def __init__(self, cursor, **kwargs):
        """ Initialise class with query cursor and arguments.

        Parameters:
            cursor (list): cursor of documents to swap.

        """
        self.args = kwargs
        # define some swap macros
        self.periodic_table = get_periodic_table()
        self.swap_dict_list = None
        del self.periodic_table['X']
        self.template_structure = None
        self.cursor = list(cursor)
        if self.args.get('top') is not None:
            self.cursor = self.cursor[:self.args.get('top')]

        # parse new parameters
        self.cell_dict, self.param_dict = {}, {}
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

        if self.args.get('uniq'):
            from matador.similarity.similarity import get_uniq_cursor
            print('Filtering for unique structures...')
            unique_set, _, _, _ = get_uniq_cursor(self.cursor,
                                                  debug=self.args.get('debug'),
                                                  sim_tol=self.args.get('uniq'))
            print('Filtered {} down to {}'.format(len(self.cursor), len(unique_set)))
            self.cursor = [self.cursor[ind] for ind in unique_set]

    def parse_swaps(self, swap_args=None):
        """ Parse command line options into valid atomic species swaps.

        e.g. --swap LiP:NaAs

            ==> [[['Li'], ['P']], [['Na'], ['P']].

        Handles multiple many-to-many swaps, macros for groups of the
        periodic table, and wildcards.

        Keyword arguments:
            swap_args (str): overrides command-line swap args.


        """

        self.swap_pairs = []
        if swap_args is None:
            swap_args = self.args.get('swap')
        if isinstance(swap_args, str):
            swap_args = [swap_args.strip()]
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
                tmp_list = [x for x in swap.split('][') if x != '']
            # check if only first option is group
            elif swap[0] == '[':
                tmp_list = [x for x in swap.split(']') if x != '']
            # check if only last option is group
            elif swap[-1] == ']':
                tmp_list = [x for x in swap.split('[') if x != '']
            # check if no groups
            else:
                tmp_list = [x for x in re.split(r'([A-Z][a-z]*)', swap) if x != '']
            for ind, tmp in enumerate(tmp_list):
                tmp_list[ind] = self._atoms_to_list(tmp)
            if len(tmp_list) != 2:
                raise RuntimeError('Unable to parse swap! {} should contain only two entries'.format(tmp_list))
            self.swap_pairs.append(tmp_list)
            self.construct_swap_options()

    def _atoms_to_list(self, atom_string):
        """ For a given set of atoms in a string, parse any macros and
        return a list of options.

        e.g. '[V' -> [<all group V atoms>],
        and 'V' -> ['V'].

        Parameters:
            atom_string (str): formula string with macros.

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
        """ Iterate over possible combinations of multiple many-to-many
        swaps and create a dict for each swap.

        """
        self.swap_dict_list = []
        from itertools import product
        for branch in product(*([pair[1] for pair in self.swap_pairs])):
            self.swap_dict_list.append(dict())
            for ind, pair in enumerate(self.swap_pairs):
                for swap_from in pair[0]:
                    if swap_from != branch[ind]:
                        self.swap_dict_list[-1][swap_from] = branch[ind]

    def atomic_swaps(self, source_doc):
        """ Swap atomic species according to parsed options.

        Parameters:
            source_doc (dict): matador doc to swap from.

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
                swapped_doc['stoichiometry'] = get_stoich(swapped_doc['atom_types'])
                if 'source' in swapped_doc:
                    del swapped_doc['source']
                swapped_doc['source'] = [get_formula_from_stoich(swapped_doc['stoichiometry']) + '-swap-' + get_root_source(doc)]
                swapped_docs.append(swapped_doc)
        return swapped_docs, len(swapped_docs)
