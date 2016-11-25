# coding: utf-8
""" This file implements the creation of
new input files from a query and a desired
level of accuracy and atomic swaps.
"""
from __future__ import print_function

import spglib as spg
from utils.chem_utils import get_atomic_number


class Refiner:
    """ A class to crudely assess the quality of
    symmetry info of a structure in the database, with spglib.
    """

    def __init__(self, cursor):
        wrong_count = dict()
        try:
            try:
                total = len(cursor)
            except:
                total = cursor.count()
            tested_count = 0
            failed_count = 0
            symprec = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
            print('Checking space groups...')
            for sym in symprec:
                wrong_count[sym] = 0
            for ind, doc in enumerate(cursor):
                try:
                    tested_count += 1
                    cell = (doc['lattice_cart'],
                            doc['positions_frac'],
                            [get_atomic_number(elem) for elem in doc['atom_types']])
                    for sym in symprec:
                        sg = spg.get_spacegroup(cell, symprec=sym)
                        if sg.split()[0] != doc['space_group']:
                            print((sym, ':', sg + ' ' + doc['space_group']), end='\t')
                            wrong_count[sym] += 1
                            # print((sym, ':', sg + ' ' + doc['space_group']), end='\t')
                    print('\n')
                except:
                    failed_count += 1
                    pass
                print('# tested:', tested_count, '# total:', total, end='\r')
            print('\n')
            for sym in symprec:
                print('Precision: ', sym, '\t# wrong:', wrong_count[sym])
            print('\n')
            print('# tested:', tested_count)
            print('# failed:', failed_count)
            print('# total:', total)
        except:
            for sym in symprec:
                print('Precision: ', sym, '# wrong:', wrong_count[sym])
            print('# tested:', tested_count)
        try:
            cursor.close()
        except:
            pass
