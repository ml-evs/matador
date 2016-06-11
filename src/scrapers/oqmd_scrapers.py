#!/usr/bin/python
# coding: utf-8
""" This file scrapes structures from the OQMD SQL dump.
You probably don't want to be using this unless you have
MySQL set up exactly as I do.
Matthew Evans 2016
"""
from __future__ import print_function
# external libraries
import MySQLdb
import MySQLdb.cursors
# import bson.json_util as json
# standard library
from collections import defaultdict
# from os.path import isfile


class OQMDConverter:
    """ The OQMDConverter class implements methods to scrape
    a MySQL table of OQMD structures which can be found at:
    http://oqmd.org/static/docs/getting_started.html.
    """
    def __init__(self):
        """ Connect to the database, initialise the query,
        and begin scraping.
        """
        self.db = MySQLdb.connect(host='localhost',
                                  user='root',
                                  cursorclass=MySQLdb.cursors.DictCursor,
                                  db='oqmd')
        # start by scraping all structures marked with 'fine_relax'
        cursor = self.db.cursor()
        cursor.execute("select * from calculations where label in ('fine_relax')")
        while True:
            doc = cursor.fetchone()
            if doc is None:
                break
            elif doc['converged'] == 0:
                break
            else:
                calculation_dict, success = self.oqmd_calculation2dict(doc)
            if success:
                entry_id = calculation_dict['entry_id']
                execute_string = "select * from structures where label in ('fine_relax') \
                                  and entry_id in ('" + entry_id + "')"
                cursor.execute(execute_string)
            else:
                continue

    def oqmd_calculation2dict(self, doc):
        """ Take a calculation from oqmd.calculations and
        scrape its settings, returning the ID to its structure
        if possible.
        """
        calculation = defaultdict(list)
        # try to get output ID first
        calculation['source'].append(dict(('entry_id', doc['entry_id'])))
        calculation['source'].append(dict(('input_id', doc['input_id'])))
        calculation['output_id'] = doc['output_id']
        calculation['composition'] = doc['composition_id'].strip()
        # grab energies
        calculation['energy'] = float(doc['energy'])
        calculation['energy_per_atom'] = float(doc['energy_pa'])
        # grab settings
        calculation['cut_off_energy'] = float(doc['encut'])
        calculation['xc_functional'] = doc['xc'].upper()
        calculation['elec_energy_tol'] = float(doc['ediff'])
        if doc['ispin'] == 2:
            calculation['spin_polarized'] = True
        calculation['species_pot'] = doc['potentials']
        return calculation, True



        

