#!/usr/bin/python
# coding: utf-8
""" This file scrapes structures from the OQMD SQL dump.
You probably don't want to be using this unless you have
MySQL set up exactly as I do.
Matthew Evans 2016

TO-DO
* fix pspots
* find k-points
"""
from __future__ import print_function
# external libraries
import MySQLdb
import MySQLdb.cursors
import pymongo as pm
# standard library
from collections import defaultdict
from random import randint
from os.path import dirname, realpath
import argparse
import re
from traceback import print_exc
from ast import literal_eval


class OQMDConverter:
    """ The OQMDConverter class implements methods to scrape
    a MySQL table of OQMD structures which can be found at:
    http://oqmd.org/static/docs/getting_started.html.
    """
    def __init__(self, dryrun=False, debug=False, verbosity=0, scratch=False, label=None):
        """ Connect to the relevant databases and
        set off the scraper.
        """
        self.import_count = 0
        self.dryrun = dryrun
        self.debug = debug
        self.verbosity = verbosity
        self.scratch = scratch
        self.label = label
        # set up I/O for text_id
        if not self.dryrun:
            try:
                wordfile = open(dirname(realpath(__file__)) + '/words', 'r')
                nounfile = open(dirname(realpath(__file__)) + '/nouns', 'r')
                self.wlines = wordfile.readlines()
                self.num_words = len(self.wlines)
                self.nlines = nounfile.readlines()
                self.num_nouns = len(self.nlines)
                wordfile.close()
                nounfile.close()
            except Exception as oops:
                exit(oops)
        # connect to SQL database as root, with "use oqmd"
        print('Connecting to OQMD SQL...')
        self.oqmd = MySQLdb.connect(host='localhost',
                                    user='root',
                                    cursorclass=MySQLdb.cursors.DictCursor,
                                    db='oqmd')
        print('Successfully connected.')
        # if not dryrunning, connect to OQMD or scratch MongoDB
        if not self.dryrun:
            self.client = pm.MongoClient()
            self.db = self.client.crystals
            if self.scratch:
                self.repo = self.db.scratch
            else:
                self.repo = self.db.oqmd
        # scrape structures from SQL database
        self.sql2db()

    def oqmd_struct2db(self, struct):
        """ Insert completed Python dictionary into chosen
        database, with generated text_id.
        """
        plain_text_id = [self.wlines[randint(0, self.num_words-1)].strip(),
                         self.nlines[randint(0, self.num_nouns-1)].strip()]
        struct['text_id'] = plain_text_id
        struct_id = self.repo.insert_one(struct).inserted_id
        if self.debug:
            print('Inserted', struct_id)
        return 1

    def sql2db(self):
        """ Perform SQL query to scrape, and hop around the tables
        to collect all data associated with one structure.
        """
        # start by scraping all converged structures with chosen label
        cursor = self.oqmd.cursor()
        cursor.execute("select count(label) from calculations where label in ('" +
                       self.label + "') and converged in ('1')")
        count = cursor.fetchone()['count(label)']
        if count == 0:
            exit('No structures found with that label.')
        print(count, 'converged structures found.')
        # select only converged structures
        cursor.execute("select * from calculations where label in ('" +
                       self.label + "') and converged in ('1')")
        success_count = 0
        for row in cursor:
            calc_doc = row
            if calc_doc is None:
                continue
            calculation_dict, success = self.oqmd_calculation2dict(calc_doc)
            if not success:
                continue
            entry_id = calculation_dict['source']['entry_id']
            sql_query = "select * from structures where label in ('" + \
                        self.label + "') and entry_id in ('" + str(entry_id) + "')"
            cursor.execute(sql_query)
            # should only be one matching; revise at later date
            struct_doc = cursor.fetchone()
            structure_dict, success = self.oqmd_structure2dict(struct_doc)
            if not success:
                continue
            structure_id = structure_dict['source']['structure_id']
            # grab spacegroup symbol from ID
            spacegroup_id = structure_dict['spacegroup_id']
            sql_query = "select * from spacegroups where number in \
                              ('" + str(spacegroup_id) + "')"
            cursor.execute(sql_query)
            try:
                space_group = cursor.fetchone()
                structure_dict['space_group'] = space_group['hm']
            except:
                structure_dict['space_group'] = 'xxx'
                if self.debug:
                    print('Failed to get space group.')
            sql_query = "select * from atoms where structure_id in \
                              ('" + str(structure_id) + "')"
            cursor.execute(sql_query)
            atom_docs = cursor.fetchall()
            if len(atom_docs) != structure_dict['num_atoms']:
                print('Incorrect number of atoms! Skipping!')
                continue
            atoms_dict, success = self.oqmd_atoms2dict(atom_docs)
            if not success:
                continue
            final_struct = calculation_dict.copy()
            final_struct.update(structure_dict)
            final_struct.update(atoms_dict)
            success_count += 1
            if not self.dryrun:
                self.import_count += self.oqmd_struct2db(final_struct)
        if self.dryrun:
            print('Successfully scraped', success_count, '/',
                  count, 'structures.')
        if not self.dryrun:
            print('Successfully imported', self.import_count, '/',
                  count, 'structures.')
        return

    def oqmd_calculation2dict(self, doc):
        """ Take a calculation from oqmd.calculations and
        scrape its settings, returning the ID to its structure
        if possible.
        """
        try:
            calculation = defaultdict(list)
            # try to get output ID first
            calculation['source'] = dict()
            calculation['source']['entry_id'] = doc['entry_id']
            calculation['source']['input_id'] = doc['input_id']
            # convert stoich from string to list to tuple
            temp_stoich_list = [elem for elem in
                                re.split(r'([A-Z][a-z]*)',
                                         doc['composition_id']) if elem]
            for ind, item in enumerate(temp_stoich_list):
                if ind % 2 == 0:
                    calculation['stoichiometry'].append([str(temp_stoich_list[ind]),
                                                         int(temp_stoich_list[ind+1])])
            # grab energies and pretend they are enthalpies
            calculation['enthalpy'] = float(doc['energy'])
            calculation['enthalpy_per_atom'] = float(doc['energy_pa'])
            # grab settings
            settings = literal_eval(doc['settings'])
            calculation['cut_off_energy'] = settings['encut']
            calculation['xc_functional'] = settings['potentials'][0]['xc'].upper()
            calculation['elec_energy_tol'] = settings['ediff']
            if settings['ispin'] == 2:
                calculation['spin_polarized'] = True
            # calculation['species_pot'] = settings['potentials']
        except Exception:
            if self.debug:
                print('Scraping calculation failed.')
                print_exc()
                print(80*'─')
            return dict(), False
        return calculation, True

    def oqmd_structure2dict(self, doc):
        """ Create a dict containing structural information
        from oqmd.structures.
        """
        try:
            structure = defaultdict(list)
            # append ID to source
            structure['source'] = dict()
            structure['source']['structure_id'] = doc['id']
            structure['spacegroup_id'] = doc['spacegroup_id']
            structure['num_atoms'] = doc['natoms']
            # get pressure from -1/3 Tr(stress)
            structure['pressure'] = -(doc['sxx'] + doc['syy'] + doc['szz'])/3.0
            structure['lattice_cart'].append([doc['x1'], doc['x2'], doc['x3']])
            structure['lattice_cart'].append([doc['y1'], doc['y2'], doc['y3']])
            structure['lattice_cart'].append([doc['z1'], doc['z2'], doc['z3']])
            structure['cell_volume'] = doc['volume']
        except Exception:
            if self.debug:
                print('Scraping structure failed.')
                print_exc()
                print(80*'─')
            return dict(), False
        return structure, True

    def oqmd_atoms2dict(self, atom_docs):
        """ Take list of atom documents matching the
        structure_id and scrape them into a dict.
        """
        try:
            atoms = defaultdict(list)
            max_force = 0
            for doc in atom_docs:
                atoms['atom_types'].append(doc['element_id'])
                atoms['positions_frac'].append([doc['x'], doc['y'], doc['z']])
                force_tmp = doc['fx']**2 + doc['fy']**2 + doc['fz']**2
                if force_tmp > max_force:
                    max_force = force_tmp
                atoms['spins'].append(doc['magmom'])
                atoms['charges'].append(doc['charge'])
            atoms['max_force_on_atom'] = max_force
        except Exception:
            if self.debug:
                print('Scraping atoms failed.')
                print_exc()
                print(80*'─')
            return dict(), False
        return atoms, True

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Import OQMD (http://oqmd.org) structures into MongoDB database.',
            epilog='Written by Matthew Evans (2016)')
    parser.add_argument('-d', '--dryrun', action='store_true',
                        help='run the importer without connecting to the database')
    parser.add_argument('-v', '--verbosity', action='count',
                        help='enable verbose output')
    parser.add_argument('-l', '--label', type=str,
                        help='choose which OQMD calculation label to query, best options \
                        are fine_relax or relaxation.')
    parser.add_argument('--debug', action='store_true',
                        help='enable debug output to print every dict')
    parser.add_argument('-s', '--scratch', action='store_true',
                        help='import to junk collection called scratch')
    args = parser.parse_args()
    if args.label is None:
        exit('Please choose a label to scrape.')
    importer = OQMDConverter(dryrun=args.dryrun,
                             debug=args.debug,
                             verbosity=args.verbosity,
                             scratch=args.scratch,
                             label=args.label)
