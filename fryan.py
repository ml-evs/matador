#!/usr/bin/python
# coding: utf-8
from __future__ import print_function
import pymongo as pm
import numpy as np
import argparse
import bson.json_util as json
import re

class DBQuery:
    ''' Class that implements queries to MongoDB
    structure database.
    '''

    def __init__(self, **kwargs):
        ''' Initialise the query with command line
        arguments.
        '''
        
        client = pm.MongoClient()
        self.repo = client.crystals.repo
        self.args = args
        self.top = args.top if args.top != None else 10
        self.details = args.details
        if args.stoichiometry != None:
            cursor = self.query_stoichiometry()
        if args.id != None:
            cursor = self.repo.find({'text_id': args.id})
        elif args.composition != None:
            cursor = self.query_composition()
        print(cursor.count(), 'structures found.')
        if cursor.count() != 0:
            self.display_results(cursor)
        # for doc in cursor:
            # print(json.dumps(doc, indent=2))

    def display_results(self, cursor):
        ''' Print query results in a cryan-like fashion. '''
        struct_string = []
        detail_string = []
        gs_enthalpy = 0
        header_string = "{:<20}".format('   Id')
        # header_string +=  "{:^24}".format('ObjectId')
        header_string += "{:^12}".format('Pressure')
        header_string += "{:^12}".format('Volume')
        header_string += "{:^18}".format('Enthalpy')
        header_string += "{:^12}".format('Space group')
        header_string += "{:^12}".format('Formula')
        header_string += "{:^10}".format('# f.u.')
        for ind, doc in enumerate(cursor):
            struct_string.append(
                    "{:<20}".format(doc['text_id'][0]+' '+doc['text_id'][1])
                    + "{:^ 12.3f}".format(doc['pressure'])
                    + "{:^12.3f}".format(doc['cell_volume']/doc['num_atoms'])
                    + "{:^18.5f}".format(doc['enthalpy_per_atom']-gs_enthalpy)
                    + "{:^12}".format(doc['space_group']))
            if self.details:
                detail_string.append(18 * ' ' + u"└──── "
                        + doc['xc_functional'] 
                        + ', ' + 'cutoff: '+ "{:4.2f}".format(doc['cut_off_energy']) + ' eV'
                        + ', ' + 'external pressure: ' + "{:4.2f}".format(doc['external_pressure'][0][0]) + ' GPa'
                        )
                try:
                    detail_string[-1] += ', ' + doc['kpoints_mp_spacing']
                except: 
                    pass
            sub_string = ''
            atom_per_fu = 0
            for item in doc['stoichiometry']:
                for item_ind, subitem in enumerate(item):
                    if item_ind == 0:
                        atom_per_fu += 1
                    sub_string += str(subitem)
            struct_string[-1] += "{:^12}".format(sub_string)
            struct_string[-1] += "{:^10}".format(doc['num_atoms']/atom_per_fu)
            if ind == 0:
                gs_enthalpy = doc['enthalpy_per_atom']
        import os
        print(len(header_string)*'─')
        print(header_string)
        print(len(header_string)*'─')
        for ind, string in enumerate(struct_string):
            print(string)
            if self.details:
                print(detail_string[ind])
        
    def query_stoichiometry(self):
        ''' Query DB for particular stoichiometry. '''
        # alias stoichiometry
        stoich = self.args.stoichiometry
        # if there's only one string, try split it by caps
        if len(stoich) == 1:
            stoich = [elem for elem in re.split(r'([A-Z][a-z]*)', stoich[0]) if elem]
        elements = []
        fraction = []
        for i in range(0, len(stoich), 1):
            if not bool(re.search(r'\d', stoich[i])):
                elements.append(stoich[i])
                try:
                    fraction.append(float(stoich[i+1]))
                except:
                    fraction.append(1.0)
        fraction = np.asarray(fraction)
        fraction /= np.min(fraction)
        # pyMongo doesn't like generators... could patch pyMongo?
        # cursor = self.repo.find({'stoichiometry.'+[element for element in elements]: {'$exists' : True}})
        if len(elements) == 1:
            cursor = self.repo.find({'stoichiometry' : {'$in' : [[elements[0], fraction[0]]]}})
        elif len(elements) == 2:
            cursor = self.repo.find({ '$and': [ 
                                        {'stoichiometry' : {'$in' : [[elements[0], fraction[0]]]}},
                                        {'stoichiometry' : {'$in' : [[elements[1], fraction[1]]]}}
                                    ]})
        elif len(elements) == 3:
            cursor = self.repo.find({ '$and': [ 
                                        {'stoichiometry' : {'$in' : [[elements[0], fraction[0]]]}},
                                        {'stoichiometry' : {'$in' : [[elements[1], fraction[1]]]}},
                                        {'stoichiometry' : {'$in' : [[elements[2], fraction[2]]]}}
                                    ]})
        elif len(elements) == 4:
            cursor = self.repo.find({ '$and': [ 
                                        {'stoichiometry' : {'$in' : [[elements[0], fraction[0]]]}},
                                        {'stoichiometry' : {'$in' : [[elements[1], fraction[1]]]}},
                                        {'stoichiometry' : {'$in' : [[elements[2], fraction[2]]]}},
                                        {'stoichiometry' : {'$in' : [[elements[3], fraction[3]]]}}
                                    ]})
        cursor.sort('enthalpy_per_atom', pm.ASCENDING)
        if cursor.count() < self.top:
            return cursor
        else:
            return cursor[:self.top]
    
    def query_composition(self):
        ''' Query DB for all structures containing 
        all the elements taken as input.
        '''
        elements = self.args.composition
        # if there's only one string, try split it by caps
        if len(elements) == 1:
            elements = [elem for elem in re.split(r'([A-Z][a-z]*)', elements[0]) if elem]
        try:
            for elem in elements:
                if bool(re.search(r'\d', elem)):
                    raise RuntimeError('Composition string cannot contain a number.')
        except Exception as oops:
            print(oops)
            return EmptyCursor()
        # pyMongo doesn't like generators... could patch pyMongo?
        # cursor = self.repo.find({'stoichiometry.'+[element for element in elements]: {'$exists' : True}})
        if len(elements) == 1:
            cursor = self.repo.find({'atom_types' : {'$in' : [elements[0]]}})
        elif len(elements) == 2:
            cursor = self.repo.find({ '$and': [ 
                                        {'atom_types' : {'$in' : [elements[0]]}},
                                        {'atom_types' : {'$in' : [elements[1]]}}
                                    ]})
        elif len(elements) == 3:
            cursor = self.repo.find({ '$and': [ 
                                        {'atom_types' : {'$in' : [elements[0]]}},
                                        {'atom_types' : {'$in' : [elements[1]]}},
                                        {'atom_types' : {'$in' : [elements[2]]}}
                                    ]})
        elif len(elements) == 4:
            cursor = self.repo.find({ '$and': [ 
                                        {'atom_types' : {'$in' : [elements[0]]}},
                                        {'atom_types' : {'$in' : [elements[1]]}},
                                        {'atom_types' : {'$in' : [elements[2]]}},
                                        {'atom_types' : {'$in' : [elements[3]]}}
                                    ]})
        cursor.sort('enthalpy_per_atom', pm.ASCENDING)
        if cursor.count() < self.top:
            return cursor
        else:
            return cursor[:self.top]

class EmptyCursor:
    ''' Empty cursor class for failures. '''
    def count(self):
        return 0 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-s', '--stoichiometry', nargs='+', type=str,
        help='choose a stoichiometry, e.g. Ge 1 Te 1 Si 3, or GeTeSi3')
    group.add_argument('-c', '--composition', nargs='+', type=str,
        help='find all structures containing the given elements, e.g. GeTeSi.')
    group.add_argument('-i', '--id', type=str, nargs='+',
            help='specify a particular structure from its text_id')
    parser.add_argument('-t', '--top', type=int,
            help='number of structures to show (DEFAULT: 10)')
    parser.add_argument('-d', '--details', action='store_true',
            help='show as much detail about calculation as possible')
    parser.add_argument('-p', '--pressure', type=float,
            help='specify an isotropic external pressure to search for, e.g. 10 (GPa)')
    args = parser.parse_args()
    query = DBQuery(stoichiometry=args.stoichiometry,
                    composition=args.composition,
                    id=args.id,
                    top=args.top,
                    details=args.details,
                    pressure=args.pressure)
