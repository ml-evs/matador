#!/usr/bin/python
from __future__ import print_function
import pymongo as pm
import numpy as np
import argparse
import bson.json_util as json

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
        print(args)
        if 'stoichiometry' in args:
            cursor = self.query_stoichiometry()
        self.display_results(cursor)
        # for doc in cursor:
            # print(json.dumps(doc, indent=2))

    def display_results(self, cursor):
        ''' Print query results in a cryan-like fashion. '''
        struct_string = []
        gs_enthalpy = 0
        header_string = 'ObjectId\t\t\tPressure\tVolume\t\t\tEnthalpy\t\tSpace group\t\tComposition'
        for ind, doc in enumerate(cursor):
            struct_string.append(
                      str(doc['_id']) +                           '\t' \
                    + str(doc['pressure']) +                      '\t\t' \
                    + str(doc['cell_volume']) +                   '\t\t' \
                    + str(doc['enthalpy_per_atom']-gs_enthalpy) + '\t\t' \
                    + doc['space_group'] +                        '\t\t' \
                    + str((doc['stoichiometry'])))
            if ind == 0:
                gs_enthalpy = doc['enthalpy_per_atom']
        print(header_string)
        for string in struct_string:
            print(string)
        
    def query_stoichiometry(self):
        ''' Query DB for particular stoichiometry. '''
        # alias stoichiometry
        stoich = self.args.stoichiometry
        elements = []
        fraction = []
        for i in range(0, len(stoich), 2):
            elements.append(stoich[i])
            fraction.append(float(stoich[i+1]))
        fraction = np.asarray(fraction)
        fraction /= np.min(fraction)
        print(zip(elements,fraction))
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
        print(cursor.count())
        cursor.sort('enthalpy_per_atom', pm.ASCENDING)
        return cursor[0:10]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--stoichiometry', nargs='+', type=str)
    args = parser.parse_args()
    
    query = DBQuery(stoichiometry=args.stoichiometry)



