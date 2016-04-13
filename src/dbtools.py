#!/usr/bin/python
from __future__ import print_function
import pymongo as pm
import argparse
import bson.json_util as json

def dbstats():

    client = pm.MongoClient()
    db = client.crystals
    repo = db.repo
    print(json.dumps(db.command('dbstats'),indent=2))
    for name in db.collection_names():
        print(json.dumps(db.command('collstats', name),indent=2))
    
if __name__ == '__main__':
    
    client = pm.MongoClient()
    db = client.crystals

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--drop', type=str,
            help='drop named collection')
    args = parser.parse_args() 
    if args.drop != None:
        print(args.drop)
        try:
            db[str(args.drop)].dropIndexes()
        except:
            pass
        db[str(args.drop)].drop()
        print('dropped', str(args.drop))
    dbstats()

