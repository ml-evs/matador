#!/usr/bin/python
from __future__ import print_function
import pymongo as pm
import bson.json_util as json

def dbstats():

    client = pm.MongoClient()
    db = client.crystals
    repo = db.repo
    print(json.dumps(db.command('dbstats'),indent=2))
    print(json.dumps(db.command('collstats', 'repo'),indent=2))
    
if __name__ == '__main__':

    dbstats()

