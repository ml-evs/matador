# encoding: utf-8
""" Some simple utilities for DB connections. """
from os import uname
import pymongo as pm


def make_connection_to_collection(coll_names):
    """ Connect to database of choice.

    Input:

        | coll_names: str, name of collection.

    Returns:

        | client      : MongoClient, the connection to the database
        | db          : Database, the database to query
        | collections : dict(Collection), collections stored by name

    """

    local = uname()[1]
    if local == 'cluster2':
        remote = 'node1'
    else:
        remote = None
    client = pm.MongoClient(remote)
    db = client.crystals
    collections = dict()
    if coll_names is not None:
        if not isinstance(coll_names, list):
            coll_names = [coll_names]
        for database in coll_names:
            if database == 'ajm':
                database = 'repo'
                collections['ajm'] = db['repo']
            else:
                collections[database] = db[database]
    else:
        collections['ajm'] = db['repo']

    return client, db, collections
