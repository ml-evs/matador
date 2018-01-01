# encoding: utf-8
""" Some simple utilities for DB connections. """
from os import uname
import pymongo as pm


def make_connection_to_collection(coll_names, allow_changelog=False):
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
    possible_collections = db.collection_names()
    collections = dict()
    if coll_names is not None:
        if not isinstance(coll_names, list):
            coll_names = [coll_names]
        for database in coll_names:
            if database not in possible_collections:
                client.close()
                exit('Database {} not found!'.format(database))
            if not allow_changelog and database.startswith('__'):
                exit('Queries to database prefixed with __ are VERBOTEN!')
            collections[database] = db[database]
    else:
        collections['repo'] = db['repo']

    return client, db, collections
