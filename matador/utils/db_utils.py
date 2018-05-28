# coding: utf-8
# Distributed under the terms of the MIT License.

""" Some simple utilities for making DB connections. """


import pymongo as pm
from matador.config import load_custom_settings


def make_connection_to_collection(coll_names, check_collection=False, allow_changelog=False,
                                  mongo_settings=None, testing=False):
    """ Connect to database of choice.

    Parameters:
        coll_names (str): name of collection.

    Keyword Arguments:
        check_collection (bool): check whether collections exist (forces connection)
        allow_changelog (bool): allow queries to collections with names prefixed by __
        mongo_settings (dict): dict containing mongo and related config
        testing (bool): if testing, then don't ask for user input from stdin

    Returns:
        client (MongoClient): the connection to the database
        db (Database): the database to query
        collections (dict): Collection objects indexed by name

    """

    if mongo_settings is None:
        settings = load_custom_settings()
    else:
        settings = mongo_settings

    print('Trying to connect to {host}:{port}/{db}'.format(**settings['mongo']))

    client = pm.MongoClient(
        host=settings['mongo']['host'],
        port=settings['mongo']['port'],
        connect=False,
        maxIdleTimeMS=600000,  # disconnect after 10 minutes idle
        socketTimeoutMS=10000,  # give up on database after 20 seconds without results
        serverSelectionTimeoutMS=10000,  # give up on server after 2 seconds without results
        connectTimeoutMS=10000)  # give up trying to connect to new database after 2 seconds

    try:
        if settings['mongo']['db'] not in client.database_names():
            if testing:
                response = 'y'
            else:
                response = input('Database {db} does not exist at {host}:{port}/{db}, '
                                 'would you like to create it? (y/n) '
                                 .format(**settings['mongo']))

            if response.lower() != 'y':
                exit('Exiting...')
            else:
                print('Creating database {}'.format(settings['mongo']['db']))
    except pm.errors.ServerSelectionTimeoutError as exc:
        print('{}: {}'.format(type(exc).__name__, exc))
        raise SystemExit('Unable to connect to {host}:{port}/{db}, exiting...'.format(**settings['mongo']))

    db = client[settings['mongo']['db']]
    possible_collections = db.collection_names()
    collections = dict()
    # allow lists of collections for backwards-compat, though normally
    # we only want to connect to one at a time
    if coll_names is not None:
        if not isinstance(coll_names, list):
            coll_names = [coll_names]
        for collection in coll_names:
            if collection not in possible_collections:
                if check_collection:
                    client.close()
                    exit('Collection {} not found!'.format(collection))
                else:
                    print('Creating new collection {}...'.format(collection))
            if not allow_changelog and collection.startswith('__'):
                exit('Queries to collections prefixed with __ are VERBOTEN!')
            collections[collection] = db[collection]
    else:
        default_collection = settings['mongo']['default_collection']
        if default_collection not in possible_collections:
            if check_collection:
                client.close()
                exit('Default collection {} not found!'.format(default_collection))
            else:
                print('Creating new collection {}...'.format(default_collection))
        collections['repo'] = db[default_collection]

    return client, db, collections
