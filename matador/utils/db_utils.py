# encoding: utf-8
""" Some simple utilities for DB connections. """
import os
import pymongo as pm

MONGO_DEFAULTS = {'mongo': {'db': 'crystals',
                            'host': 'node1',
                            'default_collection': 'repo'}}


def load_custom_settings(config_fname=None):
    """ Load mongodb settings dict from file given by fname, or from defaults.

    Keyword Arguments:

        fname (str): filename of custom settings file.

    """
    import json
    from traceback import print_exc
    if config_fname is None:
        config_fname = '/'.join(__file__.split('/')[:-1]) + '/../../config/matador_conf.json'
        print('Loading settings from {}'.format(config_fname))

    custom_settings = {}
    if os.path.isfile(config_fname):
        try:
            with open(config_fname, 'r') as fp:
                custom_settings = json.load(fp)
        except:
            print_exc()
            exit('Failed to read custom settings file {}'.format(config_fname))

    settings = {}
    settings.update(MONGO_DEFAULTS)
    settings.update(custom_settings)
    return settings


def make_connection_to_collection(coll_names, check_collection=False, allow_changelog=False, mongo_settings=None):
    """ Connect to database of choice.

    Parameters:

        coll_names (str): name of collection.

    Keyword Arguments:

        check_collection (bool): check whether collections exist (forces connection)
        allow_changelog (bool): allow queries to collections with names prefixed by __
        mongo_settings (dict): dict containing mongo and related config

    Returns:

        client (MongoClient): the connection to the database
        db (Database): the database to query
        collections (dict): Collection objects indexed by name

    """

    if mongo_settings is None:
        settings = load_custom_settings()
    else:
        settings = mongo_settings

    client = pm.MongoClient(settings['mongo']['host'])
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
