""" This submodule includes functionality to create a
template .matadorrc file.

"""

__all__ = ['quickstart_settings']

import os
from matador.utils.print_utils import print_warning, print_success


def quickstart_settings():
    """ Call this if no user-specified configuration file is found. A series
    of questions will be asked to create a basic matadorrc.

    """
    print('No user-specified input found at either ~/.matadorrc or '
          '~/.config/matadorrc')
    response = input('Would you like to use quickstart to create one now? [y/n] ')
    if not response.lower().startswith('y'):
        return None

    hline = 80*'-'
    print('Okay, let\'s get started!')
    print('Firstly, where would you like me to put your config file: ')
    valid = False
    while not valid:
        response = input('[1] ~/.matadorrc, [2] ~/.config/matadorrc? ')
        if response not in ['1', '2']:
            print('I didn\'t understand, sorry... please enter 1 or 2')
        else:
            valid = True

    if response == '1':
        fname = os.path.expanduser('~/.matadorrc')
    else:
        fname = os.path.expanduser('~/.config/matadorrc')

    flines = []
    print(hline)

    response = input('Great! Would you like to use matador with a MongoDB database? [y/n] ')
    if response.lower().startswith('y'):
        flines.append('mongo:')
        host_response = input('Where is this MongoDB hosted? [default: localhost] ')
        if host_response == '':
            host_response = 'localhost'
        flines.append('  host: {}'.format(host_response))
        print(hline)
        port_response = input('and at what port? [default: 27017] ')
        if port_response == '':
            port_response = 27017
        flines.append('  port: {}'.format(port_response))
        print(hline)
        print('Please set your firewall and MongoDB settings accordingly: '
              'by default MongoDB\nneeds no authentication so be careful '
              'if your database is running on an open port!')
        test_conn_response = input('Would you like me to test the connection to this MongoDB instance now? [y/n] ')
        keep_trying = True
        if test_conn_response.lower() == 'y':
            while keep_trying:
                import pymongo as pm
                cli = pm.MongoClient(host_response,
                                     port=int(port_response),
                                     maxIdleTimeMS=600000,
                                     socketTimeoutMS=3600000,
                                     serverSelectionTimeoutMS=10000,
                                     connectTimeoutMS=10000)
                try:
                    cli.database_names()
                    keep_trying = False
                except pm.errors.ServerSelectionTimeoutError:
                    print_warning('Failed to connect to {}:{}, are you sure it is running?'
                                  .format(host_response, port_response))
                    response = input('Would you like me to try again? [y/n] ')
                    keep_trying = response.lower().startswith('y')
                except Exception:
                    keep_trying = False

        print(hline)
        db_response = input('What is the name of the database (not collection) '
                            'that you want to query/create? [default: crystals] ')
        if db_response == '':
            db_response = 'crystals'
        flines.append('  db: {}'.format(db_response))

        print(hline)
        coll_response = input('What about the default collection name to query? [default: repo] ')
        if coll_response == '':
            coll_response = 'repo'
        flines.append('  default_collection: {}'.format(coll_response))

        print(hline)
        print('Would you like to set a protected folder for this collection?')
        file_path_response = input('Any import to this collection must happen from this path, and the '
                                   '--force flag must be used if so. [default: None] ')
        if file_path_response != '':
            flines.append('  default_file_collection_path: {}'.format(file_path_response))

        print(hline)
        print_success('Mongo section of matadorrc complete! Saving as {} and exiting...'.format(fname))
        print(hline)

    if not flines:
        print('Okay, there is no config file to be made and I will keep using the defaults.')
        return None

    with open(fname, 'a') as f:
        for line in flines:
            f.write(line + '\n')

    return fname
