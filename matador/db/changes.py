# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements an interface for querying the
__changelog_<collection> collections to allow for display and reversion
of particular database changes.

"""


from matador.db import make_connection_to_collection
from matador.utils.print_utils import print_warning, print_notify


class DatabaseChanges:
    """ Class to view and undo particular
    database changesets.

    """
    def __init__(self, collection_name: str, changeset_ind=0, action='view', mongo_settings=None):
        """ Parse arguments and run changes interface.

        Parameters:
            collection_name (str): the base collection name to act upon

        Keyword arguments:
            changset_ind (int): the number of the changset to act upon (1 is oldest)
            action (str): either 'view' or 'undo'

        """
        self.changelog_name = '__changelog_{}'.format(collection_name)
        _, _, self.collections = make_connection_to_collection(self.changelog_name,
                                                               allow_changelog=True,
                                                               mongo_settings=mongo_settings)
        self.repo = [self.collections[key] for key in self.collections][0]
        curs = list(self.repo.find())

        if len(curs) == 0:
            exit('No changesets found for {}'.format(collection_name))

        # if no changeset specified, print summary
        if changeset_ind == 0:
            self.print_change_summary(curs)

        # otherwise, try to act on particular changeset
        elif changeset_ind <= len(curs):
            self.change = curs[changeset_ind-1]
            self.view_changeset(self.change, changeset_ind-1)
            if action == 'undo':
                count = curs[changeset_ind-1]['count']
                print_warning('An attempt will now be made to remove {} structures from {}.'.format(count, collection_name))
                print_notify('Are you sure you want to do that? (y/n)')
                response = input()
                if response.lower() == 'y':
                    print_notify('You don\'t have any doubts at all? (y/n)')
                    reponse = input()
                    if reponse.lower() == 'n':
                        print('You\'re the boss, deleting structures now...')
                    else:
                        exit('As I thought...')
                else:
                    return

                # proceed with deletion
                client, _, collections = make_connection_to_collection(collection_name, allow_changelog=False)
                collection_to_delete_from = [collections[key] for key in collections][0]
                result = collection_to_delete_from.remove({'_id': {'$in': self.change['id_list']}})
                print('Deleted {}/{} successfully.'.format(result['n'], self.change['count']))
                print('Tidying up changelog database...')
                self.repo.remove({'_id': self.change['_id']})
                if self.repo.count() == 0:
                    print('No structures left remaining, deleting database...')
                    collection_to_delete_from.drop()
                    self.repo.drop()
                print('Success!')

    def view_changeset(self, changeset, index):
        """ Prints all details about a particular changeset.

        Parameters:
            changeset (dict): changeset stored in changelog database
            index (int): changeset index

        """
        print('Files added by changeset:')
        for src in changeset['src_list']:
            print('(+) {}'.format(src))
        print('({}) {} {:>7d} additions'.format(index+1, changeset['date'], changeset['count']))

    def print_change_summary(self, curs):
        """ Prints a summary of changes.

        Parameters:
            curs (list): cursor from changelog database

        """
        for index, change in enumerate(curs):
            print('({}) {} {:>7d} additions'.format(index+1, change['date'], change['count']))
