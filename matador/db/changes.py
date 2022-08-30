# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements an interface for querying the
__changelog_<collection> collections to allow for display and reversion
of particular database changes.

"""


from matador.db import make_connection_to_collection
from matador.utils.print_utils import print_warning, print_notify


class DatabaseChanges:
    """Class to view and undo particular
    database changesets.

    """

    def __init__(
        self,
        collection_name: str,
        changeset_ind=0,
        action="view",
        override=False,
        mongo_settings=None,
    ):
        """Parse arguments and run changes interface.

        Parameters:
            collection_name (str): the base collection name to act upon

        Keyword arguments:
            changset_ind (int): the number of the changset to act upon (1 is oldest)
            action (str): either 'view' or 'undo'
            override (bool): override all options to positive answers for testing
            mongo_settings (dict): dictionary of already-sources mongo settings

        """
        self.changelog_name = "__changelog_{}".format(collection_name)
        _, _, self.collections = make_connection_to_collection(
            self.changelog_name,
            allow_changelog=True,
            override=override,
            mongo_settings=mongo_settings,
        )
        self.repo = [self.collections[key] for key in self.collections][0]
        curs = list(self.repo.find())

        if not curs:
            exit("No changesets found for {}".format(collection_name))

        # if no changeset specified, print summary
        if changeset_ind == 0:
            self.print_change_summary(curs)

        elif changeset_ind > len(curs):
            exit(
                'No changeset {} found for collection called "{}".'.format(
                    changeset_ind, collection_name
                )
            )

        # otherwise, try to act on particular changeset
        elif changeset_ind <= len(curs):
            self.change = curs[changeset_ind - 1]
            self.view_changeset(self.change, changeset_ind - 1)
            if action == "undo":
                count = curs[changeset_ind - 1]["count"]
                print_warning(
                    "An attempt will now be made to remove {} structures from {}.".format(
                        count, collection_name
                    )
                )
                print_notify("Are you sure you want to do that? (y/n)")
                if override:
                    response = "y"
                else:
                    response = input()
                if response.lower() == "y":
                    print_notify("You don't have any doubts at all? (y/n)")
                    if override:
                        next_response = "n"
                    else:
                        next_response = input()
                    if next_response.lower() == "n":
                        print("You're the boss, deleting structures now...")
                    else:
                        exit("As I thought...")
                else:
                    exit()

                # proceed with deletion
                _, _, collections = make_connection_to_collection(
                    collection_name, allow_changelog=False, override=override
                )
                collection_to_delete_from = [collections[key] for key in collections][0]
                result = collection_to_delete_from.delete_many(
                    {"_id": {"$in": self.change["id_list"]}}
                )
                print(
                    "Deleted {}/{} successfully.".format(
                        result.deleted_count, self.change["count"]
                    )
                )
                print("Tidying up changelog database...")
                self.repo.delete_one({"_id": self.change["_id"]})
                if not self.repo.find_one():
                    print("No structures left remaining, deleting database...")
                    collection_to_delete_from.drop()
                    self.repo.drop()
                print("Success!")

    @staticmethod
    def view_changeset(changeset, index):
        """Prints all details about a particular changeset.

        Parameters:
            changeset (dict): changeset stored in changelog database
            index (int): changeset index

        """
        print("Files added by changeset:")
        for src in changeset["src_list"]:
            print("(+) {}".format(src))
        print(
            "({}) {} {:>7d} additions".format(
                index + 1, changeset["date"], changeset["count"]
            )
        )

    @staticmethod
    def print_change_summary(curs):
        """Prints a summary of changes.

        Parameters:
            curs (list): cursor from changelog database

        """
        for index, change in enumerate(curs):
            print(
                "({}) {} {:>7d} additions".format(
                    index + 1, change["date"], change["count"]
                )
            )
