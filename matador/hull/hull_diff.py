# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the diffing of phase diagrams as new structures
are added.

Intended usage:

- Compare current hull with 5 days ago

    matador hulldiff -c KSnP -int --compare 5

- Compare hull a month ago to 2 days ago

    matador hulldiff -c KSnP -int --compare 30 2

"""


from matador.utils.cursor_utils import display_results


class HullDiff:
    """ Takes two QueryConvexHull objects and 'diffs' them, printing a
    +/- line for each new/deleted structure within the specified hull
    cutoff (i.e. compares hull.hull_cursor only).

    """

    def __init__(self, hull_old, hull_new, **args):
        """ Initialise difference from two hulls.

        Parameters:
            hull_old (QueryConvexHull): older hull object.
            hull_new (QueryConvexHull): newer hull object.

        Keyword arguments:
            args: dict, command-line arguments to print results.

        """
        self.id_new = [doc['text_id'] for doc in hull_new.hull_cursor]
        self.id_old = [doc['text_id'] for doc in hull_old.hull_cursor]
        self.additions = []
        self.deletions = []
        self.cursor = []
        # self.cursor = [doc for doc in hull_new.hull_cursor]
        for doc in hull_new.hull_cursor:
            if doc['text_id'] not in self.id_old:
                self.cursor.append(doc)
                self.additions.append(doc['text_id'])
        for doc in hull_old.hull_cursor:
            if doc['text_id'] not in self.id_new:
                self.cursor.append(doc)
                self.deletions.append(doc['text_id'])
        self.cursor = sorted(self.cursor, key=lambda x: (x['enthalpy_per_atom'], x['concentration']))
        self.print_diff(**args)

    def print_diff(self, **args):
        """ Use the display_results function to print the diff. """
        display_results(self.cursor, hull=True, additions=self.additions, deletions=self.deletions)
        if self.additions is not None:
            if args.get('summary'):
                print('{} additions and {} replacements as best at stoichiometry.'
                      .format(len(self.additions), len(self.deletions)))
            else:
                print('{} additions and {} deletions.'.format(len(self.additions), len(self.deletions)))

    def plot_diff(self):
        """ Plot the hull with additions ONLY highlighted in green. """
        raise NotImplementedError


def diff_hulls(client, collections, **args):
    """ Helper function to diff hulls from cmd-line args.

    Parameters:
        client (MongoClient): connection to database
        collections (dict of Collection): dict of collections to query

    Keyword arguments:
        args: dict, command-line arguments to select hulls.

    """
    from matador.query import DBQuery
    from matador.hull import QueryConvexHull
    diff_args = args['compare']
    del args['compare']
    args['no_plot'] = True
    if args.get('verbosity') in [0, None] and not args.get('debug'):
        args['quiet'] = True

    diff_args = sorted(diff_args)
    args['time'] = diff_args[0]
    print('Calculating hull {} days ago...'.format(args['time']))
    query_old = DBQuery(client, collections, **args)
    hull_old = QueryConvexHull(query_old, **args)
    if len(diff_args) == 1 or isinstance(diff_args, str):
        args['time'] = '0'
    else:
        args['time'] = diff_args[1]
    if args['time'] == '0':
        print('... to compare with up-to-date hull.')
    else:
        print('... to compare with hull {} days ago.'.format(args['time']))
    query_new = DBQuery(client, collections, **args)
    hull_new = QueryConvexHull(query_new, **args)
    HullDiff(hull_old, hull_new, **args)
