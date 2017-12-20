# coding: utf-8
""" This file implements the diffing of phase
diagrams as new structures are added.

Intended usage:

- Compare current hull with 5 days ago

    matador hulldiff -c KSnP -int --compare 5

- Compare hull a month ago to 2 days ago

    matador hulldiff -c KSnP -int --compare 30 2

"""
from matador.utils.cursor_utils import display_results


class HullDiff:
    """ Takes two QueryConvexHull objects and 'diffs' them,
    printing a +/- line for each new/deleted structure within
    the specified hull cutoff (i.e. compares hull.hull_cursor only).
    """
    def __init__(self, hull_old, hull_new):
        """ Initialise difference from two hulls. """
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
        self.print_diff()

    def print_diff(self):
        """ Use the display_results function to print the diff. """
        display_results(self.cursor, hull=True, additions=self.additions, deletions=self.deletions)

    def plot_diff(self):
        """ Plot the hull with additions ONLY highlighted in green. """
        raise NotImplementedError
