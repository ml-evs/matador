# coding: utf-8
""" This file defines some useful generic cursor methods. """
import numpy as np

def set_cursor_from_array(cursor, array, key):
    """ Updates the key-value pair for documents in
    internal cursor from a numpy array.
    """
    assert(len(array) == len(cursor) or len(array) - 1 == len(cursor))
    for ind, doc in enumerate(cursor):
        cursor[ind][key] = array[ind]
    return

def get_array_from_cursor(cursor, key):
    """ Returns a numpy array of the values of a key
    in a cursor.
    """
    array = []
    try:
        for doc in cursor:
            array.append(doc[key])
    except:
        print_exc()
    array = np.asarray(array)
    assert(len(array) == len(cursor))
    return array

def filter_cursor(cursor, key, min, max):
    """ Returns a cursor obeying the filter on the given key. """
    filtered_cursor = list()
    print('Filtering', key, min, max)
    for doc in cursor:
        try:
            if doc[key] < max and doc[key] > min:
                filtered_cursor.append(doc)
        except:
            pass
    return filtered_cursor
