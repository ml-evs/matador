#!/usr/bin/env python
# coding: utf-8
""" This file implements the Crystal class, a wrapper
to the raw dictionary stored in MongoDB that allows for validation,
manipulation and analysis of the lattice.
"""
from matador.similarity.pdf_similarity import PDF
from matador.utils.cell_utils import cart2abc, cart2volume


class Crystal(object):
    """ Class that wraps the MongoDB document, providing useful
    interfaces for cell manipulation and validation.
    """
    def __init__(self, doc):
        self.doc = doc
        self.__dict__.update(self.doc)

    def __getitem__(self, key):
        if key not in self.doc:
            return self.__getattribute__(key)
        else:
            return self.doc[key]

    def __setitem__(self, key, item):
        self.doc[key] = item

    @property
    def lattice_cart(self):
        return self['lattice_cart']

    @property
    def lattice_abc(self):
        if 'lattice_abc' not in self.__dict__:
            return cart2abc(self['lattice_cart'])
        else:
            return self['lattice_abc']

    @property
    def cell_volume(self):
        if 'cell_volume' not in self.__dict__:
            print('Recomputing cell volume')
            return cart2volume(self.lattice_cart)
        else:
            return self['cell_volume']

    @property
    def pdf(self, **kwargs):
        self.doc['pdf'] = PDF(self.doc, dr=0.01, gaussian_width=0.01, **kwargs)
        return self.doc['pdf']


if __name__ == '__main__':

    from matador.query import DBQuery
    doc = DBQuery(composition=['LiP'], top=1).cursor[0]
    crystal = Crystal(doc)
    print(crystal.cell_volume)
    del crystal.__dict__['cell_volume']
    print(crystal.cell_volume)
    print(crystal['cell_volume'])
