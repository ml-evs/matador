# coding: utf-8
""" This file implements the Crystal class, a wrapper
to the raw dictionary stored in MongoDB that allows for validation,
manipulation and analysis of the lattice.
"""
from matador.similarity.pdf_similarity import PDF


class Crystal(object):
    """ Class that wraps the MongoDB document, providing useful
    interfaces for cell manipulation and validation.
    """
    def __init__(self, doc, a, b, d):
        self.doc = doc
        self.__dict__.update(self.doc)
        self.a = a
        self.b = b
        self.d = d
        # self.__c = None

    def __getitem__(self, key):
        if key not in self.doc:
            return self.__getattribute__(key)
        else:
            return self.doc[key]

    def __setitem__(self, key, item):
        self.doc[key] = item
    @property
    def c(self):
        print('Recomputing c..')
        return self.a + self.b
    @property
    def pdf(self, **kwargs):
        print('Recomputing pdf..')
        self.doc['pdf'] = PDF(self.doc, **kwargs)
        return self.doc['pdf']
    # @c.setter
    # def c(self, **kwargs):
        # print('Recomputing c..')
        # self.__c = self.a + self.b

if __name__ == '__main__':
    from matador.query import DBQuery
    doc = DBQuery(formula='Li3P', top=1).cursor[0]
    cryst = Crystal(doc, 5, 3, 10)
    print(cryst.c)
    print(cryst['pdf'])
    # cryst.pdf.plot_pdf()
    print(cryst['lattice_cart'])
    print(cryst.lattice_cart)

    cryst.lattice_cart[0][0] = cryst.lattice_cart[0][0] * 0.5
    print(cryst['lattice_cart'])
    cryst['lattice_cart'][0][0] = cryst.lattice_cart[0][0] * 0.5
    print(cryst.lattice_cart)
    print(cryst['pdf'])
    print(cryst.geom_method)
    cryst.geom_method = 'lbfgs'
    print(cryst['geom_method'])
    # cryst.d += 1
    # print(cryst.c)
    # print(cryst.pdf)
    # assert cryst['geom_method'] == 'lbfgs'
    # print(cryst.pdf)
    # print(cryst['lattice_cart'])
    # # cryst['pdf'].plot_pdf()
    # print(cryst.pdf)
    # cryst['lattice_cart'][0][0] = cryst['lattice_cart'][0][0] + 10
    # print(cryst['lattice_cart'])
    # print(cryst['pdf'])
    # cryst['pdf'].plot_pdf()
