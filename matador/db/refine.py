# coding: utf-8
# Distributed under the terms of the MIT License.

""" This module contains functionality to update and overwrite database
entries with specific tasks, e.g. symmetry and substructure analysis.
This can be especially useful for "patching" the data in databases
created with older matador versions.

"""


from traceback import print_exc
import os

import pymongo as pm

from matador.utils.print_utils import print_notify, print_warning, print_failure


class Refiner:
    """ Refiner implements methods to alter certain parts of the
    database in place, either in overwrite, set or compare/display mode.
    Current modifiables are space groups, substructures,
    the set of elements, tags and DOIs.

    """
    def __init__(self, cursor, collection=None, task=None, mode='display', **kwargs):
        """ Parses args and initiates modification.

        Parameters:
            cursor (list of dicts): matador cursor to refine.

        Keyword arguments:
            collection (Collection): mongodb collection to query/edit.
            task (str/callable): one of 'sym', 'spg', 'elem_set', 'tag', 'doi' or 'source',
                or a custom function that takes in and returns a cursor and field to modify.
            mode (str): one of 'display', 'overwrite', 'set'.

        """
        possible_tasks = ['sym', 'spg', 'elem_set', 'tag', 'doi', 'source', 'pspot', 'raw']
        possible_modes = ['display', 'overwrite', 'set']

        if mode not in possible_modes:
            print('Mode not understood, defaulting to "display".')
            mode = 'display'
        if collection is None and mode in ['overwrite', 'set']:
            raise SystemExit('Impossible to overwite or set without db collection, exiting...')
        if task is None:
            raise SystemExit('No specified task, exiting...')
        if task not in possible_tasks and not callable(task):
            raise SystemExit('Did not understand task, please choose one of ' + ', '.join(possible_tasks))
        if task == 'tag' and mode == 'set':
            raise SystemExit('Task "tags" and mode "set" will not alter the database, please use mode "overwrite".')

        self.cursor = list(cursor)
        self.diff_cursor = []
        self.collection = collection
        self.mode = mode
        self.changed_count = 0
        self.failed_count = 0
        self.args = kwargs

        if task in ['spg', 'sym']:
            if kwargs.get('symprec'):
                self.symmetry(symprec=kwargs.get('symprec'))
            else:
                self.symmetry()
            self.field = 'space_group'
        elif task == 'elem_set':
            self.elem_set()
            self.field = 'elems'
        elif task == 'tag':
            self.field = 'tags'
            self.tag = self.args.get('new_tag')
            if self.tag is None:
                print_warning('No new tag defined, nothing will be done.')
            else:
                self.add_tag()
        elif task == 'doi':
            self.field = 'doi'
            self.doi = self.args.get('new_doi')
            if self.doi is None:
                print_warning('No new DOI defined, nothing will be done.')
            else:
                self.add_doi()
        elif task == 'source':
            self.field = 'root_source'
            self.add_root_source()
        elif task == 'pspot':
            self.field = 'species_pot'
            self.tidy_pspots()
        elif task == 'raw':
            self.field = '_raw'
            self.add_raw_data()
        elif callable(task):
            print('Using custom task function: {}'.format(task))
            self.diff_cursor, self.field = task(self.cursor)

        self.changed_count = len(self.diff_cursor)
        print(self.changed_count, '/', len(self.cursor), 'to be changed.')
        print(self.failed_count, '/', len(self.cursor), 'failed.')

        if self.mode in ['set', 'overwrite'] and self.changed_count > 0:
            self.update_docs()

    def update_docs(self):
        """ Updates documents in database with correct priority. """
        requests = []
        # if in "set" mode, do not overwrite, just apply
        if self.mode == 'set':
            for _, doc in enumerate(self.diff_cursor):
                requests.append(pm.UpdateOne({'_id': doc['_id'], self.field: {'$exists': False}},
                                             {'$set': {self.field: doc[self.field]}}))
        # else if in overwrite mode, overwrite previous field
        elif self.mode == 'overwrite':
            for _, doc in enumerate(self.diff_cursor):
                requests.append(pm.UpdateOne({'_id': doc['_id']}, {'$set': {self.field: doc[self.field]}}))
        if self.args.get('debug'):
            for request in requests:
                print(request)

        result = self.collection.bulk_write(requests)
        print_notify(str(result.modified_count) + ' docs modified.')

    def symmetry(self, symprec=1e-3):
        """ Compute space group with spglib. """
        from matador.utils.cell_utils import doc2spg
        import spglib as spg
        print('Refining symmetries...')
        if self.mode == 'display':
            print_warning('{}'.format('At symprec: ' + str(symprec)))
            print_warning("{:^36}{:^16}{:^16}".format('text_id', 'new sg', 'old sg'))
        for _, doc in enumerate(self.cursor):
            try:
                spg_cell = doc2spg(doc)
                sg = spg.get_spacegroup(spg_cell, symprec=symprec).split(' ')[0]
                if sg != doc['space_group']:
                    self.changed_count += 1
                    self.diff_cursor.append(doc)
                    if self.mode == 'display':
                        print_notify("{:^36}{:^16}{:^16}"
                                     .format(doc['text_id'][0]+' '+doc['text_id'][1], sg, doc['space_group']))
                    doc['space_group'] = sg
                else:
                    if self.mode == 'display':
                        print("{:^36}{:^16}{:^16}"
                              .format(doc['text_id'][0]+' '+doc['text_id'][1], sg, doc['space_group']))
            except Exception:
                self.failed_count += 1
                if self.args.get('debug'):
                    print_exc()
                    print_failure('Failed for' + ' '.join(doc['text_id']))

    def elem_set(self):
        """ Imbue documents with the set of elements,
        i.e. set(doc['atom_types']), for quicker look-up.

        """
        for doc in self.cursor:
            try:
                doc['elems'] = list(set(doc['atom_types']))
                self.diff_cursor.append(doc)
                self.changed_count += 1
            except Exception as oops:
                if self.args.get('debug'):
                    print(repr(oops))
                self.failed_count += 1

    def add_tag(self):
        """ Add a tag to each document. """
        for _, doc in enumerate(self.cursor):
            try:
                if 'tags' in doc:
                    if doc['tags'] is None:
                        doc['tags'] = list()
                    else:
                        doc['tags'] = list(doc['tags'])
                    doc['tags'].append(self.tag)
                else:
                    doc['tags'] = [self.tag]
                self.diff_cursor.append(doc)
                self.changed_count += 1
            except Exception as error:
                print(repr(error))
                self.failed_count += 1

    def add_doi(self):
        """ Add a doi to each document. """
        if self.doi.count('/') != 1:
            raise SystemExit('Malformed DOI... please use xxxxx/xxxxx format.')
        for _, doc in enumerate(self.cursor):
            try:
                if 'doi' in doc:
                    doc['doi'].append(self.doi)
                else:
                    doc['doi'] = [self.doi]
                self.diff_cursor.append(doc)
                self.changed_count += 1
            except Exception as error:
                print(repr(error))
                self.failed_count += 1

    def add_root_source(self):
        """ Add the "root_source" key to a document in the database,
        i.e. the name of the structure, minus file extension.

        """
        from matador.utils.chem_utils import get_root_source
        for _, doc in enumerate(self.cursor):
            try:
                if 'root_source' in doc:
                    continue
                else:
                    doc['root_source'] = get_root_source(doc['source'])
                    self.diff_cursor.append(doc)
                    self.changed_count += 1
            except Exception as error:
                print(repr(error))
                self.failed_count += 1

    def tidy_pspots(self):
        """ Loop over all documents and make sure they only have pspots
        for the elements that exist in the structure.

        """
        for _, doc in enumerate(self.cursor):
            try:
                marked_for_del = []
                atoms = set(doc['atom_types'])
                for elem in doc['species_pot']:
                    if elem not in atoms:
                        marked_for_del.append(elem)
                for elem in marked_for_del:
                    del doc['species_pot'][elem]
                    self.diff_cursor.append(doc)
                    self.changed_count += 1
            except Exception as error:
                print(repr(error))
                self.failed_count += 1

    def add_raw_data(self):
        """ Loop over all documents in the query and try to open the files
        listed under their `source` fields, storing them under the `_raw` key.

        """

        for _, doc in enumerate(self.cursor):
            try:
                sources = doc['source']
                raw_files = []
                for source in sources:
                    if os.path.isfile(source):
                        with open(source, 'r') as f:
                            raw_files.append((source, f.readlines()))

                doc['_raw'] = raw_files
                self.diff_cursor.append(doc)
                self.changed_count += 1
            except Exception as error:
                print(repr(error))
                self.failed_count += 1
