#!/usr/bin/env python2
# coding: utf-8
# Distributed under the terms of the MIT License.


""" This script uses qmpy to loop over all entries on a
self-hosted MySQL OQMD database and creates a matador
MongoDB collection containing entries in the matador format.

Instructions for self-hosting the OQMD can be found at:

http://oqmd.org/static/docs/getting_started.html#setting-up-the-database

This script has only been tested with MariaDB 10 (with Centos 7 packages
mariadb-devel, mariadb-server), which are required to install qmpy.

Note:
    As of writing this script, qmpy is Python 2 only, so is not compatible
    with the main matador core.

This script can balloon in memory usage in restart mode. You may wish to
invoke it as `while ./spatulate_oqmd.py --append --no_restart; do :; done;`
to clear the memory.

"""

from random import randint
from os.path import dirname, realpath
import argparse
from traceback import print_exc
import numpy as np
import pymongo as pm
import tqdm
import qmpy


def qmpy_parse_kpoints(string):
    """ Parse VASP kpoint block into grid (ignores offset). """
    return [int(i) for i in string.split('\n')[-2].split(' ')]


def qmpy_entry_to_doc(entry):
    """ Parse a qmpy entry object into a matador document,
    with associated VASP parameters converted into the CASTEP
    keywords.

    Parameters:
        entry (:obj:`qmpy.Entry`): the entry in the OQMD.

    Returns:
        dict: matador document.

    """

    doc = dict()

    if entry is None:
        raise RuntimeError('Entry is None.')
    if entry.duplicate_of is not None and entry.duplicate_of.id != entry.id:
        raise RuntimeError('Structure is a duplicate.')

    doc['task'] = 'geometryoptimization'

    try:
        doc['num_atoms'] = int(entry.natoms)
    except TypeError:
        doc['num_atoms'] = int(entry.structure.natoms)

    doc['elems'] = sorted(list({elem.symbol for elem in entry.elements}))
    doc['stoichiometry'] = sorted([[elem, entry.comp[elem]] for elem in entry.comp], key=lambda x: x[0])
    doc['num_fu'] = doc['num_atoms'] / int(sum(doc['stoichiometry'][i][1] for i in range(len(doc['stoichiometry']))))
    try:
        doc['formation_energy_per_atom'] = entry.energy
        doc['formation_energy'] = entry.energy * doc['num_atoms']
    except TypeError:
        raise RuntimeError('entry.id={} has no energy'.format(entry.id))
    doc['enthalpy_per_atom'] = entry.energy
    doc['enthalpy'] = entry.energy * doc['num_atoms']
    doc['source'] = ['OQMD {}'.format(int(entry.id))]
    doc['_oqmd_entry_id'] = int(entry.id)
    if 'icsd' in entry.keywords:
        doc['icsd'] = int(entry.label.split('-')[-1])
        try:
            doc['reference'] = str(entry.reference.citation)
        except Exception:
            pass

    doc['space_group'] = entry.spacegroup.hm
    doc['lattice_cart'] = entry.structure.cell.tolist()
    doc['positions_abs'] = entry.structure.cartesian_coords.tolist()
    doc['positions_frac'] = entry.structure.coords.tolist()
    doc['atom_types'] = []
    for atom in entry.structure.atoms:
        doc['atom_types'].append(atom.species)
    doc['cell_volume'] = entry.structure.volume
    doc['task'] = 'geometryoptimization'
    doc['species_pot'] = {elem: 'OQMD' for elem in doc['elems']}
    doc['stable'] = entry.stable

    # some entries have multiple E_f, choose the last one,
    # which should match entry.energy: OR SO I THOUGHT
    # instead, we need to loop over the formation energies and find
    # the matching one. This is potentially the lowest one every time.
    efs = entry.formationenergy_set.all()
    for ef in efs:
        if abs(ef.delta_e - entry.energy) < 1e-6:
            eform = ef
            break
    else:
        print('Unable to match formation energy for entry.id={}'.format(entry.id))
        raise RuntimeError
    doc['hull_distance'] = max(0, eform.stability)
    ef_calc = eform.calculation
    oqmd_calc_settings = ef_calc.settings
    doc['cut_off_energy'] = oqmd_calc_settings['encut']
    doc['spin_polarized'] = bool(oqmd_calc_settings['ispin'] - 1)
    doc['kpoints_mp_grid'] = qmpy_parse_kpoints(ef_calc.KPOINTS)
    doc['kpoints_mp_spacing'] = calc_mp_spacing(doc['lattice_cart'], doc['kpoints_mp_grid'])
    doc['xc_functional'] = 'PBE'

    doc['forces'] = ef_calc.output.forces
    doc['pressure'] = 0.0
    doc['stress'] = -0.1 * (ef_calc.output.sxx + ef_calc.output.syy + ef_calc.output.szz)/3.0
    try:
        doc['max_force_on_atom'] = np.max(np.linalg.norm(doc['forces'], axis=-1))
    except Exception:
        print(entry.id, doc['forces'])
    doc['forces'] = doc['forces'].tolist()

    return doc


def calc_mp_spacing(real_lat, mp_grid, prec=3):
    """ Convert real lattice in Cartesian basis and the
    kpoint_mp_grid into a grid spacing. Copied from matador
    utils version.

    Parameters:
        real_lat (list): Cartesian lattice vectors.
        mp_grid (:obj:`list` of :obj:`int`): 3 integers defining the MP grid.

    Keyword arguments:
        prec (int): desired decimal precision of output.

    Returns:
        float: mp_spacing rounded to `prec`.

    """
    recip_lat = real2recip(real_lat)
    recip_len = np.zeros((3))
    recip_len = np.sqrt(np.sum(np.power(recip_lat, 2), axis=1))
    spacing = recip_len / (2*np.pi*np.asarray(mp_grid))
    max_spacing = np.max(spacing)
    exponent = round(np.log10(max_spacing) - prec)
    return round(max_spacing + 0.5*10**exponent, prec)


def real2recip(real_lat):
    """ Convert the real lattice in Cartesian basis to
    the reciprocal space lattice. Copy of matador.utils
    version.

    Parameters:
        real_lat (list): Cartesian lattice vectors.

    Returns:
        list: Cartesian lattice vectors of reciprocal lattice.

    """
    real_lat = np.asarray(real_lat)
    recip_lat = np.zeros((3, 3))
    recip_lat[0] = (2*np.pi)*np.cross(real_lat[1], real_lat[2]) / \
        (np.dot(real_lat[0], np.cross(real_lat[1], real_lat[2])))
    recip_lat[1] = (2*np.pi)*np.cross(real_lat[2], real_lat[0]) / \
        (np.dot(real_lat[1], np.cross(real_lat[2], real_lat[0])))
    recip_lat[2] = (2*np.pi)*np.cross(real_lat[0], real_lat[1]) / \
        (np.dot(real_lat[2], np.cross(real_lat[0], real_lat[1])))
    return recip_lat.tolist()


class DBConverter:
    def __init__(
        self,
        host=None, client=None, dryrun=False, debug=False,
        verbosity=0, db_name='oqmd_1.2', append=False,
        start_id=0, chunk_size=1000, restart=True
    ):
        """ Connect to the relevant databases and
        set off the scraper.
        """
        self.import_count = 0
        self.success_count = 0
        self.limit = 50000
        self.num_scraped = 0
        self.dryrun = dryrun
        self.debug = debug
        self.chunk_size = chunk_size
        self.db_name = db_name
        self.client = client
        self.verbosity = verbosity
        self.append = append
        self.restart = restart
        self.start_id = start_id
        # set up I/O for text_id
        if not self.dryrun:
            print(__file__)
            wordfile = open(dirname(realpath(__file__)) + '/../../matador/scrapers/words', 'r')
            nounfile = open(dirname(realpath(__file__)) + '/../../matador/scrapers/nouns', 'r')
            self.wlines = wordfile.readlines()
            self.num_words = len(self.wlines)
            self.nlines = nounfile.readlines()
            self.num_nouns = len(self.nlines)
            wordfile.close()
            nounfile.close()

            if self.client is None:
                print('connecting to {}'.format(host))
                self.client = pm.MongoClient(host)

            self.db = self.client.crystals
            current_collections = self.db.list_collection_names()
            if self.db_name in current_collections and not self.append:
                raise SystemExit('Desired db_name already exists!')

            self.repo = self.db[self.db_name]
            self.create_indices()
            self.build_mongo()

    def create_indices(self):
        # create unique index for oqmd ID to allow for repeated imports
        self.repo.create_index([('enthalpy_per_atom', pm.ASCENDING)])
        self.repo.create_index([('stoichiometry', pm.ASCENDING)])
        self.repo.create_index([('elems', pm.ASCENDING)])
        self.repo.create_index([('source', pm.ASCENDING)])
        self.create_extra_indices()

    def struct2db(self, struct):
        """ Insert completed Python dictionary into chosen
        database, with generated text_id.
        """
        plain_text_id = [self.wlines[randint(0, self.num_words-1)].strip(),
                         self.nlines[randint(0, self.num_nouns-1)].strip()]
        struct['text_id'] = plain_text_id
        if '_id' in struct:
            raise RuntimeError('{}'.format(struct))
        try:
            _ = self.repo.insert_one(struct)
        except pm.errors.DuplicateKeyError:
            return 0

        return 1


class MPConverter(DBConverter):

    def __init__(self, cursor, *args, **kwargs):
        self.cursor = cursor
        super().__init__(*args, **kwargs)

    def create_extra_indices(self):
        # create unique index for MP ID to allow for repeated imports
        self.repo.create_index('_mp_id', unique=True)

    def build_mongo(self):
        for doc in self.cursor:
            doc = self.doc2entry(doc)
            self.import_count += self.struct2db(doc)

        print(self.import_count, '/', len(self.cursor))

    def doc2entry(self, doc):
        if 'task' not in doc:
            doc['task'] = 'geometryoptimization'
        if 'species_pot' not in doc:
            doc['species_pot'] = {elem: 'OQMD' for elem in doc['elems']}
        for key in doc:
            if isinstance(doc[key], np.ndarray):
                doc[key] = doc[key].tolist()
            elif isinstance(doc[key], set):
                doc[key] = sorted(list(doc[key]))

        return doc


class OQMDConverter(DBConverter):
    """ The OQMDConverter class implements methods to scrape
    the OQMD SQL database for all entries using the qmpy interface.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def create_extra_indices(self):
        # create unique index for oqmd ID to allow for repeated imports
        self.repo.create_index('_oqmd_entry_id', unique=True)

    def build_mongo(self):
        """ Perform QMPY query for all entries, and scrape them into a MongoDB. """
        # start by scraping all converged structures with chosen label
        chunk_iter = 0
        all_structures = qmpy.Entry.objects.all().count()
        print('Expecting {} structures total'.format(all_structures))

        try:
            self.start_id = self.repo.find().sort('_oqmd_entry_id', pm.DESCENDING).limit(1)[0]['_oqmd_entry_id']
            print('Restarting from entry.id={}'.format(self.start_id))
        except IndexError:
            self.start_id = 0

        while self.num_scraped < all_structures and self.import_count < self.limit:
            chunk_min = self.start_id + chunk_iter * self.chunk_size
            chunk_max = chunk_min + self.chunk_size
            chunk_iter += 1

            # check that there are even any structures left to find
            remaining = qmpy.Entry.objects.filter(id__gte=chunk_min).count()
            print('{} structures remaining with ID > {}'.format(remaining, chunk_min))
            if remaining == 0:
                break

            cursor = (qmpy.Entry.objects
                      .filter(id__gte=chunk_min)
                      .filter(id__lt=chunk_max))

            num_structures = cursor.count()
            if num_structures < 1:
                continue

            print('Chunk {} -> {} contains {} entries'.format(chunk_min, chunk_max, num_structures))
            for entry in tqdm.tqdm(cursor, total=num_structures):
                try:
                    doc = {}
                    doc = qmpy_entry_to_doc(entry)
                except RuntimeError:
                    print_exc()
                    continue
                except Exception:
                    print_exc()
                    continue

                self.success_count += 1
                if not self.dryrun:
                    self.import_count += self.struct2db(doc)

            # seems like this keeps memory down in Python 2.7, otherwise it can balloon
            del cursor

            self.num_scraped = self.repo.count_documents({})
            print('Reached {} de-duplicated entiries out of {} including duplicates'.format(self.num_scraped, all_structures))

            if not self.restart:
                break

        if self.dryrun:
            print('Successfully scraped', self.success_count, '/',
                  'structures.')
        if not self.dryrun:
            if self.import_count == 0:
                raise RuntimeError('Nothing imported.')
            print('Successfully imported', self.import_count, '/',
                  'structures.')


if __name__ == '__main__':
    # importer = QMPYConverter()
    parser = argparse.ArgumentParser(
        description='Import OQMD (http://oqmd.org) structures into MongoDB database.',
        epilog='Written by Matthew Evans (2016)')
    parser.add_argument('-d', '--dryrun', action='store_true',
                        help='run the importer without connecting to the database')
    parser.add_argument('-v', '--verbosity', action='count',
                        help='enable verbose output')
    parser.add_argument('--no_restart', action='store_true',
                        help='don\'t restart script if there are missing entries')
    parser.add_argument('--debug', action='store_true',
                        help='enable debug output to print every dict')
    parser.add_argument('--append', action='store_true',
                        help='add to existing collection')
    args = parser.parse_args()
    importer = OQMDConverter(dryrun=args.dryrun,
                             debug=args.debug,
                             db_name='oqmd_1.2',
                             restart=not args.no_restart,
                             append=args.append,
                             verbosity=args.verbosity)
