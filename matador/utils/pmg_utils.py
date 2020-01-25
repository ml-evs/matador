# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements some light wrappers to
the pymatgen, via ASE.

"""

import pickle
import os
import getpass

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen import MPRester
from matador.utils.ase_utils import ase2dict, doc2ase
from matador.utils.cell_utils import calc_mp_spacing

REQ_PROPERTIES = [
    'formula',
    'input',
    'doi',
    'e_above_hull',
    'snl',
    'material_id',
    'icsd_ids',
    'pf_ids',
    'formation_energy_per_atom',
    'structure'
]


def get_chemsys(elements, dumpfile=None):
    """ Scrape the Materials Project for the chemical system specified
    by elements, e.g. for elements `['A', 'B', 'C']` the query performed
    is `chemsys='A-B-C' & nelements=3`. Requires interactive user input
    of their MP API key (unless set by environment variable
    $PMG_MAPI_KEY).


    Parameters:
        elements (list): list of chemical symbols.

    Keyword arguments:
        dumpfile (str): optional filename to dump the pickled response
            after conversion to matador documents.

    """
    cursor = []
    count = {}
    api_key = os.environ.get('PMG_MAPI_KEY', None)
    if api_key is None:
        api_key = getpass.getpass('MP API key:').strip()

    with MPRester(api_key=api_key) as m:
        data = {
            'criteria': {
                'chemsys': '-'.join(elements),
                'nelements': len(elements),
            },
            'properties': REQ_PROPERTIES
        }
        response = m.query(**data)
        for result in response:
            cursor.append(mp2dict(result))

    if dumpfile:
        with open(dumpfile, 'wb') as f:
            pickle.dump(cursor, f)

    print('ICSD structures: {}'.format(sum(1 for doc in cursor if 'icsd' in doc)))
    print('PF structures: {}'.format(sum(1 for doc in cursor if '_mp_pf_ids' in doc)))

    return cursor, count


def doc2pmg(doc):
    """ Converts matador document to a pymatgen structure,
    via ASE.

    """
    ase_atoms = doc2ase(doc)
    return AseAtomsAdaptor.get_structure(ase_atoms)


def mp2dict(response):
    """ Convert a response from pymatgen.MPRester into a matador document,
    via an ASE atoms object. Expects certain properties to be requested
    in order to construct a full matador document, e.g. `structure` & `input`.

    Parameters:
        response (dict): containing one item of the MPRester response.

    """

    if 'structure' not in response:
        raise RuntimeError('`structure` key missing, nothing to scrape...')

    ase_struct = AseAtomsAdaptor.get_atoms(response['structure'])
    doc = ase2dict(ase_struct)

    doc['source'] = []

    if 'material_id' in response:
        doc['_mp_id'] = response['material_id']
        doc['source'].append('{}'.format(doc['_mp_id']))

    if 'pf_ids' in response:
        if response['pf_ids']:
            doc['_mp_pf_ids'] = response['pf_ids']
            for _id in doc['_mp_pf_ids']:
                doc['source'].append('pf-{}'.format(_id))

    if 'icsd_ids' in response:
        if response['icsd_ids']:
            doc['icsd'] = response['icsd_ids'][0]
        doc['_mp_icsd_ids'] = response['icsd_ids']

    if 'formation_energy_per_atom' in response:
        doc['formation_energy_per_atom'] = response['formation_energy_per_atom']
        doc['enthalpy_per_atom'] = response['formation_energy_per_atom']
        doc['enthalpy'] = doc['enthalpy_per_atom'] * doc['num_atoms']
        doc['formation_energy'] = doc['formation_energy_per_atom'] * doc['num_atoms']

    if 'e_above_hull' in response:
        doc['hull_distance'] = response['e_above_hull']

    if 'doi' in response:
        doc['_mp_doi'] = response['doi']

    if 'snl' in response:
        doc['_mp_references'] = response['snl'].references

    if 'input' in response:
        incar = response['input']['incar']
        if not isinstance(incar, dict):
            incar = incar.as_dict()
        doc['cut_off_energy'] = incar['ENCUT']
        doc['xc_functional'] = 'PBE'
        kpts = response['input']['kpoints']
        doc['kpoints_mp_grid'] = kpts.kpts[0]
        doc['kpoints_mp_offset'] = kpts.kpts_shift
        doc['pressure'] = 0.0
        doc['stress'] = 0.0
        doc['kpoints_mp_spacing'] = calc_mp_spacing(doc['lattice_cart'], doc['kpoints_mp_grid'])
        doc['spin_polarized'] = bool(incar['ISPIN'] - 1)

    return doc
