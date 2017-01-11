# coding: utf-8
""" This file defines some useful generic cursor methods. """
import numpy as np
from traceback import print_exc
from time import strftime
from pkg_resources import require
try:
    __version__ = require("matador")[0].version
except:
    from version import __version__
__version__ = __version__.strip()


def display_results(cursor, args=None, hull=False, markdown=False):
    """ Print query results in a cryan-like fashion, optionally
    in a markdown format.
    """
    details = args.get('details')
    if args is None:
        args = dict()
    struct_string = []
    detail_string = []
    detail_substring = []
    source_string = []
    formula_string = []
    last_formula = ''
    header_string = ''

    if markdown:
        markdown_string = ('Query generated with matador ' + __version__ +
                           ' at ' + strftime("%I:%M%S") + ' on ' + strftime("%d/%m/%Y") + '\n\n')

    if not markdown:
        header_string += "{:^24}".format('ID')
        header_string += "{:^5}".format('!?!')
    else:
        header_string += "{:^30}".format('Root')
    header_string += "{:^12}".format('Pressure')
    if args.get('per_atom'):
        header_string += "{:^12}".format('Volume/atom')
    else:
        header_string += "{:^12}".format('Volume/fu')
    if hull:
        header_string += "{:^18}".format('Hull dist./atom')
    elif args.get('per_atom'):
        header_string += "{:^18}".format('Enthalpy/atom')
    else:
        header_string += "{:^18}".format('Enthalpy/fu')
    header_string += "{:^12}".format('Space group')
    header_string += "{:^10}".format('Formula')
    header_string += "{:^8}".format('# fu')
    header_string += "{:^8}".format('Prov.')

    for ind, doc in enumerate(cursor):
        formula_substring = ''
        if 'phase' in doc:
            if 'alpha' in doc['phase']:
                formula_substring += 'α-'
            elif 'beta' in doc['phase']:
                formula_substring += 'β-'
            elif 'gamma' in doc['phase']:
                formula_substring += 'γ-'
            elif 'theta' in doc['phase']:
                formula_substring += 'θ-'
        atom_per_fu = 0
        for item in doc['stoichiometry']:
            for item_ind, subitem in enumerate(item):
                if item_ind == 0:
                    formula_substring += str(subitem)
                if item_ind == 1:
                    if subitem != 1:
                        formula_substring += str(subitem)
                    atom_per_fu += subitem
        if 'encapsulated' in doc:
            formula_substring += '+CNT'
        if last_formula != formula_substring:
            gs_enthalpy = 0.0
        formula_string.append(formula_substring)
        if not markdown:
            if hull and np.abs(doc.get('hull_distance')) <= 0.0 + 1e-12:
                struct_string.append(
                    '* ' + "{:^22}".format(doc['text_id'][0]+' '+doc['text_id'][1]))
            else:
                struct_string.append(
                    "{:^24}".format(doc['text_id'][0]+' '+doc['text_id'][1]))
            try:
                if doc['quality'] == 0:
                    struct_string[-1] += "{:^5}".format('!!!')
                else:
                    struct_string[-1] += "{:^5}".format((5-doc['quality'])*'?')
            except:
                struct_string[-1] += "{:5}".format(' ')
        else:
            struct_string.append("{:30}".format(next(source.split('/')[-1].split('.')[0] for source in doc['source']
                                                if (source.endswith('.res') or source.endswith('.castep') or
                                                    source.endswith('.history') or source.endswith('.history.gz')))))
        try:
            struct_string[-1] += "{:^12.3f}".format(doc['pressure'])
        except:
            struct_string[-1] += "{:^12}".format('xxx')
        try:
            if args.get('per_atom'):
                struct_string[-1] += "{:^12.3f}".format(doc['cell_volume']/doc['num_atoms'])
            else:
                struct_string[-1] += "{:^12.3f}".format(doc['cell_volume']/doc['num_fu'])
        except:
            struct_string[-1] += "{:^12}".format('xxx')
        try:
            if hull:
                struct_string[-1] += "{:^18.5f}".format(doc.get('hull_distance'))
                # struct_string[-1] += "{:^18.5f}".format(doc.get('formation_enthalpy_per_atom'))
            elif args.get('per_atom'):
                struct_string[-1] += "{:^18.5f}".format(doc['enthalpy_per_atom'])
            else:
                struct_string[-1] += "{:^18.5f}".format(doc['enthalpy']/doc['num_fu'] -
                                                        gs_enthalpy)
        except:
            struct_string[-1] += "{:^18}".format('xxx')
        try:
            struct_string[-1] += "{:^12}".format(doc['space_group'])
        except:
            struct_string[-1] += "{:^12}".format('xxx')
        struct_string[-1] += "{:^10}".format(formula_substring)
        try:
            struct_string[-1] += "{:^8}".format(doc['num_fu'])
        except:
            struct_string[-1] += "{:^8}".format('xxx')
        try:
            prov = get_guess_doc_provenance(doc['source'], doc.get('icsd'))
            struct_string[-1] += "{:^8}".format(prov)
        except:
            struct_string[-1] += "{:^8}".format('xxx')

        if last_formula != formula_substring:
            gs_enthalpy = doc['enthalpy'] / doc['num_fu']
        last_formula = formula_substring

        if details:
            detail_string.append(11 * ' ' + u"├╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ ")
            if args.get('source'):
                detail_substring.append(11 * ' ' + u"├╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ ")
            else:
                detail_substring.append(11 * ' ' + u"└╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ ")
            if 'spin_polarized' in doc:
                if doc['spin_polarized']:
                    detail_string[-1] += 'S-'
            if 'sedc_scheme' in doc:
                detail_string[-1] += doc['sedc_scheme'].upper()+'+'
            if 'xc_functional' in doc:
                detail_string[-1] += doc['xc_functional']
            else:
                detail_string[-1] += 'xc-functional unknown!'
            if 'cut_off_energy' in doc:
                detail_string[-1] += ', ' + "{:4.2f}".format(doc['cut_off_energy']) + ' eV'
            else:
                detail_string[-1] += 'cutoff unknown'
            if 'external_pressure' in doc:
                detail_string[-1] += (', ' +
                                      "{:4.2f}".format(doc['external_pressure'][0][0]) +
                                      ' GPa')
            if 'kpoints_mp_spacing' in doc:
                detail_string[-1] += ', ~' + str(doc['kpoints_mp_spacing']) + ' 1/A'
            if 'species_pot' in doc:
                try:
                    for species in doc['species_pot']:
                        detail_substring[-1] += doc['species_pot'][species] + ', '
                except:
                    pass
            if 'icsd' in doc:
                detail_substring[-1] += 'ICSD-CollCode' + doc['icsd'] + ', '
            if 'tags' in doc:
                try:
                    for tag in doc['tags']:
                        detail_substring[-1] += tag + ', '
                except:
                    pass
            if 'user' in doc:
                detail_substring[-1] += doc['user']
            if 'encapsulated' in doc:
                try:
                    detail_string[-1] += (', (n,m)=(' + str(doc['cnt_chiral'][0]) +
                                          ',' + str(doc['cnt_chiral'][1]) + ')')
                    detail_string[-1] += ', r=' + "{:4.2f}".format(doc['cnt_radius']) + ' A'
                    detail_string[-1] += ', z=' + "{:4.2f}".format(doc['cnt_length']) + ' A'
                except:
                    pass
            detail_string[-1] += ' ' + (len(header_string)-len(detail_string[-1])-1)*u"╌"
            detail_substring[-1] += ' ' + (len(header_string)-len(detail_substring[-1])-1)*u"╌"

        if args.get('source'):
            if len(doc['source']) == 1:
                source_string.append(11*' ' + u"└──────────────────")
            else:
                source_string.append(11*' ' + u"└───────────────┬──")
            for num, file in enumerate(doc['source']):
                if len(doc['source']) == 1:
                    source_string[-1] += ''
                elif num == len(doc['source'])-1:
                    source_string[-1] += (len(u"└────────────── ")+11)*' ' + u'└──'
                elif num != 0:
                    source_string[-1] += (len(u"└────────────── ")+11)*' ' + u'├──'
                source_string[-1] += ' ' + file.split('structure_repository')[-1]
                if num != len(doc['source'])-1:
                    source_string[-1] += '\n'

    if not markdown:
        print(len(header_string)*'─')
        print(header_string)
        print(len(header_string)*'─')
    else:
        markdown_string += header_string + '\n'
        markdown_string += len(header_string)*'-' + '\n'
    if args.get('summary'):
        current_formula = ''
        formula_list = []
        count = 0
        for ind, substring in enumerate(formula_string):
            if substring != current_formula and substring not in formula_list:
                count += 1
                if markdown:
                    markdown_string += struct_string[ind] + '\n'
                else:
                    print(struct_string[ind])
                if details and not markdown:
                    print(detail_string[ind])
                    print(detail_substring[ind])
                if args.get('source') and not markdown:
                    print(source_string[ind])
                current_formula = substring
                formula_list.append(substring)
    else:
        for ind, substring in enumerate(struct_string):
            if markdown:
                markdown_string += struct_string[ind] + '\n'
            else:
                print(substring)
            if details and not markdown:
                print(detail_string[ind])
                print(detail_substring[ind])
            if args.get('source') and not markdown:
                print(source_string[ind])
            if details or args.get('source') and not markdown:
                print(len(header_string) * '─')
    if markdown:
        return markdown_string


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


def get_guess_doc_provenance(sources, icsd=None):
    """ Returns a guess at the provenance of a structure
    from its source list.

    Return possiblities are 'ICSD', 'SWAP', 'OQMD' or
    'AIRSS'.
    """
    prov = 'AIRSS'
    if sources is str:
        sources = [sources]
    for fname in sources:
        if (fname.endswith('.castep') or fname.endswith('.res') or
                fname.endswith('.history') or 'OQMD' in fname):
            if 'swap' in fname.lower():
                prov = 'SWAPS'
            elif icsd is not None:
                prov = 'ICSD'
            elif 'oqmd' in fname.lower():
                prov = 'OQMD'
            elif 'collcode' in fname.lower():
                if fname.split('/')[-1].count('-') == 2:
                    prov = 'SWAPS'
                else:
                    prov = 'ICSD'
            elif '-icsd' in fname.lower():
                prov = 'ICSD'
    return prov


def get_spg_uniq(cursor, symprec=1e-2, latvecprec=1e-3, posprec=1e-3):
    """ Use spglib to find duplicate structures in a cursor.
    Returns uniq_list and same_list.

    * cursor     : list of matador structure docs.
    * symprec    : spglib symmetry precision for cell standardisation.
    * latvecprec : tolerance on lattice vectors.
    * posprec    : tolerance on fractional atomic positions.
    * uniq_list  : list of indices of the unique structures in cursor.
    * same_list  : list of pairs indices of duplicate structures.
    """

    from utils.cell_utils import doc2spg
    import spglib as spg

    spg_cursor = list()
    for doc in cursor:
        spg_cursor.append(doc2spg(doc))

    refined_list = []
    for crystal in spg_cursor:
        refined_list.append(spg.standardize_cell(crystal, to_primitive=False, no_idealize=False, symprec=symprec))
    for i in range(len(refined_list)):
        for j in range(len(refined_list[i][1])):
            for k in range(len(refined_list[i][1][j])):
                if refined_list[i][1][j][k] > 1-1e-10:
                    refined_list[i][1][j][k] = 0.0
    for i in range(len(refined_list)):
        refined_list[i] = (refined_list[i][0],
                           refined_list[i][1][np.argsort(refined_list[i][1][:, 0])],
                           refined_list[i][2][np.argsort(refined_list[i][1][:, 0])])
    uniq_list = np.arange(0, len(spg_cursor))
    same_list = []
    shift_list = []
    for i in range(len(spg_cursor)):
        for j in range(i+1, len(spg_cursor)):
            if cursor[i]['stoichiometry'] == cursor[j]['stoichiometry']:
                if np.allclose(refined_list[i][0], refined_list[j][0], atol=latvecprec, rtol=0):
                    if np.allclose(refined_list[i][1], refined_list[j][1], atol=posprec, rtol=0):
                        same_list.append((i, j))
                    else:
                        for dim in range(3):
                            if not rigid_shift(refined_list[i], refined_list[j], dim, posprec):
                                    break
                            elif dim == 3:
                                same_list.append((i, j))
                                shift_list.append((i, j, dim))
                            break
    dupes = list(set([pair[1] for pair in same_list]))
    uniq_list = np.delete(uniq_list, dupes)
    print(len(dupes), 'duplicates found and removed.')
    print(len(shift_list), 'of which were shifted cells.')
    return uniq_list


def rigid_shift(structA, structB, dim, posprec):
    assert(len(structA[2]) == len(structB[2]))
    shift_array = structA[1][:, dim] - structB[1][:, dim]
    # if trivial zero shift, return True
    if np.all((np.abs(shift_array)) < 1e-4):
        return True
    shift_array[np.where(shift_array < 0)] += 1
    # if nontrivial shift, return True
    return np.all((np.abs(np.diff(shift_array)) < 1e-4))


def filter_cursor(cursor, key, min, max):
    """ Returns a cursor obeying the filter on the given key. """
    filtered_cursor = list()
    print('Filtering', key, min, max)
    for doc in cursor:
        try:
            if doc[key] < max and doc[key] >= min:
                filtered_cursor.append(doc)
        except:
            pass
    return filtered_cursor
