# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule defines some useful generic cursor methods for
displaying, extracting and refining results from a Mongo cursor/list.

"""


from time import strftime

import numpy as np
from matador.utils.cell_utils import get_spacegroup_spg
from matador.utils.chem_utils import get_formula_from_stoich, get_root_source, get_stoich_from_formula
from matador import __version__

EPS = 1e-12


def recursive_get(data, keys, top=True):
    """ Recursively slice a nested dictionary by a
    list of keys.

    Parameters:
        data (dict): nested dictionary to get from.
        keys (list): list of keys/indices to delve into.

    Raises:
        KeyError: if any in chain keys are missing,
        IndexError: if any element of a sublist is
            missing.

    """
    if not isinstance(keys, (list, tuple)):
        return data[keys]

    if isinstance(keys, (list, tuple)) and len(keys) == 1:
        return data[keys[0]]
    try:
        return recursive_get(data[keys[0]], keys[1:], top=False)
    except (KeyError, IndexError) as exc:
        if top:
            raise type(exc)('Recursive keys {} missing'.format(keys))
        else:
            raise exc


def recursive_set(data, keys, value):
    """ Recursively slice a nested dictionary by a
    list of keys and set the value.

    Parameters:
        data (dict): nested dictionary to get from.
        keys (list): list of keys/indices to delve into.
        value: value to store under key.

    Raises:
        KeyError: if any intermediate keys are missing.

    """
    if isinstance(keys, (list, tuple)):
        if len(keys) == 1:
            data[keys[0]] = value
        else:
            return recursive_set(data[keys[0]], keys[1:], value)
    else:
        data[keys] = value


def display_results(cursor,
                    energy_key='enthalpy_per_atom', args=None, argstr=None, additions=None, deletions=None, no_sort=False,
                    hull=False, markdown=False, latex=False, colour=True, return_str=False):
    """ Print query results in a cryan-like fashion, optionally
    in a markdown format.

    Parameters:
        cursor (list of dict): list of matador documents

    Keyword arguments:
        argstr (str): string to store matador initialisation command
        additions (list): list of string text_ids to be coloured green with a (+)
            or, list of indices referring to those structures in the cursor.
        deletions (list): list of string text_ids to be coloured red with a (-)
            or, list of indices referring to those structures in the cursor.
        hull (bool): whether or not to print hull-style (True) or query-style
        markdown (bool): whether or not to write a markdown file containing results
        latex (bool): whether or not to create a LaTeX table
        energy_key (str or list): key (or recursive key) to print as energy (per atom)
        colour (bool): colour on-hull structures
        return_str (bool): return string instead of printing.

    Returns:
        str or None: markdown or latex string, if markdown or latex is True, else None.

    """
    if args is None:
        args = dict()

    details = args.get('details')
    use_source = args.get('use_source')

    add_index_mode = False
    del_index_mode = False
    if additions:
        if isinstance(additions[0], int):
            add_index_mode = True
    if deletions:
        if isinstance(deletions[0], int):
            del_index_mode = True

    if add_index_mode:
        assert max(additions) <= len(cursor) and min(additions) >= 0
    if del_index_mode:
        assert max(deletions) <= len(cursor) and min(deletions) >= 0

    if markdown and latex:
        raise RuntimeError('Cannot specify both latex and markdown output at once.')

    struct_string = []
    detail_string = []
    detail_substring = []
    source_string = []
    formula_string = []
    last_formula = ''
    units_string = ''

    if not cursor:
        raise RuntimeError('No structures found in cursor.')

    if markdown:
        markdown_string = 'Date: {}  \n'.format(strftime("%H:%M %d/%m/%Y"))
        if argstr is not None:
            markdown_string += 'Command: matador {}  \n'.format(' '.join(argstr))
        markdown_string += 'Version: {}  \n\n'.format(__version__)

    if latex:
        latex_string = (
            "\\begin{tabular}{l r r c l l}\n"
            "\\rowcolor{gray!20}\n"
            "formula & \\thead{$\\Delta E$\\\\(meV/atom)} & \\thead{$Q_m$\\\\(mAh/g)} & "
            "spacegroup & provenance & description \\\\ \n\n"
        )
        latex_struct_string = []

    header_string = ''
    if additions is not None or deletions is not None:
        header_string += '   '
    if not markdown:
        if use_source:
            header_string += "{:^40}".format('Source')
            units_string += "{:^40}".format('')
        else:
            header_string += "{:^28}".format('ID')
            units_string += "{:^28}".format('')

        header_string += "{:^5}".format('!?!')
        units_string += "{:^5}".format('')
    else:
        header_string += '```\n'
        header_string += "{:^40}".format('Root')
        units_string += "{:^40}".format('')
    header_string += "{:^10}".format('Pressure')
    units_string += "{:^10}".format('(GPa)')

    if args.get('per_atom'):
        header_string += "{:^11}".format('Volume/atom')
    else:
        header_string += "{:^11}".format('Volume/fu')
    units_string += "{:^11}".format('(Ang^3)')

    if hull and args.get('eform'):
        header_string += "{:^18}".format('Formation energy')
        units_string += "{:^18}".format('(eV/atom)')
    elif hull:
        header_string += "{:^18}".format('Hull dist.')
        units_string += "{:^18}".format('(meV/atom)')
    elif args.get('per_atom'):
        header_string += "{:^18}".format(' '.join(energy_key.replace('_per_atom', '').split('_')).title())
        units_string += "{:^18}".format('(eV/atom)')
    else:
        header_string += "{:^18}".format(' '.join(energy_key.replace('_per_atom', '').split('_')).title())
        units_string += "{:^18}".format('(meV/fu)')

    header_string += "{:^13}".format('Space group')
    header_string += "{:^15}".format('Formula')
    header_string += "{:^8}".format('# fu')
    header_string += "{:^8}".format('Prov.')

    # ensure cursor is sorted by enthalpy
    if not no_sort:
        cursor = sorted(cursor, key=lambda doc: recursive_get(doc, energy_key), reverse=False)

    if latex:
        latex_sub_style = r'\text'
    else:
        latex_sub_style = ''

    for ind, doc in enumerate(cursor):
        postfix = ''
        prefix = ''
        formula_substring = ''

        formula_substring = get_formula_from_stoich(doc['stoichiometry'],
                                                    tex=latex,
                                                    latex_sub_style=latex_sub_style)

        if 'encapsulated' in doc:
            formula_substring += '+CNT'
        if last_formula != formula_substring:
            gs_enthalpy = 0.0
        formula_string.append(formula_substring)
        if not markdown:
            if use_source:
                src = get_root_source(doc['source'])
                max_len = 34
                struct_string.append("  {:<38.{max_len}}".format(src if len(src) < max_len else src[:max_len-4]+'[..]', max_len=max_len))
            else:
                struct_string.append("  {:^26.22}".format(' '.join(doc.get('text_id', ['xxx', 'yyy']))))
            if hull and np.abs(doc.get('hull_distance')) <= 0.0 + 1e-12:
                if colour:
                    prefix = '\033[92m\033[1m'
                    postfix = '\033[0m'
                struct_string[-1] = '*' + struct_string[-1][1:]

            if additions is not None:
                if (add_index_mode and ind in additions) or doc.get('text_id', '_') in additions:
                    struct_string[-1] = '+' + struct_string[-1][1:]
                    if colour:
                        prefix = '\033[92m\033[1m'
                        postfix = '\033[0m'
            if deletions is not None:
                if (del_index_mode and ind in deletions) or doc.get('text_id', '_') in deletions:
                    struct_string[-1] = '-' + struct_string[-1][1:]
                    if colour:
                        prefix = '\033[91m\033[1m'
                        postfix = '\033[0m'
            try:
                if doc.get('prototype'):
                    struct_string[-1] += "{:^5}".format('*p*')
                elif doc['quality'] == 0:
                    struct_string[-1] += "{:^5}".format('!!!')
                else:
                    struct_string[-1] += "{:^5}".format((5 - doc['quality']) * '?')
            except KeyError:
                struct_string[-1] += "{:^5}".format(' ')
        else:
            struct_string.append("{:40}".format(get_root_source(doc['source'])))

        if 'pressure' in doc and doc['pressure'] != 'xxx':
            struct_string[-1] += "{: >9.2f}".format(doc['pressure'])
        else:
            struct_string[-1] += "{:^9}".format('xxx')
        try:
            if args.get('per_atom') and 'cell_volume' in doc and 'num_atoms' in doc:
                struct_string[-1] += "{:>11.1f}".format(doc['cell_volume'] / doc['num_atoms'])
            elif 'cell_volume' in doc and 'num_fu' in doc:
                struct_string[-1] += "{:>11.1f}".format(doc['cell_volume'] / doc['num_fu'])
            else:
                struct_string[-1] += "{:^11}".format('xxx')
        except Exception:
            struct_string[-1] += "{:^11}".format('xxx')
        try:
            if hull and args.get('eform'):
                struct_string[-1] += "{:>13.3f}".format(
                    doc['formation_' + energy_key]
                )
            elif hull:
                struct_string[-1] += "{:>13.1f}".format(
                    1000 * doc['hull_distance']
                )
            elif args.get('per_atom'):
                struct_string[-1] += "{:>18.5f}".format(recursive_get(doc, energy_key) - gs_enthalpy)
            else:
                struct_string[-1] += "{:>18.5f}".format(recursive_get(doc, energy_key) * doc['num_atoms'] / doc['num_fu'] - gs_enthalpy)
        except KeyError:
            struct_string[-1] += "{:^18}".format('xxx')

        if 'space_group' in doc:
            struct_string[-1] += "{:^13}".format(doc['space_group'])
        else:
            struct_string[-1] += "{:^13}".format('xxx')

        struct_string[-1] += "{:^15}".format(formula_substring)

        if 'num_fu' in doc:
            struct_string[-1] += "{:^8}".format(int(doc['num_fu']))
        else:
            struct_string[-1] += "{:^8}".format('xxx')

        if 'source' in doc:
            prov = get_guess_doc_provenance(doc['source'], doc.get('icsd'))
            struct_string[-1] += "{:^8}".format(prov)
        else:
            struct_string[-1] += "{:^8}".format('xxx')

        struct_string[-1] = prefix + struct_string[-1] + postfix

        if latex:
            latex_struct_string.append("{:^30} {:^10} & ".format(formula_substring, '$\\star$'
                                                                 if doc.get('hull_distance', 0.1) == 0 else ''))
            latex_struct_string[-1] += ("{:^20.0f} & ".format(doc.get('hull_distance') * 1000)
                                        if doc.get('hull_distance', 0) > 0 else '{:^20} &'.format('-'))
            latex_struct_string[-1] += ("{:^20.0f} & ".format(doc.get('gravimetric_capacity', '-'))
                                        if doc.get('hull_distance', 0.1) == 0 else '{:^20} &'.format('-'))
            latex_struct_string[-1] += "{:^20} & ".format(get_spacegroup_spg(doc))
            prov = get_guess_doc_provenance(doc['source'], doc.get('icsd'))
            if doc.get('icsd'):
                prov += ' {}'.format(doc['icsd'])
            latex_struct_string[-1] += "{:^30} & ".format(prov)
            latex_struct_string[-1] += "{:^30} \\\\".format('')

        if last_formula != formula_substring:
            if args.get('per_atom'):
                gs_enthalpy = recursive_get(doc, energy_key)
            else:
                gs_enthalpy = recursive_get(doc, energy_key) * doc['num_atoms'] / doc['num_fu']

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
                detail_string[-1] += doc['sedc_scheme'].upper() + '+'
            if 'xc_functional' in doc:
                detail_string[-1] += doc['xc_functional']
            else:
                detail_string[-1] += 'xc-functional unknown!'
            if 'cut_off_energy' in doc:
                detail_string[-1] += ', ' + "{:4.2f}".format(doc['cut_off_energy']) + ' eV'
            else:
                detail_string[-1] += 'cutoff unknown'
            if 'external_pressure' in doc:
                detail_string[-1] += (', ' + "{:4.2f}".format(doc['external_pressure'][0][0]) + ' GPa')
            if 'kpoints_mp_spacing' in doc:
                detail_string[-1] += ', ~' + str(doc['kpoints_mp_spacing']) + ' 1/A'
            if 'geom_force_tol' in doc:
                detail_string[-1] += ', {:.2f} eV/A, '.format(doc['geom_force_tol'])
            if 'species_pot' in doc:
                try:
                    for species in doc['species_pot']:
                        detail_substring[-1] += doc['species_pot'][species] + ', '
                except KeyError:
                    pass
            if 'icsd' in doc:
                detail_substring[-1] += 'ICSD-CollCode {}, '.format(doc['icsd'])
            if 'tags' in doc:
                try:
                    if isinstance(doc['tags'], list):
                        for tag in doc['tags']:
                            detail_substring[-1] += tag + ', '
                except KeyError:
                    pass
            if 'user' in doc:
                detail_substring[-1] += doc['user']
            if 'encapsulated' in doc:
                try:
                    detail_string[-1] += (
                        ', (n,m)=(' + str(doc['cnt_chiral'][0]) + ',' + str(doc['cnt_chiral'][1]) + ')'
                    )
                    detail_string[-1] += ', r=' + "{:4.2f}".format(doc['cnt_radius']) + ' A'
                    detail_string[-1] += ', z=' + "{:4.2f}".format(doc['cnt_length']) + ' A'
                except KeyError:
                    pass
            detail_string[-1] += ' ' + (len(header_string) - len(detail_string[-1]) - 1) * u"╌"
            detail_substring[-1] += ' ' + (len(header_string) - len(detail_substring[-1]) - 1) * u"╌"

        if args.get('source'):
            if len(doc['source']) == 1:
                source_string.append(11 * ' ' + u"└──────────────────")
            else:
                source_string.append(11 * ' ' + u"└───────────────┬──")
            for num, file in enumerate(doc['source']):
                if len(doc['source']) == 1:
                    source_string[-1] += ''
                elif num == len(doc['source']) - 1:
                    source_string[-1] += (len(u"└────────────── ") + 11) * ' ' + u'└──'
                elif num != 0:
                    source_string[-1] += (len(u"└────────────── ") + 11) * ' ' + u'├──'
                # source_string[-1] += ' ' + file.split('structure_repository')[-1]
                source_string[-1] += ' ' + file
                if num != len(doc['source']) - 1:
                    source_string[-1] += '\n'

    total_string = ''

    if not markdown and not latex:
        total_string += len(header_string) * '─' + '\n'
        total_string += header_string + '\n'
        total_string += units_string + '\n'
        total_string += len(header_string) * '─' + '\n'

    if markdown:
        markdown_string += len(header_string) * '-' + '\n'
        markdown_string += header_string + '\n'
        markdown_string += units_string + '\n'
        markdown_string += len(header_string) * '-' + '\n'

    if args.get('summary'):
        current_formula = ''
        formula_list = []
        count = 0
        for ind, substring in enumerate(formula_string):
            if substring != current_formula and substring not in formula_list:
                count += 1
                if markdown:
                    markdown_string += struct_string[ind] + '\n'
                elif latex:
                    latex_string += latex_struct_string[ind]
                else:
                    total_string += ('{}\n'.format(struct_string[ind]))
                if details and not markdown:
                    total_string += detail_string[ind] + '\n'
                    total_string += detail_substring[ind] + '\n'
                if args.get('source') and not markdown:
                    total_string += source_string[ind] + '\n'
                current_formula = substring
                formula_list.append(substring)
    else:
        for ind, substring in enumerate(struct_string):
            if markdown:
                markdown_string += struct_string[ind] + '\n'
            elif latex:
                latex_string += latex_struct_string[ind] + '\n'
            else:
                total_string += substring + '\n'
                if details:
                    total_string += detail_string[ind] + '\n'
                    total_string += detail_substring[ind] + '\n'
                if args.get('source'):
                    total_string += source_string[ind] + '\n'
                if details or args.get('source'):
                    total_string += len(header_string) * '─' + '\n'
    if markdown:
        markdown_string += '```'
        return markdown_string
    if latex:
        latex_string += '\\end{tabular}'
        return latex_string
    if return_str:
        return total_string

    print(total_string)


def loading_bar(iterable, width=80):
    """ Checks if tqdm exists and makes a loading bar, otherwise
    just returns initial iterable.

    Parameters:
        iterable (iterable): the thing to be iterated over.

    Keyword arguments:
        width (int): maximum number of columns to use on screen.

    Returns:
        iterable: the decorated iterator.

    """
    try:
        import tqdm
        return tqdm.tqdm(iterable, ncols=width)
    except ImportError:
        return iterable


def set_cursor_from_array(cursor, array, key):
    """ Updates the key-value pair for documents in
    internal cursor from a numpy array.
    """
    if len(array) != len(cursor):
        raise RuntimeError('Trying to fit array of shape {} into cursor of length {}'
                           .format(np.shape(array), len(cursor)))
    for ind, _ in enumerate(cursor):
        recursive_set(cursor[ind], key, array[ind])


def get_array_from_cursor(cursor, key, pad_missing=False):
    """ Returns a numpy array of the values of a key
    in a cursor, where the key can be defined as list
    of keys to use with `recursive_get`.

    Parameters:
        cursor (list): list of matador dictionaries.
        key (str or list): the key to extract, or list
            of keys/subkeys/indices to extract with
            recursive_get.

    Keyword arguments:
        pad_missing (bool): whether to fill array with NaN's
            where data is missing.
    Raises:
        KeyError: if any document is missing that key,
            unless pad_missing is True.

    Returns:
        np.ndarray: numpy array containing results, padded
            with np.nan if key is missing and pad_missing is True.
    """
    array = []
    for ind, doc in enumerate(cursor):
        try:
            if isinstance(key, (tuple, list)):
                array.append(recursive_get(doc, key))
            else:
                array.append(doc[key])
        except KeyError as exc:
            print('{} missing  in entry {}, with source {}'.format(key, ind, doc.get('source')))
            if pad_missing:
                array.append(np.NaN)
            else:
                raise exc
    array = np.asarray(array)
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
        fname_with_folder = fname
        fname = fname.split('/')[-1].lower()
        if (fname.endswith('.castep') or fname.endswith('.res') or fname.endswith('.history') or
                ('oqmd' in fname and fname.count('.') == 0)):
            if any(substr in fname for substr in ['collcode', 'colcode', 'collo']):
                if fname.count('-') == 2 + fname.count('oqmd') or 'swap' in fname:
                    prov = 'SWAPS'
                else:
                    prov = 'ICSD'
            elif 'swap' in fname_with_folder:
                prov = 'SWAPS'
            elif '-ga-' in fname:
                prov = 'GA'
            elif icsd is not None:
                prov = 'ICSD'
            elif 'oqmd' in fname:
                prov = 'OQMD'
            elif '-icsd' in fname:
                prov = 'ICSD'
            elif '-mp-' in fname:
                prov = 'MP'
            elif '-doi-' in fname:
                prov = 'DOI'

    return prov


def filter_cursor(cursor, key, vals, verbosity=0):
    """ Returns a cursor obeying the filter on the given key. Any
    documents that are missing the key will not be returned. Any
    documents with values that cannot be compared to floats will also
    not be returned.

    Parameters:
        cursor (list): list of dictionaries to filter.
        key (str): key to filter.
        vals (list): either 1 value to 2 values to use as a range.
            The values are interpreted as floats for comparison.

    Returns:
        list: list of dictionaries that pass the filter.

    """
    filtered_cursor = list()
    orig_cursor_len = len(cursor)
    if not isinstance(vals, list):
        vals = [vals]
    if len(vals) == 2:
        min_val = float(vals[0])
        max_val = float(vals[1])
        if verbosity > 0:
            print('Filtering {} <= {} < {}'.format(min_val, key, max_val))
        for doc in cursor:
            try:
                if doc[key] < max_val and doc[key] >= min_val:
                    filtered_cursor.append(doc)
            except (TypeError, ValueError, KeyError):
                pass
    else:
        min_val = float(vals[0])
        if verbosity > 0:
            print('Filtering {} >= {}'.format(key, min_val))
        for doc in cursor:
            try:
                if doc[key] >= min_val:
                    filtered_cursor.append(doc)
            except (TypeError, ValueError, KeyError):
                pass

    if verbosity > 0:
        print(orig_cursor_len, 'filtered to', len(filtered_cursor), 'documents.')
    return filtered_cursor


def filter_cursor_by_chempots(species, cursor):
    """ For the desired chemical potentials, remove any incompatible structures
    from cursor.

    Parameters:
        species (list): list of chemical potential formulae.
        cursor (list): list of matador documents to filter.

    Returns:
        list: the filtered cursor.

    """
    # filter out structures with any elements with missing chem pots
    chempot_stoichiometries = []
    for label in species:
        chempot_stoichiometries.append(get_stoich_from_formula(label))

    inds_to_remove = set()
    for ind, doc in enumerate(cursor):
        from matador.utils.chem_utils import get_number_of_chempots
        try:
            cursor[ind]['num_chempots'] = get_number_of_chempots(doc, chempot_stoichiometries)
            cursor[ind]['concentration'] = (cursor[ind]['num_chempots'][:-1] /
                                            np.sum(cursor[ind]['num_chempots'])).tolist()
            for idx, conc in enumerate(cursor[ind]['concentration']):
                if conc < 0 + EPS:
                    cursor[ind]['concentration'][idx] = 0.0
                elif conc > 1 - EPS:
                    cursor[ind]['concentration'][idx] = 1.0

        except RuntimeError:
            inds_to_remove.add(ind)

    return [doc for ind, doc in enumerate(cursor) if ind not in inds_to_remove]
