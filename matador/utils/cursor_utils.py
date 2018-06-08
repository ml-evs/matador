# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule defines some useful generic cursor methods for
displaying, extracting and refining results from a Mongo cursor/list.

"""


from traceback import print_exc
from time import strftime

import numpy as np
from matador.utils.cell_utils import get_spacegroup_spg
from matador.utils.chem_utils import get_formula_from_stoich, get_root_source
from matador import __version__


def display_results(cursor,
                    args=None, argstr=None, additions=None, deletions=None,
                    hull=False, markdown=False, latex=False, use_source=False, colour=True):
    """ Print query results in a cryan-like fashion, optionally
    in a markdown format.

    Parameters:
        cursor (list of dict): list of matador documents

    Keyword arguments:
        args (dict): extra keyword arguments
        argstr (str): string to store matador initialisation command
        additions (list): list of string text_ids to be coloured green with a (+)
        deletions (list): list of string text_ids to be coloured red with a (-)
        hull (bool): whether or not to print hull-style (True) or query-style
        markdown (bool): whether or not to write a markdown file containing results
        latex (bool): whether or not to create a LaTeX table
        use_source (bool): print source instead of text_id
        colour (bool): colour on-hull structures

    Returns:
        str or None: markdown or latex string, if markdown or latex is True, else None.

    """
    if args is None:
        args = dict()

    details = args.get('details')

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
        markdown_string = ('GenDate: ' + strftime("%H:%M %d/%m/%Y") + '\n')
        if argstr is not None:
            markdown_string += ('Command: matador ' + ' '.join(argstr) + '\n')
        markdown_string += ('Version: ' + __version__ + '\n\n')

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
            header_string += "{:^40}".format('ID')
            units_string += "{:^40}".format('')
        else:
            header_string += "{:^28}".format('ID')
            units_string += "{:^28}".format('')
        header_string += "{:^5}".format('!?!')
        units_string += "{:^5}".format('')
    else:
        header_string += "{:^40}".format('Root')
        units_string += "{:^40}".format('')
    header_string += "{:^10}".format('Pressure')
    units_string += "{:^10}".format('(GPa)')
    if args.get('per_atom'):
        header_string += "{:^11}".format('Volume/atom')
    else:
        header_string += "{:^11}".format('Volume/fu')
    units_string += "{:^11}".format('(Ang^3)')
    if hull:
        header_string += "{:^13}".format('Hull dist.')
        units_string += "{:^13}".format('(meV/atom)')
    elif args.get('per_atom'):
        header_string += "{:^18}".format('Enthalpy')
        units_string += "{:^18}".format('(eV/atom)')
    else:
        header_string += "{:^18}".format('Enthalpy')
        units_string += "{:^18}".format('(meV/fu)')
    header_string += "{:^13}".format('Space group')
    header_string += "{:^15}".format('Formula')
    header_string += "{:^8}".format('# fu')
    header_string += "{:^8}".format('Prov.')

    # ensure cursor is sorted by enthalpy
    cursor = sorted(cursor, key=lambda doc: doc['enthalpy_per_atom'], reverse=False)

    if latex:
        latex_sub_style = r'\text'
    else:
        latex_sub_style = ''

    for ind, doc in enumerate(cursor):
        postfix = ''
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

        formula_substring = get_formula_from_stoich(doc['stoichiometry'],
                                                    tex=latex,
                                                    latex_sub_style=latex_sub_style)

        if 'encapsulated' in doc:
            formula_substring += '+CNT'
        if last_formula != formula_substring:
            gs_enthalpy = 0.0
        formula_string.append(formula_substring)
        if not markdown:
            if hull and np.abs(doc.get('hull_distance')) <= 0.0 + 1e-12:
                if use_source:
                    src = get_root_source(doc['source'])
                    if colour:
                        struct_string.append("\033[92m\033[1m* {:<38}".format(src))
                        postfix = '\033[0m'
                    else:
                        struct_string.append("* {:<38}".format(src))
                else:
                    if colour:
                        struct_string.append(
                            "\033[92m\033[1m* {:^26}".format(doc['text_id'][0]+' '+doc['text_id'][1]))
                        postfix = '\033[0m'
                    else:
                        struct_string.append("* {:^26}".format(doc['text_id'][0] + ' ' + doc['text_id'][1]))
            else:
                if use_source:
                    src = get_root_source(doc['source'])
                    struct_string.append("  {:<38}".format(src))
                else:
                    struct_string.append("  {:^26}".format(doc['text_id'][0] + ' ' + doc['text_id'][1]))
            if additions is not None and doc['text_id'] in additions:
                struct_string[-1] = '\033[92m\033[1m' + ' + ' + struct_string[-1]
                postfix = '\033[0m'
            elif deletions is not None and doc['text_id'] in deletions:
                struct_string[-1] = '\033[91m\033[1m' + ' - ' + struct_string[-1]
                postfix = '\033[0m'
            try:
                if doc.get('prototype'):
                    struct_string[-1] += "{:^5}".format('*p*')
                elif doc['quality'] == 0:
                    struct_string[-1] += "{:^5}".format('!!!')
                else:
                    struct_string[-1] += "{:^5}".format((5 - doc['quality']) * '?')
            except KeyError:
                struct_string[-1] += "{:5}".format(' ')
        else:
            struct_string.append("{:40}".format(get_root_source(doc['source'])))

        if 'pressure' in doc and doc['pressure'] != 'xxx':
            struct_string[-1] += "{: >9.2f}".format(doc['pressure'])
        else:
            struct_string[-1] += "{:^10}".format('xxx')
        if args.get('per_atom') and 'cell_volume' in doc and 'num_atoms' in doc:
            struct_string[-1] += "{:>11.1f}".format(doc['cell_volume'] / doc['num_atoms'])
        elif 'cell_volume' in doc and 'num_fu' in doc:
            struct_string[-1] += "{:>11.1f}".format(doc['cell_volume'] / doc['num_fu'])
        else:
            struct_string[-1] += "{:^11}".format('xxx')
        try:
            if hull:
                struct_string[-1] += "{:>13.1f}".format(
                    0 if doc.get('hull_distance') <= 1e-12 else 1000 * doc.get('hull_distance')
                )
            elif args.get('per_atom'):
                struct_string[-1] += "{:>18.5f}".format(doc['enthalpy_per_atom'] - gs_enthalpy)
            else:
                struct_string[-1] += "{:>18.5f}".format(doc['enthalpy'] / doc['num_fu'] - gs_enthalpy)
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

        struct_string[-1] += postfix

        if latex:
            latex_struct_string.append("{:^30} {:^10} & ".format(formula_substring, '$\\star$'
                                                                 if doc['hull_distance'] == 0 else ''))
            latex_struct_string[-1] += ("{:^20.0f} & ".format(doc.get('hull_distance') * 1000)
                                        if doc.get('hull_distance') > 0 else '{:^20} &'.format('-'))
            latex_struct_string[-1] += ("{:^20.0f} & ".format(doc['gravimetric_capacity'])
                                        if doc.get('hull_distance') == 0 else '{:^20} &'.format('-'))
            latex_struct_string[-1] += "{:^20} & ".format(get_spacegroup_spg(doc))
            prov = get_guess_doc_provenance(doc['source'], doc.get('icsd'))
            if doc.get('icsd'):
                prov += ' {}'.format(doc['icsd'])
            latex_struct_string[-1] += "{:^30} & ".format(prov)
            latex_struct_string[-1] += "{:^30} \\\\".format('')

        if last_formula != formula_substring:
            if args.get('per_atom'):
                gs_enthalpy = doc['enthalpy'] / doc['num_atoms']
            else:
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
                source_string[-1] += ' ' + file.split('structure_repository')[-1]
                if num != len(doc['source']) - 1:
                    source_string[-1] += '\n'

    if not markdown or latex:
        print(len(header_string) * '─')
        print(header_string)
        print(units_string)
        print(len(header_string) * '─')

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
            elif latex:
                latex_string += latex_struct_string[ind] + '\n'
            else:
                print(substring)
                if details:
                    print(detail_string[ind])
                    print(detail_substring[ind])
                if args.get('source'):
                    print(source_string[ind])
                if details or args.get('source'):
                    print(len(header_string) * '─')
    if markdown:
        return markdown_string
    elif latex:
        latex_string += '\\end{tabular}'
        return latex_string

    return None


def set_cursor_from_array(cursor, array, key):
    """ Updates the key-value pair for documents in
    internal cursor from a numpy array.
    """
    assert len(array) == len(cursor) or len(array) - 1 == len(cursor)
    for ind, _ in enumerate(cursor):
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
    except KeyError:
        print_exc()
    array = np.asarray(array)
    if len(array) != len(cursor):
        raise RuntimeError('Some keys were missing.')
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
