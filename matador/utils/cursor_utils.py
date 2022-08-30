# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule defines some useful generic cursor methods for
displaying, extracting and refining results from a Mongo cursor/list.

"""

from time import strftime

import numpy as np
import pymongo as pm
from matador.utils.chem_utils import (
    get_formula_from_stoich,
    get_root_source,
    get_stoich_from_formula,
    get_subscripted_formula,
)
from matador import __version__

EPS = 1e-12


def recursive_get(data, keys, _top=True):
    """Recursively slice a nested dictionary by a
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
        return recursive_get(data[keys[0]], keys[1:], _top=False)
    except (KeyError, IndexError) as exc:
        if _top:
            raise type(exc)("Recursive keys {} missing".format(keys))
        raise exc


def recursive_set(data, keys, value):
    """Recursively slice a nested dictionary by a
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


def display_results(
    cursor,
    energy_key="enthalpy_per_atom",
    summary=False,
    args=None,
    argstr=None,
    additions=None,
    deletions=None,
    sort=True,
    hull=False,
    markdown=False,
    latex=False,
    colour=True,
    return_str=False,
    use_source=True,
    details=False,
    per_atom=False,
    eform=False,
    source=False,
    **kwargs,
):
    """Print query results in a table, with many options for customisability.

    TODO: this function has gotten out of control and should be rewritten.

    Parameters:
        cursor (list of dict or pm.cursor.Cursor): list of matador documents

    Keyword arguments:
        summary (bool): print a summary per stoichiometry, that uses the lowest
            energy phase (requires `sort=True`).
        argstr (str): string to store matador initialisation command
        eform (bool): prepend energy key with "formation_".
        sort (bool): sort input cursor by the value of energy key.
        return_str (bool): return string instead of printing.
        details (bool): print extra details as an extra line per structure.
        per_atom (bool): print quantities per atom, rather than per fu.
        source (bool): print all source files associated with the structure.
        use_source (bool): use the source instead of the text id when displaying a structure.
        hull (bool): whether or not to print hull-style (True) or query-style
        energy_key (str or list): key (or recursive key) to print as energy (per atom)
        markdown (bool): whether or not to write a markdown file containing results
        latex (bool): whether or not to create a LaTeX table
        colour (bool): colour on-hull structures
        additions (list): list of string text_ids to be coloured green with a (+)
            or, list of indices referring to those structures in the cursor.
        deletions (list): list of string text_ids to be coloured red with a (-)
            or, list of indices referring to those structures in the cursor.
        kwargs (dict): any extra args are ignored.

    Returns:
        str or None: markdown or latex string, if markdown or latex is True, else None.

    """

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
        latex = False

    # lists in which to accumulate the table
    struct_strings = []
    detail_strings = []
    detail_substrings = []
    source_strings = []
    formulae = []

    # tracking the last formula
    last_formula = ""

    if not cursor:
        raise RuntimeError("No structures found in cursor.")

    if markdown:
        markdown_string = "Date: {}  \n".format(strftime("%H:%M %d/%m/%Y"))
        if argstr is not None:
            markdown_string += "Command: matador {}  \n".format(" ".join(argstr))
        markdown_string += "Version: {}  \n\n".format(__version__)

    if latex:
        latex_string = (
            "\\begin{tabular}{l r r c l l}\n"
            "\\rowcolor{gray!20}\n"
            "formula & "
            "\\thead{$\\Delta E$\\\\(meV/atom)} & "
            "spacegroup & "
            "provenance & "
            "description \\\\ \n\n"
        )
        latex_struct_strings = []
        latex_sub_style = r"\mathrm"
    else:
        latex_sub_style = ""

    header_string, units_string = _construct_header_string(
        markdown, use_source, per_atom, eform, hull, summary, energy_key
    )

    if summary and isinstance(cursor, pm.cursor.Cursor):
        raise RuntimeError("Unable to provide summary when displaying cursor object.")

    # ensure cursor is sorted by enthalpy
    if sort and isinstance(cursor, pm.cursor.Cursor):
        print("Unable to check sorting of cursor, assuming it is already sorted.")
    elif sort:
        sorted_inds = sorted(
            enumerate(cursor), key=lambda element: recursive_get(element[1], energy_key)
        )
        cursor = [ind[1] for ind in sorted_inds]
        sorted_inds = [ind[0] for ind in sorted_inds]
        if additions is not None and add_index_mode:
            additions = [sorted_inds.index(ind) for ind in additions]
        if deletions is not None and del_index_mode:
            deletions = [sorted_inds.index(ind) for ind in deletions]

    # loop over structures and create pretty output
    for ind, doc in enumerate(cursor):

        # use the formula to see if we need to update gs_enthalpy for this formula
        formula_substring = get_formula_from_stoich(
            doc["stoichiometry"], tex=latex, latex_sub_style=latex_sub_style
        )

        if not latex:
            formula_substring = get_subscripted_formula(formula_substring)

        if "encapsulated" in doc:
            formula_substring += "+CNT"
        if last_formula != formula_substring:
            gs_enthalpy = 0.0
        formulae.append(formula_substring)

        struct_strings.append(
            _construct_structure_string(
                doc,
                ind,
                formula_substring,
                gs_enthalpy,
                use_source,
                colour,
                hull,
                additions,
                deletions,
                add_index_mode,
                del_index_mode,
                energy_key,
                per_atom,
                eform,
                markdown,
                latex,
            )
        )

        if latex:
            latex_struct_strings.append(
                "{:^30} {:^10} & ".format(
                    formula_substring,
                    "$\\star$" if doc.get("hull_distance") == 0 else "",
                )
            )
            latex_struct_strings[-1] += (
                "{:^20.0f} & ".format(doc.get("hull_distance") * 1000)
                if doc.get("hull_distance", 0) > 0
                else "{:^20} &".format("-")
            )
            latex_struct_strings[-1] += "{:^20} & ".format(
                doc.get("space_group", "xxx")
            )
            prov = get_guess_doc_provenance(doc["source"], doc.get("icsd"))
            if doc.get("icsd"):
                prov += " {}".format(doc["icsd"])
            latex_struct_strings[-1] += "{:^30} & ".format(prov)
            latex_struct_strings[-1] += "{:^30} \\\\".format("")

        if last_formula != formula_substring:
            if per_atom:
                gs_enthalpy = recursive_get(doc, energy_key)
            else:
                gs_enthalpy = (
                    recursive_get(doc, energy_key) * doc["num_atoms"] / doc["num_fu"]
                )

        last_formula = formula_substring

        if details:
            detail_string, detail_substring = _construct_detail_strings(
                doc, padding_length=len(header_string), source=source
            )
            detail_strings.append(detail_string)
            detail_substrings.append(detail_substring)

        if source:
            source_strings.append(_construct_source_string(doc["source"]))

    total_string = ""
    total_string += len(header_string) * "─" + "\n"
    total_string += header_string + "\n"
    total_string += units_string + "\n"
    total_string += len(header_string) * "─" + "\n"

    if markdown:
        markdown_string += len(header_string) * "-" + "\n"
        markdown_string += header_string + "\n"
        markdown_string += units_string + "\n"
        markdown_string += len(header_string) * "-" + "\n"

    summary_inds = []
    # filter for lowest energy phase per stoichiometry
    if summary:
        current_formula = ""
        formula_list = {}
        for ind, substring in enumerate(formulae):
            if substring != current_formula and substring not in formula_list:
                current_formula = substring
                formula_list[substring] = 0
                summary_inds.append(ind)
            formula_list[substring] += 1
    else:
        summary_inds = range(len(struct_strings))

    # construct final string containing table
    if markdown:
        markdown_string += "\n".join(struct_strings[ind] for ind in summary_inds)
    elif latex:
        latex_string += "\n".join(latex_struct_strings[ind] for ind in summary_inds)
    else:
        for ind in summary_inds:
            total_string += struct_strings[ind] + "\n"
            if details:
                total_string += detail_strings[ind] + "\n"
                total_string += detail_substrings[ind] + "\n"
            if source:
                total_string += source_strings[ind] + "\n"
            if details or source:
                total_string += len(header_string) * "─" + "\n"

    if markdown:
        markdown_string += "```"
        return markdown_string

    if latex:
        latex_string += "\\end{tabular}"
        return latex_string

    if return_str:
        return total_string

    print(total_string)


def loading_bar(iterable, width=80, verbosity=0):
    """Checks if tqdm exists and makes a loading bar, otherwise
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

        if verbosity < 1:
            raise RuntimeError
        return tqdm.tqdm(iterable, ncols=width)
    except (ImportError, RuntimeError):
        return iterable


def set_cursor_from_array(cursor, array, key):
    """Updates the key-value pair for documents in
    internal cursor from a numpy array.
    """
    if len(array) != len(cursor):
        raise RuntimeError(
            "Trying to fit array of shape {} into cursor of length {}".format(
                np.shape(array), len(cursor)
            )
        )
    for ind, _ in enumerate(cursor):
        recursive_set(cursor[ind], key, array[ind])


def get_array_from_cursor(cursor, key, pad_missing=False):
    """Returns a numpy array of the values of a key
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
            print(
                "{} missing  in entry {}, with source {}".format(
                    key, ind, doc.get("source")
                )
            )
            if pad_missing:
                array.append(np.NaN)
            else:
                raise exc
    array = np.asarray(array)
    return array


def get_guess_doc_provenance(sources, icsd=None):
    """Returns a guess at the provenance of a structure
    from its source list.

    Return possiblities are 'ICSD', 'SWAP', 'OQMD' or
    'AIRSS', 'MP' or 'PF'.
    """
    prov = "AIRSS"
    if isinstance(sources, dict):
        sources = sources["source"]
    elif isinstance(sources, str):
        sources = [sources]
    for fname in sources:
        fname_with_folder = fname
        fname = fname.split("/")[-1].lower()
        if (
            fname.endswith(".castep")
            or fname.endswith(".res")
            or fname.endswith(".history")
            or fname.endswith(".phonon")
            or fname.count(".") == 0
        ):
            if any(substr in fname for substr in ["collcode", "colcode", "collo"]):
                if fname.count("-") == 2 + fname.count("oqmd") or "swap" in fname:
                    prov = "SWAPS"
                else:
                    prov = "ICSD"
            elif "swap" in fname_with_folder:
                prov = "SWAPS"
            elif "-ga-" in fname:
                prov = "GA"
            elif icsd is not None:
                prov = "ICSD"
            elif "oqmd" in fname:
                prov = "OQMD"
            elif "-icsd" in fname:
                prov = "ICSD"
            elif "pf-" in fname and prov is None:
                prov = "PF"
            elif any(s in fname for s in ["mp-", "mp_"]) and prov != "PF":
                prov = "MP"
            elif "-sm-" in fname:
                prov = "SM"
            elif "-doi-" in fname:
                prov = "DOI"
            elif "-config_enum" in fname:
                prov = "ENUM"

    return prov


def filter_unique_structures(cursor, quiet=False, **kwargs):
    """Wrapper for `matador.fingerprints.similarity.get_uniq_cursor` that
    displays the results and returns the filtered cursor.

    """
    from matador.fingerprints.similarity import get_uniq_cursor

    uniq_inds, dupe_dict, _, _ = get_uniq_cursor(cursor, **kwargs)
    filtered_cursor = [cursor[ind] for ind in uniq_inds]

    if not quiet:
        display_cursor = []
        additions = []
        deletions = []
        for key in dupe_dict:
            additions.append(len(display_cursor))
            display_cursor.append(cursor[key])
            if dupe_dict[key]:
                for _, jnd in enumerate(dupe_dict[key]):
                    deletions.append(len(display_cursor))
                    display_cursor.append(cursor[jnd])

        if not display_cursor:
            display_cursor = filtered_cursor

        display_results(
            display_cursor,
            additions=additions,
            deletions=deletions,
            sort=True,
            use_source=True,
            **kwargs,
        )

    print("Filtered {} down to {}".format(len(cursor), len(uniq_inds)))

    return filtered_cursor


def filter_cursor(cursor, key, vals, verbosity=0):
    """Returns a cursor obeying the filter on the given key. Any
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
            print("Filtering {} <= {} < {}".format(min_val, key, max_val))
        for doc in cursor:
            try:
                if doc[key] < max_val and doc[key] >= min_val:
                    filtered_cursor.append(doc)
            except (TypeError, ValueError, KeyError):
                pass
    else:
        min_val = float(vals[0])
        if verbosity > 0:
            print("Filtering {} >= {}".format(key, min_val))
        for doc in cursor:
            try:
                if doc[key] >= min_val:
                    filtered_cursor.append(doc)
            except (TypeError, ValueError, KeyError):
                pass

    if verbosity > 0:
        print(orig_cursor_len, "filtered to", len(filtered_cursor), "documents.")
    return filtered_cursor


def filter_cursor_by_chempots(species, cursor):
    """For the desired chemical potentials, remove any incompatible structures
    from cursor.

    Parameters:
        species (list): list of chemical potential formulae.
        cursor (list): list of matador documents to filter.

    Returns:
        list: the filtered cursor.

    """
    from matador.utils.chem_utils import get_number_of_chempots

    # filter out structures with any elements with missing chem pots
    chempot_stoichiometries = []
    for label in species:
        chempot_stoichiometries.append(get_stoich_from_formula(label))

    inds_to_remove = set()
    for ind, doc in enumerate(cursor):
        try:
            cursor[ind]["num_chempots"] = get_number_of_chempots(
                doc, chempot_stoichiometries
            )
        except RuntimeError:
            inds_to_remove.add(ind)
        else:
            cursor[ind]["concentration"] = (
                cursor[ind]["num_chempots"][:-1] / np.sum(cursor[ind]["num_chempots"])
            ).tolist()
            for idx, conc in enumerate(cursor[ind]["concentration"]):
                if conc < 0 + EPS:
                    cursor[ind]["concentration"][idx] = 0.0
                elif conc > 1 - EPS:
                    cursor[ind]["concentration"][idx] = 1.0

    return [doc for ind, doc in enumerate(cursor) if ind not in inds_to_remove]


def _construct_structure_string(
    doc,
    ind,
    formula_substring,
    gs_enthalpy,
    use_source,
    colour,
    hull,
    additions,
    deletions,
    add_index_mode,
    del_index_mode,
    energy_key,
    per_atom,
    eform,
    markdown,
    latex,
):
    """Construct the pretty output for an individual structure.

    Options passed from `matador.utils.cursor_utils.display_results.`

    Returns:
        str: the pretty output.

    """

    # start with two spaces, replaced by the prefix from hull/add/del
    this_struct_string = "  "
    prefix = ""
    suffix = ""
    # apply appropriate prefices and suffices to structure
    if hull and np.abs(doc.get("hull_distance")) <= 0.0 + 1e-12:
        if colour:
            prefix = "\033[92m"
            suffix = "\033[0m"
        this_struct_string = "* "

    if additions is not None:
        if (add_index_mode and ind in additions) or doc.get(
            "text_id", "_"
        ) in additions:
            this_struct_string = "+ "
            if colour:
                prefix = "\033[92m"
                suffix = "\033[0m"
    if deletions is not None:
        if (del_index_mode and ind in deletions) or doc.get(
            "text_id", "_"
        ) in deletions:
            this_struct_string = "- "
            if colour:
                prefix = "\033[91m"
                suffix = "\033[0m"

    # display the canonical name for the structure
    if use_source:
        src = get_root_source(doc["source"])
        max_len = 34
        this_struct_string += "{:<36.{max_len}}".format(
            src if len(src) <= max_len else src[: max_len - 4] + "[..]", max_len=max_len
        )
    else:
        this_struct_string += "{:^24.22}".format(
            " ".join(doc.get("text_id", ["xxx", "yyy"]))
        )

    # again, if we're not outputting to markdown, then flag warnings in the quality column
    try:
        if doc.get("prototype"):
            this_struct_string += "{:^5}".format("*p*")
        elif doc.get("quality", 5) == 0:
            this_struct_string += "{:^5}".format("!!!")
        else:
            this_struct_string += "{:^5}".format((5 - doc.get("quality", 5)) * "?")
    except KeyError:
        this_struct_string += "{:^5}".format(" ")

    # loop over header names and print the appropriate values
    if "pressure" in doc and doc["pressure"] != "xxx":
        this_struct_string += "{: >9.2f} ".format(doc["pressure"])
    else:
        this_struct_string += "{:^9} ".format("xxx")

    try:
        if per_atom and "cell_volume" in doc and "num_atoms" in doc:
            this_struct_string += "{:>12.1f}  ".format(
                doc["cell_volume"] / doc["num_atoms"]
            )
        elif "cell_volume" in doc and "num_fu" in doc:
            this_struct_string += "{:>12.1f}  ".format(
                doc["cell_volume"] / doc["num_fu"]
            )
        else:
            this_struct_string += "{:^12}  ".format("xxx")
    except Exception:
        this_struct_string += "{:^10} ".format("xxx")

    try:
        if hull and eform:
            this_struct_string += "{:>12.3f}      ".format(
                doc["formation_" + energy_key]
            )
        elif hull:
            this_struct_string += "{:>12.1f}      ".format(1000 * doc["hull_distance"])
        elif per_atom:
            this_struct_string += "{:>16.4f}  ".format(
                recursive_get(doc, energy_key) - gs_enthalpy
            )
        else:
            this_struct_string += "{:>16.4f}  ".format(
                recursive_get(doc, energy_key) * doc["num_atoms"] / doc["num_fu"]
                - gs_enthalpy
            )
    except KeyError:
        this_struct_string += "{:^18}".format("xxx")

    if latex:
        from matador.utils.cell_utils import get_space_group_label_latex

        this_struct_string += " {:^13} ".format(
            get_space_group_label_latex(doc.get("space_group", "xxx"))
        )
    else:
        this_struct_string += " {:^13} ".format(doc.get("space_group", "xxx"))

    # now we add the formula column
    this_struct_string += " {:^13} ".format(formula_substring)

    if "num_fu" in doc:
        this_struct_string += " {:^6} ".format(int(doc["num_fu"]))
    else:
        this_struct_string += " {:^6} ".format("xxx")

    if "source" in doc:
        prov = get_guess_doc_provenance(doc["source"], doc.get("icsd"))
        this_struct_string += "{:^8}".format(prov)
    else:
        this_struct_string += "{:^8}".format("xxx")

    this_struct_string = prefix + this_struct_string + suffix

    return this_struct_string


def _construct_source_string(sources):
    """From a list of sources, return a fancy string output
    displaying them as a list.

    """
    num_sources = len(sources)
    if num_sources == 1:
        this_source_string = 11 * " " + "└──────────────────"
    else:
        this_source_string = 11 * " " + "└───────────────┬──"

    for num, _file in enumerate(sources):
        if num_sources == 1:
            this_source_string += ""
        elif num == num_sources - 1:
            this_source_string += (len("└────────────── ") + 11) * " " + "└──"
        elif num != 0:
            this_source_string += (len("└────────────── ") + 11) * " " + "├──"

        this_source_string += " " + _file
        if num != num_sources - 1:
            this_source_string += "\n"

    return this_source_string


def _construct_detail_strings(doc, padding_length=0, source=False):
    """From a document, return a fancy string output
    displaying the desired field names and details of
    the structure.

    Parameters:
        doc (dict): matador document to print.

    Keyword arguments:
        source (bool): whether to allow adjust output so it
            can be interleaved with the source pretty output.
        padding_length (int): how much to pad the detail string
            output.

    Returns:
        (str, str): containing pretty output over two lines.

    """
    detail_string = 11 * " " + "├╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ "
    if source:
        detail_substring = 11 * " " + "├╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ "
    else:
        detail_substring = 11 * " " + "└╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌ "

    if "spin_polarized" in doc:
        if doc["spin_polarized"]:
            detail_string += "S-"

    if "sedc_scheme" in doc:
        detail_string += doc["sedc_scheme"].upper() + "+"

    if "xc_functional" in doc:
        detail_string += doc["xc_functional"]
    else:
        detail_string += "xc-functional unknown!"

    if "cut_off_energy" in doc:
        detail_string += ", {:4.2f} eV".format(doc["cut_off_energy"])
    else:
        detail_string += "cutoff unknown"

    if "external_pressure" in doc:
        detail_string += ", {:4.2f} GPa".format(
            sum(doc["external_pressure"][i][0] for i in range(3)) / 3
        )

    if "kpoints_mp_spacing" in doc:
        detail_string += ", ~{:0.3f} 1/A".format(doc["kpoints_mp_spacing"])

    if "geom_force_tol" in doc:
        detail_string += ", {:.2f} eV/A, ".format(doc["geom_force_tol"])

    if "species_pot" in doc:
        for species in doc["species_pot"]:
            detail_substring += "{}: {}, ".format(species, doc["species_pot"][species])

    if "icsd" in doc:
        detail_substring += "ICSD-CollCode {}, ".format(doc["icsd"])

    if "tags" in doc:
        if isinstance(doc["tags"], list):
            detail_substring += ", ".join(doc["tags"])

    if "user" in doc:
        detail_substring += doc["user"]

    if "encapsulated" in doc and all(
        key in doc for key in ["cnt_chiral", "cnt_radius", "cnt_length"]
    ):
        detail_string += (
            ", (n,m)=("
            + str(doc["cnt_chiral"][0])
            + ","
            + str(doc["cnt_chiral"][1])
            + ")"
        )
        detail_string += ", r={:4.2f} A".format(doc["cnt_radius"])
        detail_string += ", z={:4.2f} A".format(doc["cnt_length"])

    detail_string += " " + (padding_length - len(detail_string) - 1) * "╌"
    detail_substring += " " + (padding_length - len(detail_substring) - 1) * "╌"

    return detail_string, detail_substring


def _construct_header_string(
    markdown, use_source, per_atom, eform, hull, summary, energy_key
):
    """Construct the header of the table, from the passed options.
    For arguments, see docstring for `matador.utils.cursor_utils.display_results`.

    Returns:
        (str, str): the header and units string.

    """
    header_string = ""
    units_string = ""
    if not markdown:
        if use_source:
            header_string += "{:^38}".format("Source")
            units_string += "{:^38}".format("")
        else:
            header_string += "{:^26}".format("ID")
            units_string += "{:^26}".format("")

        header_string += "{:^5}".format("!?!")
        units_string += "{:^5}".format("")
    else:
        header_string += "```\n"
        header_string += "{:^43}".format("Root")
        units_string += "{:^43}".format("")

    header_string += "{:^10}".format("Pressure")
    units_string += "{:^10}".format("(GPa)")

    header_string += "{:^14}".format("Cell volume")
    if per_atom:
        units_string += "{:^14}".format("(Å³/atom)")
    else:
        units_string += "{:^14}".format("(Å³/fu)")

    if eform:
        header_string += "{:^18}".format("Formation energy")
        units_string += "{:^18}".format("(eV/atom)")
    elif hull:
        header_string += "{:^18}".format("Hull dist.")
        units_string += "{:^18}".format("(meV/atom)")
    elif per_atom:
        header_string += "{:^18}".format(
            " ".join(energy_key.replace("_per_atom", "").split("_")).title()
        )
        units_string += "{:^18}".format("(eV/atom)")
    else:
        header_string += "{:^18}".format(
            " ".join(energy_key.replace("_per_atom", "").split("_")).title()
        )
        units_string += "{:^18}".format("(eV/fu)")

    header_string += "{:^15}".format("Space group")
    header_string += " {:^13} ".format("Formula")
    header_string += "{:^8}".format("# fu")
    header_string += "{:^8}".format("Prov.")

    if summary:
        header_string += "{:^12}".format("Occurrences")

    return header_string, units_string


def index_cursors_by_structure(cursors, structure_labeller=get_root_source):
    """For a dictionary of lists of structures, reindex the
    list by the root source of each structure.

    Args:
        cursors: A dictionary of input cursors. Keys will be used
            as labels in the output dictionary.
        structure_labeller: A function called on each structure,
            the result of which will be that structure's key in
            the output dictionary.

    Returns:
        A dictionary with one key per structure, with subkeys
        corresponding to the elements of the initial cursors,
        under which structures are stored from each cursor.

    """
    from collections import defaultdict

    structure_map = defaultdict(dict)
    for label in cursors:
        for s in cursors[label]:
            structure_map[structure_labeller(s)][label] = s

    return structure_map


def _compare_field(bench, other, field):
    """For a given field, compute absolute and relative differences
    between the benchmark and other structure.

    Args:
        bench: The structure to compare against.
        other: The structure to compare.
        field: The field to compare (will be accessed recursively if iterable,
            e.g., `(lattice_abc, 0, 0)`).

    Returns:
        A dictionary summarising the differences.

    """
    if isinstance(field, str):
        field = [field]

    field_label = "_".join(str(_) for _ in field)

    try:
        benchmark_field = recursive_get(bench, field)
    except KeyError:
        raise KeyError(f"Benchmark structure is missing field {field}")

    try:
        other_field = recursive_get(other, field)
    except KeyError:
        raise KeyError(f"Trial structure is missing field {field}")
    # Normalize cell volume to per-atom so different settings can be compared
    if "cell_volume" in field:
        benchmark_field /= bench["num_atoms"]
        other_field /= other["num_atoms"]
    summary = {f"abs_{field_label}": benchmark_field - other_field}
    if abs(benchmark_field) > 1e-10:
        summary[f"rel_{field_label}"] = summary[f"abs_{field_label}"] / benchmark_field

    summary[field_label] = other_field

    return summary


def compare_structures(structures, order, fields=None):
    """Compare structures across various specified or default fields.

    Intended use is to compare crystal structures/energies of the "same"
    crystal when relaxed with different parameters.

    Args:
        structures: A dictionary containing the structures to compare. Keys
            will be used to label the output.
        order: The order of the input keys to use, the first of which will be
            treated as the 'benchmark' structure.
        fields: A list of fields to compare. If None, defaults to comparing
            the lattice parameters, cell volumes and stabilities (hull distance,
            formation energy).

    Returns:
        A dictionary summarising the differences.

    """
    root_sources = set(get_root_source(structures[s]) for s in (structures))
    if len(root_sources) != 1:
        raise RuntimeError(
            f"Not comparing structures with multiple root sources: {root_sources}"
        )

    if fields is None:
        fields = [
            "cell_volume",
            "formation_enthalpy_per_atom",
            "hull_distance",
            ("lattice_abc", 0, 0),
            ("lattice_abc", 0, 1),
            ("lattice_abc", 0, 2),
            ("lattice_abc", 1, 0),
            ("lattice_abc", 1, 1),
            ("lattice_abc", 1, 2),
        ]

    if order[0] not in structures:
        raise RuntimeError(
            f"Benchmark parameter set {order[0]} not found in entry {root_sources}"
        )
    benchmark = structures[order[0]]
    results = {}
    for label in order[1:]:
        summary = {}
        if label in structures:
            for field in fields:
                summary.update(_compare_field(benchmark, structures[label], field))
            results[label] = summary

    return results


def compare_structure_cursor(cursor, order, fields=None):
    """Compare the "same" structures across different accuracies.

    Args:
        cursor: A dict of dicts keyed by structure ID storing data for each
            structure at different accuracies.
        order: An ordered list of the subkeys for each structure; the first
            will be used as the benchmark.
        fields: A list of fields to compare. If None, defaults to comparing
            the lattice parameters, cell volumes and stabilities (hull distance,
            formation energy).

    Returns:
        A dictionary of dictionaries summarising the differences.

    """
    import warnings

    structure_comparator = {}
    for entry in cursor:
        structures = cursor[entry]
        if len(structures) > 1:
            if order[0] not in structures:
                warnings.warn(
                    f"Benchmark parameter set {order[0]} not found for entry {entry}"
                )
                continue
            try:
                structure_comparator[entry] = compare_structures(
                    structures, order=order, fields=fields
                )
            except RuntimeError as exc:
                structure_comparator[entry] = exc

    return structure_comparator
