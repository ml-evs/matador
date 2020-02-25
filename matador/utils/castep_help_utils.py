#!/usr/bin/env python
# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule is essentially a script to scrape CASTEP help strings
for all possible CASTEP parameters.

"""


def scrape_castep_params(executable):
    """ Scan CASTEP help system for parameter file keywords.

    Parameters:
        executable (str): name of CASTEP executable

    Returns:
        list: list of cell keyword strings
        list: list of cell blocks names
        list: list of param keyword strings
        str: castep version number used

    """
    import string
    import subprocess
    params = set()
    cell = set()
    for letter in string.ascii_lowercase:
        output = str(subprocess.check_output('{} -s {}'.format(executable, letter), shell=True))
        for param in output.split('PARAMETERS')[-1].split('\\n')[1:]:
            try:
                if len(param) > 1:
                    params.add(str(param.split()[0]))
            except Exception:
                pass
        for keyword in output.split('PARAMETERS')[0].split('\\n')[1:]:
            try:
                if len(keyword) > 1:
                    cell.add(str(keyword.split()[0]))
            except Exception:
                pass

    output = str(subprocess.check_output('{} -s block'.format(executable), shell=True))
    blocks = set()
    for keyword in output.split('PARAMETERS')[0].split('\\n')[1:]:
        try:
            if len(keyword) > 1:
                blocks.add(str(keyword.split()[0]))
        except Exception:
            pass

    params = list(sorted(list(params)))
    cell = list(sorted(list(cell)))
    blocks = list(sorted(list(blocks)))
    cell = [keyword for keyword in cell if keyword not in blocks]

    version_string = str(subprocess.check_output('{} --version'
                         .format(executable), shell=True)).split('\\n')[0].split()[-1].strip()
    return cell, blocks, params, version_string


def update_castep_param_list(executable):
    """ Update the castep_params.py file.

    Parameters:
        executable (str): name of CASTEP executable

    """
    cell, blocks, params, version = scrape_castep_params(executable)
    with open('castep_params.py', 'w') as f:
        f.write('""" This file contains a Python list of all CASTEP parameters,\n')
        f.write('automatically generated with file_utils.scrape_castep_params().\n"""')
        f.write('\n\n')
        f.write('CASTEP_VERSION = {}'.format(version))
        f.write('\n\n')
        f.write('CASTEP_CELL_KEYWORDS = [\n')
        for keyword in cell:
            f.write('                 \'{}\',\n'.format(keyword.lower()))
        f.write('                ]\n')

        f.write('CASTEP_CELL_BLOCKS = [\n')
        for keyword in blocks:
            f.write('                 \'{}\',\n'.format(keyword.lower()))
        f.write('                ]\n')

        f.write('CASTEP_PARAMS = [\n')
        for param in params:
            f.write('                 \'{}\',\n'.format(param.lower()))
        f.write('                ]\n')


if __name__ == '__main__':
    from sys import argv
    if len(argv) > 1:
        executable = argv[1]
    else:
        executable = 'castep'
    print('Updating castep_parameters.py with executable {}'.format(executable))
    update_castep_param_list(executable)
