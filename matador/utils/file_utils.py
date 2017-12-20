#!/usr/bin/env python
""" This file implements some simple utitlies to do with
parsing/writing of files. """


def scrape_castep_params(executable):
    """ Scan CASTEP help system for parameter file keywords. """
    import string
    import subprocess
    params = set()
    for letter in string.ascii_lowercase:
        output = str(subprocess.check_output('{} -s {}'.format(executable, letter), shell=True))
        for param in output.split('PARAMETERS')[-1].split('\\n')[1:]:
            try:
                if len(param) > 1:
                    params.add(str(param.split()[0]))
            except:
                pass
    params = list(params)
    params = sorted(params)
    version_string = str(subprocess.check_output('{} --version'.format(executable), shell=True)).split('\\n')[0].split()[-1].strip()
    return params, version_string


def update_castep_param_list(executable):
    """ Update the param list file. """
    params, version = scrape_castep_params(executable)
    with open('castep_params.py', 'w') as f:
        f.write('""" This file contains a Python list of all CASTEP parameters,')
        f.write('automatically generated from file_utils.scrape_castep_params.\n"""')
        f.write('\n\n')
        f.write('CASTEP_VERSION = {}'.format(version))
        f.write('\n\n')
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
