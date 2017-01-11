# coding: utf-8
""" This file implements some prototype scrapers of
experimental .synth and .expt files.
"""

from __future__ import print_function
# external libraries
import bson.json_util as json
# standard library
from collections import defaultdict


def synth2dict(seed, **kwargs):
    """ Take a .synth file and create a dict
    suitable for database insertion.
    """
    synth = defaultdict(list)
    if seed.endswith('.synth'):
        seed = seed.replace('.synth', '')
    with open(seed+'.synth', 'r') as f:
        flines = f.readlines()
    synth['source'].append(seed+'.synth')
    # define splitters
    splitters = [':', '=']
    try:
        for line_no, line in enumerate(flines):
            # skip blank lines and comments
            if line.startswith(('#', '!')) or len(line.strip()) == 0:
                continue
            else:
                # check if line is split further by semi-colons
                if ';' in line:
                    line = line.split(';')
                    for ind, entry in enumerate(line):
                        line[ind] = line[ind].strip()
                # else set all lines to be single item lists so that
                # we can scrape lines and ;-terminated lines in the same way
                else:
                    line = [line]
                # now look for splitters between tags and data
                for entry in line:
                    for splitter in splitters:
                        if splitter in entry:
                            synth[entry.split(splitter)[0].strip()] = entry.split(splitter)[-1].replace('\n', '')
    except Exception as oopsy:
        if kwargs.get('verbosity') > 0:
            print(oopsy)
            print('Error in', seed+'.synth, skipping...')
        return seed + '\t\t' + str(oopsy), False
    if kwargs.get('debug'):
        print(json.dumps(synth, indent=2, ensure_ascii=False))
    return synth, True


def expt2dict(seed, **kwargs):
    """ Take a .expt file and create a dict
    suitable for database insertion.
    """
    expt = defaultdict(list)
    if seed.endswith('.expt'):
        seed = seed.replace('.expt', '')
    with open(seed+'.expt', 'r') as f:
        flines = f.readlines()
    expt['source'].append(seed+'.expt')
    # define splitters
    splitters = [':', '=']
    try:
        for line_no, line in enumerate(flines):
            # skip blank lines and comments
            if line.startswith(('#', '!')) or len(line.strip()) == 0:
                continue
            else:
                # check if line is split further by semi-colons
                if ';' in line:
                    line = line.split(';')
                    for ind, entry in enumerate(line):
                        line[ind] = line[ind].strip()
                # else set all lines to be single item lists so that
                # we can scrape lines and ;-terminated lines in the same way
                else:
                    line = [line]
                # now look for splitters between tags and data
                for entry in line:
                    for splitter in splitters:
                        if splitter in entry:
                            if entry.split(splitter)[0].strip() in expt:
                                expt[entry.split(splitter)[0].strip()] = [expt[entry.split(splitter)[0].strip()]]
                                expt[entry.split(splitter)[0].strip()].append(entry.split(splitter)[-1].replace('\n', ''))
                            else:
                                expt[entry.split(splitter)[0].strip()] = entry.split(splitter)[-1].replace('\n', '')
    except Exception as oopsy:
        if kwargs.get('verbosity') > 0:
            print(oopsy)
            print('Error in', seed+'.expt, skipping...')
        return seed + '\t\t' + str(oopsy), False
    if kwargs.get('debug'):
        print(json.dumps(expt, indent=2, ensure_ascii=False))
    return expt, True
