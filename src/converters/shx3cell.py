#!/usr/bin/python
# coding: utf-8
from scrapers import castep_scrapers
from export import doc2cell
from print_utils import print_notify, print_success, print_failure
from sys import argv

fnames = argv[1:]
for fname in fnames:
    print_notify('Reading ' + fname)
    res_dict, success = castep_scrapers.res2dict(fname, db=False)
    if success:
        doc2cell(res_dict, fname.replace('.res', ''), copy_pspots=False, hash_dupe=True)
        print_success('Wrote .cell file to ' + fname.replace('.cell', '') + '.cell')
    else:
        print_failure('Unable to find final structure in ' + fname + '.')
print_success('Completed!')
