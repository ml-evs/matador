# coding: utf-8
# Distributed under the terms of the MIT License.

""" The scrapers module implements the construction of Python dicts from
many common DFT/crystallography file formats.

Currently supported filetypes:

    - CASTEP filetypes: .cell, .param, .castep, .bands, .phonon, .usp
    - Custom .res type based on SHELX files, first used by AIRSS package https://www.mtg.msm.cam.ac.uk/Codes/AIRSS
    - OptaDOS output: .odo, .adaptive.dat
    - Quantum Espresso output files: .out
    - Magres files: .magres
    - Crystallographic Information File (via ASE): .cif

"""


__all__ = ['res2dict', 'cell2dict', 'param2dict', 'castep2dict', 'bands2dict', 'arbitrary2dict',
           'phonon2dict', 'phonon_dos2dict', 'optados2dict', 'usp2dict', 'pwout2dict', 'magres2dict', 'cif2dict']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'

from matador.scrapers.castep_scrapers import res2dict, cell2dict, param2dict
from matador.scrapers.castep_scrapers import castep2dict, bands2dict, arbitrary2dict
from matador.scrapers.castep_scrapers import phonon2dict, phonon_dos2dict, optados2dict, usp2dict
from matador.scrapers.qe_scrapers import pwout2dict
from matador.scrapers.magres_scrapers import magres2dict
from matador.scrapers.cif_scraper import cif2dict
