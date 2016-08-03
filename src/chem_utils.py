# coding: utf-8
""" This file defines some useful chemistry. """
import periodictable
import numpy as np

FARADAY_CONSTANT_Cpermol = 96.485332e3
Cperg_to_mAhperg = 2.778e-1


def get_periodic_table():
    """ Return some periodic table macros. """
    periodic_table = dict()
    periodic_table['I'] = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
    periodic_table['II'] = ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']
    periodic_table['III'] = ['B', 'Al', 'Ga', 'In', 'Tl']
    periodic_table['IV'] = ['C', 'Si', 'Ge', 'Sn', 'Pb']
    periodic_table['V'] = ['N', 'P', 'As', 'Sb', 'Bi']
    periodic_table['VI'] = ['O', 'S', 'Se', 'Te', 'Po']
    periodic_table['VII'] = ['F', 'Cl', 'Br', 'I', 'At']
    periodic_table['Tran'] = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                              'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                              'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
    periodic_table['Lan'] = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                             'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    periodic_table['Act'] = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
                             'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
    return periodic_table


def get_molar_mass(elem):
    """ Return molar mass of chosen element. """
    return periodictable.elements.symbol(elem).mass


def get_capacities(x, m_B):
    """ Return capacity in mAh/g from x in A_x B
    and m_B in a.m.u.
    """
    x = np.array(x)
    Q = x * FARADAY_CONSTANT_Cpermol * Cperg_to_mAhperg / m_B
    return Q
