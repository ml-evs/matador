# coding: utf-8
# Distributed under the terms of the MIT License.

""" The calculator module contains classes for use as DFT/atomistic
calculators within the compute module class.

"""


__all__ = ['Calculator', 'CastepCalculator']  # , 'QuantumEspressoCalculator', 'ASECalculator']
__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.calculators.castep import CastepCalculator
# from matador.calculators.quantum_espresso import QuantumEspressoCalculator
# from matador.calculators.ase import ASECalculator
