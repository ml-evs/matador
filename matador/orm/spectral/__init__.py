""" This submodule implements some useful classes for manipulating
DOS and dispersion data.

"""

from .dos import VibrationalDOS, ElectronicDOS, DensityOfStates
from .dispersion import VibrationalDispersion, ElectronicDispersion, Dispersion
from .spectral import Spectral


__all__ = [
    'VibrationalDOS',
    'ElectronicDOS',
    'DensityOfStates',
    'VibrationalDispersion',
    'ElectronicDispersion',
    'Dispersion',
    'Spectral'
]
