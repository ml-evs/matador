""" This submodule implements some useful classes for manipulating
DOS and dispersion data.

"""

from .dos import VibrationalDOS, ElectronicDOS, DensityOfStates
from .dispersion import VibrationalDispersion, ElectronicDispersion, Dispersion


__all__ = [
    'VibrationalDOS',
    'ElectronicDOS',
    'DensityOfStates',
    'VibrationalDispersion',
    'ElectronicDispersion',
    'Dispersion'
]
