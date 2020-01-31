""" This submodule implements some useful classes for manipulating
DOS and dispersion data.

"""

from .dos import VibrationalDOS, ElectronicDOS
from .dispersion import VibrationalDispersion, ElectronicDispersion

__all__ = [
    'VibrationalDOS',
    'ElectronicDOS',
    'VibrationalDispersion',
    'ElectronicDispersion'
]
