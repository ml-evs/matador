# coding: utf-8
# Distributed under the terms of the MIT License.

""" The workflows.castep submodule contains any workflows that rely on
CASTEP-specific implementation, for example bandstructures or phonons
(but not e.g. bulk modulus calculations, which require just changing
volumes in a non-code specific way).

"""


__all__ = ['castep_full_phonon']

__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.workflows.castep.phonons import castep_full_phonon
from matador.workflows.castep.spectral import castep_full_spectral
