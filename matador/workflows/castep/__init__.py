# coding: utf-8
# Distributed under the terms of the MIT License.

""" The workflows.castep submodule contains any workflows that rely on
CASTEP-specific implementation, for example electronic spectropscopy or phonons.

"""


__all__ = ['castep_full_phonon', 'castep_full_spectral', 'castep_elastic']

__author__ = 'Matthew Evans'
__maintainer__ = 'Matthew Evans'


from matador.workflows.castep.phonons import castep_full_phonon
from matador.workflows.castep.spectral import castep_full_spectral
from matador.workflows.castep.elastic import castep_elastic
