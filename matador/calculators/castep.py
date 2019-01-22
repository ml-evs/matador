# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the CastepCalculator class for use
with the compute module.

"""

import os
from matador.calculators import Calculator
from matador.compute.errors import InputError

VALID_PSPOT_LIBS = ['C7', 'C8', 'C9', 'C17', 'C18', 'MS', 'HARD',
                    'QC5', 'NCP', 'NCP18', 'NCP17', 'NCP9']


class CastepCalculator(Calculator):
    """ The CASTEP calculator. """
    @staticmethod
    def verify_calculation_parameters(calculation_parameters, res_dict):
        errors = []
        if 'species_pot' not in calculation_parameters:
            calculation_parameters['species_pot'] = {'library': 'C18'}
        if 'library' not in calculation_parameters['species_pot']:
            for elem in res_dict['stoichiometry']:
                if ('|' not in calculation_parameters['species_pot'][elem[0]] and
                        not os.path.isfile(os.path.expanduser(calculation_parameters['species_pot'][elem[0]])) and
                        calculation_parameters['species_pot'][elem[0]] not in VALID_PSPOT_LIBS):
                    msg = 'Unable to find pseudopotential file/string/library: {}'.format(calculation_parameters['species_pot'][elem[0]])
                    errors.append(msg)

        if 'cut_off_energy' not in calculation_parameters:
            msg = 'Unable to find cut_off_energy field in param file'
            errors.append(msg)
        if 'xc_functional' not in calculation_parameters:
            msg = 'Unable to find xc_functional field in param file'
            errors.append(msg)
        if not any([string in calculation_parameters for string in ['kpoints_list', 'kpoints_mp_grid', 'kpoints_mp_spacing']]):
            msg = 'Unable to find kpoint specifications in cell'
            errors.append(msg)

        if errors:
            raise InputError('. '.join(errors))

    @staticmethod
    def verify_simulation_cell(res_dict):
        super().verify_simulation_cell()
        # any extra checks here
        return

    def do_memcheck(self):
        return 0
