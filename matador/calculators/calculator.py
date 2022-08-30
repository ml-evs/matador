# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements the abstract Calculator class for use
with compute.

"""
import abc
from matador.utils.errors import CalculationError


class Calculator:
    """Base class for Calculator objects, with some useful methods and
    abstract methods.

    """

    def __init__(self, calculation_parameters):
        self._initial_calculation_parameters = {}
        self._task = calculation_parameters["task"]

    @staticmethod
    @abc.abstractmethod
    def verify_calculation_parameters(calculation_parameters, structure):
        """Verify if the parameters passed make sense for this calculator.
        Mutates the input dict if required.

        Parameters:
            calculation_parameters (dict): contains global parameters, i.e.
                non-structure-specific parameters.
            structure (dict): contains the structure.

        Raises:
            InputError: if verification fails, with a (hopefully) helpful message.

        """
        return

    @staticmethod
    def verify_simulation_cell(structure):
        """Run some checks on the structure.

        Parameters:
            structure (dict): dictionary containing the structure.

        Raises:
            CalculationError: if cell is pathological.

        """
        errors = []
        if all([angle < 30 for angle in structure["lattice_abc"][1]]):
            msg = "Cell is pathological (at least one angle < 30)."
            errors.append(msg)

        # TODO: check overlapping atoms
        if errors:
            raise CalculationError(". ".join(errors))

    def do_memcheck(self):
        """If possible, do a dryrun of the code to estimate memory usage.

        Returns:
            int: estimated memory usage in MB. 0 if method not implemented.

        """
