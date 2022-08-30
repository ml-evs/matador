""" This submodule contains aliases for some useful constants. """

from scipy.constants import physical_constants

FARADAY_CONSTANT_Cpermol = physical_constants["Faraday constant"][0]
HARTREE_TO_EV = physical_constants["Hartree energy in eV"][0]
Cperg_to_mAhperg = 2.778e-1
C_TO_mAh = Cperg_to_mAhperg
BOHR_TO_ANGSTROM = physical_constants["Bohr radius"][0] * 1e10
RY_TO_EV = physical_constants["Rydberg constant times hc in eV"][0]
KBAR_TO_GPA = 0.1
eV_PER_ANGSTROM_CUBED_TO_GPa = 160.21776
AVOGADROS_NUMBER = physical_constants["Avogadro constant"][0]
ANGSTROM_CUBED_TO_CENTIMETRE_CUBED = 1e-24
ELECTRON_CHARGE = physical_constants["elementary charge"][0]
KELVIN_TO_EV = physical_constants["kelvin-electron volt relationship"][0]
INVERSE_CM_TO_EV = (
    physical_constants["inverse meter-electron volt relationship"][0] * 100
)
EFG_AU_TO_SI = physical_constants["atomic unit of electric field gradient"][0]
BARN_TO_M2 = 1e-28
PLANCK_CONSTANT = physical_constants["Planck constant"][0]


__all__ = [
    "FARADAY_CONSTANT_Cpermol",
    "HARTREE_TO_EV",
    "Cperg_to_mAhperg",
    "C_TO_mAh",
    "BOHR_TO_ANGSTROM",
    "RY_TO_EV",
    "KBAR_TO_GPA",
    "eV_PER_ANGSTROM_CUBED_TO_GPa",
    "AVOGADROS_NUMBER",
    "ANGSTROM_CUBED_TO_CENTIMETRE_CUBED",
    "ELECTRON_CHARGE",
    "KELVIN_TO_EV",
    "INVERSE_CM_TO_EV",
]
