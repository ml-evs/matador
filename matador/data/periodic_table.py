"""This submodule implements some convenience dictionaries for working with the periodic table.

Data was scraped from a combination of the periodictable Python package (on which matador used to
depend) and pymatgen.

"""
from typing import Dict
from dataclasses import dataclass

__all__ = ("PERIODIC_TABLE",)


@dataclass
class Element:
    symbol: str
    number: int
    mass: float
    iupac_order: int


PERIODIC_TABLE: Dict[str, Element] = {
    "n": Element(symbol="n", number=0, mass=0, iupac_order=-1),
    "H": Element(symbol="H", number=1, mass=1.00794, iupac_order=92),
    "He": Element(symbol="He", number=2, mass=4.002602, iupac_order=5),
    "Li": Element(symbol="Li", number=3, mass=6.941, iupac_order=11),
    "Be": Element(symbol="Be", number=4, mass=9.012182, iupac_order=17),
    "B": Element(symbol="B", number=5, mass=10.811, iupac_order=81),
    "C": Element(symbol="C", number=6, mass=12.0107, iupac_order=86),
    "N": Element(symbol="N", number=7, mass=14.0067, iupac_order=91),
    "O": Element(symbol="O", number=8, mass=15.9994, iupac_order=97),
    "F": Element(symbol="F", number=9, mass=18.9984032, iupac_order=102),
    "Ne": Element(symbol="Ne", number=10, mass=20.1797, iupac_order=4),
    "Na": Element(symbol="Na", number=11, mass=22.98977, iupac_order=10),
    "Mg": Element(symbol="Mg", number=12, mass=24.305, iupac_order=16),
    "Al": Element(symbol="Al", number=13, mass=26.981538, iupac_order=80),
    "Si": Element(symbol="Si", number=14, mass=28.0855, iupac_order=85),
    "P": Element(symbol="P", number=15, mass=30.973761, iupac_order=90),
    "S": Element(symbol="S", number=16, mass=32.065, iupac_order=96),
    "Cl": Element(symbol="Cl", number=17, mass=35.453, iupac_order=101),
    "Ar": Element(symbol="Ar", number=18, mass=39.948, iupac_order=3),
    "K": Element(symbol="K", number=19, mass=39.0983, iupac_order=9),
    "Ca": Element(symbol="Ca", number=20, mass=40.078, iupac_order=15),
    "Sc": Element(symbol="Sc", number=21, mass=44.95591, iupac_order=49),
    "Ti": Element(symbol="Ti", number=22, mass=47.867, iupac_order=52),
    "V": Element(symbol="V", number=23, mass=50.9415, iupac_order=55),
    "Cr": Element(symbol="Cr", number=24, mass=51.9961, iupac_order=58),
    "Mn": Element(symbol="Mn", number=25, mass=54.938049, iupac_order=61),
    "Fe": Element(symbol="Fe", number=26, mass=55.845, iupac_order=64),
    "Co": Element(symbol="Co", number=27, mass=58.9332, iupac_order=67),
    "Ni": Element(symbol="Ni", number=28, mass=58.6934, iupac_order=70),
    "Cu": Element(symbol="Cu", number=29, mass=63.546, iupac_order=73),
    "Zn": Element(symbol="Zn", number=30, mass=65.409, iupac_order=76),
    "Ga": Element(symbol="Ga", number=31, mass=69.723, iupac_order=79),
    "Ge": Element(symbol="Ge", number=32, mass=72.64, iupac_order=84),
    "As": Element(symbol="As", number=33, mass=74.9216, iupac_order=89),
    "Se": Element(symbol="Se", number=34, mass=78.96, iupac_order=95),
    "Br": Element(symbol="Br", number=35, mass=79.904, iupac_order=100),
    "Kr": Element(symbol="Kr", number=36, mass=83.798, iupac_order=2),
    "Rb": Element(symbol="Rb", number=37, mass=85.4678, iupac_order=8),
    "Sr": Element(symbol="Sr", number=38, mass=87.62, iupac_order=14),
    "Y": Element(symbol="Y", number=39, mass=88.90585, iupac_order=48),
    "Zr": Element(symbol="Zr", number=40, mass=91.224, iupac_order=51),
    "Nb": Element(symbol="Nb", number=41, mass=92.90638, iupac_order=54),
    "Mo": Element(symbol="Mo", number=42, mass=95.94, iupac_order=57),
    "Tc": Element(symbol="Tc", number=43, mass=98, iupac_order=60),
    "Ru": Element(symbol="Ru", number=44, mass=101.07, iupac_order=63),
    "Rh": Element(symbol="Rh", number=45, mass=102.9055, iupac_order=66),
    "Pd": Element(symbol="Pd", number=46, mass=106.42, iupac_order=69),
    "Ag": Element(symbol="Ag", number=47, mass=107.8682, iupac_order=72),
    "Cd": Element(symbol="Cd", number=48, mass=112.411, iupac_order=75),
    "In": Element(symbol="In", number=49, mass=114.818, iupac_order=78),
    "Sn": Element(symbol="Sn", number=50, mass=118.71, iupac_order=83),
    "Sb": Element(symbol="Sb", number=51, mass=121.76, iupac_order=88),
    "Te": Element(symbol="Te", number=52, mass=127.6, iupac_order=94),
    "I": Element(symbol="I", number=53, mass=126.90447, iupac_order=99),
    "Xe": Element(symbol="Xe", number=54, mass=131.293, iupac_order=1),
    "Cs": Element(symbol="Cs", number=55, mass=132.90545, iupac_order=7),
    "Ba": Element(symbol="Ba", number=56, mass=137.327, iupac_order=13),
    "La": Element(symbol="La", number=57, mass=138.9055, iupac_order=47),
    "Ce": Element(symbol="Ce", number=58, mass=140.116, iupac_order=46),
    "Pr": Element(symbol="Pr", number=59, mass=140.90765, iupac_order=45),
    "Nd": Element(symbol="Nd", number=60, mass=144.24, iupac_order=44),
    "Pm": Element(symbol="Pm", number=61, mass=145, iupac_order=43),
    "Sm": Element(symbol="Sm", number=62, mass=150.36, iupac_order=42),
    "Eu": Element(symbol="Eu", number=63, mass=151.964, iupac_order=41),
    "Gd": Element(symbol="Gd", number=64, mass=157.25, iupac_order=40),
    "Tb": Element(symbol="Tb", number=65, mass=158.92534, iupac_order=39),
    "Dy": Element(symbol="Dy", number=66, mass=162.5, iupac_order=38),
    "Ho": Element(symbol="Ho", number=67, mass=164.93032, iupac_order=37),
    "Er": Element(symbol="Er", number=68, mass=167.259, iupac_order=36),
    "Tm": Element(symbol="Tm", number=69, mass=168.93421, iupac_order=35),
    "Yb": Element(symbol="Yb", number=70, mass=173.04, iupac_order=34),
    "Lu": Element(symbol="Lu", number=71, mass=174.967, iupac_order=33),
    "Hf": Element(symbol="Hf", number=72, mass=178.49, iupac_order=50),
    "Ta": Element(symbol="Ta", number=73, mass=180.9479, iupac_order=53),
    "W": Element(symbol="W", number=74, mass=183.84, iupac_order=56),
    "Re": Element(symbol="Re", number=75, mass=186.207, iupac_order=59),
    "Os": Element(symbol="Os", number=76, mass=190.23, iupac_order=62),
    "Ir": Element(symbol="Ir", number=77, mass=192.217, iupac_order=65),
    "Pt": Element(symbol="Pt", number=78, mass=195.078, iupac_order=68),
    "Au": Element(symbol="Au", number=79, mass=196.96655, iupac_order=71),
    "Hg": Element(symbol="Hg", number=80, mass=200.59, iupac_order=74),
    "Tl": Element(symbol="Tl", number=81, mass=204.3833, iupac_order=77),
    "Pb": Element(symbol="Pb", number=82, mass=207.2, iupac_order=82),
    "Bi": Element(symbol="Bi", number=83, mass=208.98038, iupac_order=87),
    "Po": Element(symbol="Po", number=84, mass=209, iupac_order=93),
    "At": Element(symbol="At", number=85, mass=210, iupac_order=98),
    "Rn": Element(symbol="Rn", number=86, mass=222, iupac_order=0),
    "Fr": Element(symbol="Fr", number=87, mass=223, iupac_order=6),
    "Ra": Element(symbol="Ra", number=88, mass=226, iupac_order=12),
    "Ac": Element(symbol="Ac", number=89, mass=227, iupac_order=32),
    "Th": Element(symbol="Th", number=90, mass=232.0381, iupac_order=31),
    "Pa": Element(symbol="Pa", number=91, mass=231.03588, iupac_order=30),
    "U": Element(symbol="U", number=92, mass=238.02891, iupac_order=29),
}