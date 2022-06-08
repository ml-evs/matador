"""This submodule implements some convenience functions for working with OPTIMADE structures."""
from typing import Dict, Any
from matador.crystal import Crystal


def optimade2dict(structure: Dict[str, Any]) -> Crystal:
    """This function takes an OPTIMADE structure and converts
    it into a matador Crystal."""
    return Crystal({
        "lattice_cart": structure["attributes"]["lattice_vectors"],
        "positions_abs": structure["attributes"]["cartesian_site_positions"],
        "atom_types": structure["attributes"]["species_at_sites"],
        "source": [structure["id"]],
        "optimade_structure": structure
    })
