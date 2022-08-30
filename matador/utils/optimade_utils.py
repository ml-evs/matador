"""This submodule implements some convenience functions for working with OPTIMADE structures."""
from typing import Dict, Any, List, Union
from matador.crystal import Crystal


def optimade2dict_from_url(url: str) -> Union[Crystal, List[Crystal]]:
    """Queries the provided OPTIMADE URL and returns
    a Crystal of list of crystals for the corresponding
    structures.

    Parameters:
        url: The URL of a single or multiple OPTIMADE structure entries.

    Returns:
        The crystal or list of crystals.

    """
    import requests

    data = requests.get(url).json()["data"]
    if isinstance(data, list):
        return [optimade2dict(data) for data in data]

    return optimade2dict(data)


def optimade2dict(structure: Dict[str, Any]) -> Crystal:
    """This function takes an OPTIMADE structure and converts
    it into a matador Crystal."""
    return Crystal(
        {
            "lattice_cart": structure["attributes"]["lattice_vectors"],
            "positions_abs": structure["attributes"]["cartesian_site_positions"],
            "atom_types": structure["attributes"]["species_at_sites"],
            "source": [structure["id"]],
            "optimade_structure": structure,
        }
    )
