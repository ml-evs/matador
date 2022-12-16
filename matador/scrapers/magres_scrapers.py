# coding: utf-8
# Distributed under the terms of the MIT License.

""" This submodule implements some scraper functions for
NMR-related inputs and outputs, e.g. .magres files.

"""


from collections import defaultdict
from os import stat

import numpy as np
from matador.utils.cell_utils import cart2abc, cart2frac
from matador.utils.chem_utils import get_stoich
from matador.scrapers.utils import scraper_function, get_flines_extension_agnostic
from matador.data.constants import (
    ELECTRON_CHARGE,
    PLANCK_CONSTANT,
    BARN_TO_M2,
    EFG_AU_TO_SI,
)
from matador.data.magres import ELECTRIC_QUADRUPOLE_MOMENTS


@scraper_function
def magres2dict(fname, **kwargs):
    """Extract available information from .magres file. Assumes units of
    Angstrom and ppm for relevant quantities.
    """
    magres = defaultdict(list)
    flines, fname = get_flines_extension_agnostic(fname, "magres")
    magres["source"] = [fname]

    # grab file owner username
    try:
        from pwd import getpwuid

        magres["user"] = getpwuid(stat(fname).st_uid).pw_name
    except Exception:
        magres["user"] = "xxx"

    magres["magres_units"] = dict()
    for line_no, line in enumerate(flines):
        line = line.lower().strip()
        if line in ["<atoms>", "[atoms]"]:
            i = 1
            while flines[line_no + i].strip().lower() not in ["</atoms>", "[/atoms]"]:
                split_line = flines[line_no + i].split()
                if not split_line:
                    i += 1
                    continue
                if i > len(flines):
                    raise RuntimeError("Something went wrong in reader loop")
                if split_line[0] == "units":
                    magres["magres_units"][split_line[1]] = split_line[2]
                elif "lattice" in split_line:
                    lattice = split_line[1:]
                    for j in range(3):
                        magres["lattice_cart"].append(
                            [float(elem) for elem in lattice[j * 3 : (j + 1) * 3]]
                        )
                    magres["lattice_abc"] = cart2abc(magres["lattice_cart"])
                elif "atom" in split_line:
                    atom = split_line
                    magres["atom_types"].append(atom[1])
                    magres["positions_abs"].append([float(elem) for elem in atom[-3:]])
                i += 1
            break

    if "atom_types" in magres:
        magres["num_atoms"] = len(magres["atom_types"])
        magres["positions_frac"] = cart2frac(
            magres["lattice_cart"], magres["positions_abs"]
        )
        magres["stoichiometry"] = get_stoich(magres["atom_types"])

    for line_no, line in enumerate(flines):
        line = line.lower().strip()
        if line in ["<magres>", "[magres]"]:
            i = 1
            while flines[line_no + i].strip().lower() not in ["</magres>", "[/magres]"]:
                split_line = flines[line_no + i].split()
                if not split_line:
                    i += 1
                    continue
                if i > len(flines):
                    raise RuntimeError("Something went wrong in reader loop")
                if split_line[0] == "units":
                    magres["magres_units"][split_line[1]] = split_line[2]
                elif "sus" in split_line:
                    magres["susceptibility_tensor"] = np.array(
                        [float(val) for val in split_line[1:]]
                    ).reshape(3, 3)

                elif "ms" in split_line:
                    ms = np.array([float(val) for val in split_line[-9:]]).reshape(3, 3)
                    s_iso = np.trace(ms) / 3

                    # find eigenvalues of symmetric part of shielding and order them to calc anisotropy eta
                    symmetric_shielding = _symmetrise_tensor(ms)
                    s_yy, s_xx, s_zz = _get_haeberlen_eigs(symmetric_shielding)
                    s_aniso = s_zz - (s_xx + s_yy) / 2.0
                    asymm = (s_yy - s_xx) / (s_zz - s_iso)

                    # convert from reduced anistropy to CSA
                    magres["magnetic_shielding_tensors"].append(ms)
                    magres["chemical_shielding_isos"].append(s_iso)
                    magres["chemical_shift_anisos"].append(s_aniso)
                    magres["chemical_shift_asymmetries"].append(asymm)

                elif "efg" in split_line:
                    efg = np.array([float(val) for val in split_line[-9:]]).reshape(
                        3, 3
                    )
                    species = split_line[1]

                    eigs = _get_haeberlen_eigs(efg)
                    v_zz, eta = eigs[2], (eigs[0] - eigs[1]) / eigs[2]

                    # calculate C_Q in MHz
                    quadrupole_moment = ELECTRIC_QUADRUPOLE_MOMENTS.get(species, 1.0)

                    C_Q = (
                        ELECTRON_CHARGE
                        * v_zz
                        * quadrupole_moment
                        * EFG_AU_TO_SI
                        * BARN_TO_M2
                    ) / (PLANCK_CONSTANT * 1e6)

                    magres["electric_field_gradients"].append(efg)
                    magres["quadrupolar_couplings"].append(C_Q)
                    magres["quadrupolar_asymmetries"].append(eta)

                i += 1

    for line_no, line in enumerate(flines):
        line = line.lower().strip()
        if line in ["<calculation>", "[calculation]"]:
            i = 1
            while flines[line_no + i].strip().lower() not in [
                "</calculation>",
                "[/calculation]",
            ]:
                if i > len(flines):
                    raise RuntimeError("Something went wrong in reader loop")
                # space important as it excludes other calc_code_x variables
                if "calc_code " in flines[line_no + i]:
                    magres["calculator"] = flines[line_no + i].split()[1]
                if "calc_code_version" in flines[line_no + i]:
                    magres["calculator_version"] = flines[line_no + i].split()[1]
                i += 1

    return dict(magres), True


def _symmetrise_tensor(t):
    """Returns the symmetrised tensor.

    Arguments:
        t: numpy array containing a NxN square tensor.

    Returns:
        np.ndarray: NxN numpy array containing the symmetric tensor.

    """
    return 0.5 * (t + t.T)


def _get_haeberlen_eigs(t: np.ndarray):
    """Return Haeberlen convention ordered eigenvalues of the passed tensor.

    ``|s_zz - s_iso| >= |s_xx - s_iso| >= |s_yy - s_iso|``

    Arguments:
        t: numpy array containing a NxN square tensor.

    Returns:
        np.ndarray: N-d numpy array containing the ordered eigenvalues.

    """
    eig_vals, eig_vecs = np.linalg.eig(t)
    eig_vals, eig_vecs = zip(
        *sorted(zip(eig_vals, eig_vecs), key=lambda eig: abs(eig[0] - np.trace(t) / 3))
    )

    return eig_vals
