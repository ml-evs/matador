# coding: utf-8
# Distributed under the terms of the MIT License.

""" This script mimics the dispersion.pl script bundled
with CASTEP. For the given bands file, a bandstructure
is created. If a <seed>.adaptive.dat file exists
then a combined BS/DOS plot will be created.
"""

from os.path import isfile
import argparse
import glob
from matador import __version__, script_epilog


def main():
    """Parse args and run the script."""
    parser = argparse.ArgumentParser(
        prog="dispersion",
        epilog=script_epilog,
        description="simple plotting script for bandstructures/DOS based on matador",
        add_help=True,
    )

    parser.add_argument(
        "--version", action="version", version="matador version " + __version__ + "."
    )
    parser.add_argument(
        "--pdf", action="store_true", help="save pdf rather than showing plot in X"
    )
    parser.add_argument(
        "--png", action="store_true", help="save png rather than showing plot in X"
    )
    parser.add_argument(
        "--svg", action="store_true", help="save svg rather than showing plot in X"
    )
    parser.add_argument(
        "--labels", type=str, nargs="*", help="list of legend labels, comma separated"
    )
    parser.add_argument("--dos_only", action="store_true", help="only plot DOS")
    parser.add_argument("--bs_only", action="store_true", help="only plot dispersion")
    parser.add_argument(
        "--preserve_kspace_distance",
        action="store_true",
        help="when linearising kpoint path, ensure distance in reciprocal space is conserved",
    )
    parser.add_argument("--cmap", type=str, help="matplotlib colourmap name to use")
    parser.add_argument(
        "--spin_only", type=str, help='either "up" or "down" to only plot one channel'
    )
    parser.add_argument(
        "-interp",
        "--pdis_interpolation_factor",
        type=float,
        help="multiple by which to interpolate pDIS bands (DEFAULT: 2)",
    )
    parser.add_argument(
        "-scale",
        "--pdis_point_scale",
        type=float,
        help="point scale in pDIS plots (DEFAULT: 25)",
    )
    parser.add_argument(
        "-p",
        "--projectors_to_plot",
        type=str,
        help=(
            "a list of projectors to plot, in the format element:orbital, "
            "e.g. -p K:s,P:p. If the orbital is omitted, all orbitals will be used for that element."
        ),
    )
    parser.add_argument(
        "--unstacked_pdos",
        action="store_true",
        help="plot PDOS as overlap rather than stack",
    )
    parser.add_argument(
        "--no_pdis",
        action="store_true",
        help="do not try to plot PDIS, even if its available",
    )
    parser.add_argument(
        "--band_reorder",
        action="store_true",
        help="try to reorder bands based on local gradients",
    )
    parser.add_argument(
        "--pdos_hide_sum", action="store_true", help="plot PDOS without sum pDOS"
    )
    parser.add_argument(
        "--projector_colours",
        type=str,
        nargs="+",
        help="override all projector colour options with a list of matplotlib-interpretable colours, "
        "e.g. --projector_colours red blue #ff00ff",
    )
    parser.add_argument(
        "--colours",
        type=str,
        nargs="+",
        help="override all colour options with a list of matplotlib-interpretable colours, "
        "e.g. --colours red blue #ff00ff",
    )
    parser.add_argument(
        "-g", "--gap", action="store_true", help="plot position and size of band gap"
    )
    parser.add_argument(
        "-ph",
        "--phonons",
        action="store_true",
        default=False,
        help="plot phonon calculation, rather than electronic",
    )
    parser.add_argument(
        "-ir",
        "--infrared",
        action="store_true",
        default=False,
        help="plot infrared spectrum from file",
    )
    parser.add_argument(
        "-ir_ss",
        "--infrared_step_size",
        type=float,
        nargs="?",
        default=1.0,
        help="step size in cm^{-1} on x-axis for IR plots; must be > 0, default is 1",
    )
    parser.add_argument(
        "-gw",
        "--gaussian_width",
        type=float,
        help=(
            "smearing width for DOS from .bands_dos (default: 0.1 eV) or "
            ".phonon_dos files (default: 10 1/cm)"
        ),
    )
    parser.add_argument(
        "--highlight_bands",
        nargs="+",
        type=int,
        help="specify band numbres to highlight in plot",
    )
    parser.add_argument(
        "-v", "--verbosity", type=int, default=0, help="control verbosity of output"
    )
    parser.add_argument(
        "-figsize",
        "--figsize",
        nargs="+",
        type=int,
        help="figure size in inches to pass to matplotlib",
    )
    parser.add_argument(
        "-pw",
        "--plot_window",
        nargs="+",
        type=float,
        help="energy window [x, y] or [-x, x] to plot either side of E_F (eV)\
                             (DEFAULT: 5 eV)",
    )
    parser.add_argument(
        "seed",
        type=str,
        nargs="+",
        help="seedname or related filename (e.g. bands or dos file)",
    )
    kwargs = vars(parser.parse_args())
    if "gap" in kwargs:
        gap = kwargs["gap"]
        del kwargs["gap"]

    seeds = kwargs.get("seed")
    exts_to_strip = [
        "bands",
        "linear.dat",
        "adaptive.dat",
        "pdis.dat",
        "bands_dos",
        "pdos.dat",
        "phonon",
        "phonon_dos",
    ]

    for ind, seed in enumerate(seeds):
        for ext in exts_to_strip:
            seeds[ind] = seeds[ind].replace("." + ext, "")

    verbosity = kwargs.get("verbosity")
    phonons = kwargs.get("phonons")
    ir = kwargs.get("infrared")
    if kwargs.get("labels"):
        labels = " ".join(kwargs.get("labels")).split(",")
    else:
        labels = None

    if kwargs["no_pdis"]:
        plot_pdis = False
    else:
        plot_pdis = True
    del kwargs["no_pdis"]

    cmap = kwargs.get("cmap")
    band_reorder = kwargs.get("band_reorder", False)

    from matador.plotting import plot_spectral, plot_ir_spectrum

    del kwargs["seed"]
    del kwargs["verbosity"]
    del kwargs["labels"]
    del kwargs["cmap"]
    del kwargs["band_reorder"]

    bandstructure = False
    dos = False

    for seed in seeds:
        if not phonons and not ir:
            bs_seed = seed + ".bands"
            bandstructure = isfile(bs_seed)
            dos_seeds = glob.glob(seed + "*.dat")
            if isfile(seed + ".bands_dos"):
                dos_seeds.append(seed + ".bands_dos")
            dos = any([isfile(dos_seed) for dos_seed in dos_seeds])

        elif phonons and not ir:
            phonon_seed = seed + ".phonon"
            bandstructure = isfile(phonon_seed)
            dos_seed = seed + ".phonon_dos"
            dos = isfile(dos_seed)

        elif ir:
            if len(seeds) > 1:
                exit("Multiple seeds not supported for IR plot.")
            ir_seed = seed + ".phonon"

    if bandstructure:
        cell_seed = seed + ".cell"
        cell = isfile(cell_seed)
    else:
        cell = False

    if kwargs.get("dos_only") and dos:
        bandstructure = False

    if kwargs.get("bs_only") and bandstructure:
        dos = False

    if not any([bandstructure, dos, ir]):
        raise SystemExit("Unable to find files for seedname {}".format(seed))

    if not ir:
        return plot_spectral(
            seeds,
            plot_bandstructure=bandstructure,
            plot_dos=dos,
            plot_pdis=plot_pdis,
            cell=cell,
            gap=gap,
            verbosity=verbosity,
            labels=labels,
            cmap=cmap,
            band_reorder=band_reorder,
            **kwargs
        )

    if ir:
        return plot_ir_spectrum(
            ir_seed, bin_width=kwargs["infrared_step_size"], **kwargs
        )

    exit("Issue plotting {}: did you specify -ph/-ir appropriately?".format(seeds))


if __name__ == "__main__":
    main()
