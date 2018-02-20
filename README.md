# **matador**

matador is an aggregator and manipulator of the results first-principles calculations, primarily geometry optimisations.

Written by [Matthew Evans](www.tcm.phy.cam.ac.uk/~me388) me388@cam.ac.uk (2016). [API documentation](www.tcm.phy.cam.ac.uk/~me388/matador/). Jupyter notebook examples can be found in ``/examples``.

### New in current version (0.7b):**
* Ternary voltage curves.
* Similarity/uniqueness filtering with element-projected PDFs.
* Updated compute engine for remote calculations (see ``compute.py`` and new binary ``oddjob``).
* Improved test suite and full pip compatiblity.
* Many bugfixes and usability changes.

### New in version (0.6b):
* Intercalation voltage curves, e.g. ```matador voltage -c Li:SnS2```.
* Ternary phase diagrams with heatmaps for structure prediction sampling, gravimetric capacity and formation enthalpy ```matador hull -c ABC --sampmap --efmap --capmap```.
* Substructural similarity interface with Can Kocer's code, as proposed by [Yang et al., PRB (2014)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.054102).

![lipzn.png](docs/img/lipzn.png)

### Core functionality: ###

1. The scraping of CASTEP/Quantum Espresso output files into flexible Python dictionaries with a sensible pragma via `matador import`.
2. The transferal of said objects into a MongoDB database.
3. Powerful CLI querying of the database, with a focus on energy storage applications using ```matador query```.
4. Calculation and presentation of binary and ternary phase diagrams with ```matador hull```.
5. 0K voltage curves for binary and intercalation electrodes using ```matador voltage```.
6. Atomic species swapping and polishing from previous calculations from ```matador swaps/polish```.
7. Automated high-throughput geometry optimisations, convergence tests and densities of states with CASTEP, with the bundled tool ```run3```. Tested on supercomputing systems Darwin (Cambridge) and ARCHER (UK-wide).
8. Refinement of structural data, powered by spglib via ```matador refine```.
9. Fitting of DFT calculations to experimental PDF data through ```matador pdfffit```.

### Usage: ###
```text
usage: matador [-h] [-v]
               
               {stats,query,import,rebuild,pdffit,hull,voltage,swaps,polish,refine}
               ...

MATerial and Atomic Database Of Refined structures.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show programs version number and exit

subcommands:
  valid sub-commands

  {stats,query,import,rebuild,pdffit,hull,voltage,swaps,polish,refine}
    stats               print some stats about the database.
    query               query and extract structures from the database
    import              import new structures in folder into database
    rebuild             rebuild whole database.
    pdffit              provide experimental .gr file and fit to calculated
                        PDF of structures in query
    hull                create a convex hull from query results (currently
                        limited to binaries)
    voltage             plot a voltage curve from query results (currently
                        limited to binaries)
    swaps               perform atomic swaps on query results
    polish              re-relax a series of structures with new parameters.
    refine              refine database structures

```
