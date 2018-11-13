---
title: 'matador: database software for reproducible *ab initio* high-throughput materials science'
tags:
  - Python
  - density-functional theory
  - ab initio
  - crystal structure prediction
  - materials discovery
  - databases
  - castep
  - quantum espresso
  - mongodb
authors:
  - name: Matthew L. Evans
    orcid:  0000-0002-1182-9098 
    affiliation: 1
  - name: Andrew J. Morris
    orcid: 0000-0001-7453-5698
    affiliation: 2
affiliations:
 - name: TCM Group, Cavendish Laboratory, University of Cambridge
   index: 1
 - name: School of Metallurgy and Materials Science, University of Birmingham
   index: 2

date: November 2018
bibliography: paper.bib
---

# Summary

`matador` is a Python 3.5+ library and set of command-line tools for performing and analysing high-throughput atomistic calculations using the CASTEP [1] or QuantumEspresso [2] codes, amongst others. It promotes the use of MongoDB databases to increase the reproducibility and reliability of the computational results, and provides tools to create publication-quality plots of phase diagrams, spectral properties and electrochemistry. The package is fully-documented at (matador-db.readthedocs.io)[matador-db.readthedocs.io], and comes with several tutorials and examples.

`matador` has been developed with high-throughput crystal structure prediction in mind [5, 6]; in this use case, a single compositional phase diagram can consist of hundreds of thousands of density-functional theory relaxations. The target audience of this package is probably already using CASTEP or QuantumEspresso and is comfortable with the command-line, yet maybe lacks the Python knowledge required to start from scratch with other packages. There are many mature packages that provide overlapping functionality with `matador`, the most prominent of which being ASE, AiiDA and pymatgen (+Fireworks). An interface is provided to ASE, such that workflows can be reused and combined. `matador` is intentionally more restrictive, more focussed and thus hopefully more accessible to those without prior programming knowledge than the above packages. 

# Functionality

There are two ways of working with `matador`, either from the command-line or through the Python API, with some features that are unique to each. The functionality of `matador` can be broadly split into three categories: 

1. *Creation and curation of databases of the results of first-principles calculations.*

`matador` currently allows for the creation of MongoDB [9] databases of CASTEP (6.0+) geometry optimisations from the command-line, using `matador import`. Only calculations that are deemed "successful", and that have fully-specified inputs (e.g. CASTEP OTF pseudopotential string rather than arbitrary filename) are stored, with errors displayed for the rest. The resulting database can be queried with `matador query`, either via the API or through the powerful CLI, which consists of handwritten macros rather than an explicit query language. The results can be filtered for structural "uniqueness" (via pair distribution function overlaps) and written to one of several supported filetypes (`.cif`, `.cell`, QuantumEspresso input format, `.shx/.res`, `.pdb`. Database writes have a weak form of version control and can be viewed and reverted using `matador changes`. Prototyping of structures directly from the database is achieved using `matador swaps`, which uses the same interface as `matador query` to return structure files with "swapped" elements [5].

2. *High-throughput calculations and automated workflows.*

The `run3` script bundled with `matador` allows for high-throughput calculations to be performed with little setup and no programming knowledge. Specialised support for CASTEP (and the related tools OptaDOS [11]) is provided to high-throughput geometry optimisations, projected bandstructures and densities of states, phonon calculations and elastic properties, but `run3` can also be used to run generic MPI programs concurrently on a set of structures. These processes are heavily automated by leveraging the open-source SeeKPath [7] and spglib [8] libraries. The bundled `dispersion` script allows for the creation of publication-quality spectral and vibrational property plots (Figure 1a). The backend behind `run3` also powers the `ilustrado` genetic algorithm code (https://bitbucket.org/ml-evs/ilustrado).

3. *Stability and structural analysis (with an emphasis on battery materials).*

The construction of reliable phase diagrams requires several independent calculations to be performed on different atomic configurations with the same set of external parameters. These can be generated from a database query using `matador hull`, which allows the user to filter between different sets of calculations, and, where relevant, `matador voltage` can provide the electrochemical properties of that same phase diagram. Structural descriptors implemented include pair distribution function fingerprinting and periodic crystal bond graphs. As more calculations are performed, changes to phase diagrams can be tracked with `matador hulldiff`.

Figure: a) Projected dispersion + DOS, b) A ternary hull and voltage curve. Created with matplotlib [10].


# Acknowledgements

M.E. would like to acknowledge the EPSRC Centre for Doctoral Training in Computational Methods for Materials Science for funding under grant number EP/L015552/1. Much of the development and testing was performed on the Cambridge Service for Data Driven Discovery (CSD3) operated by the University of Cambridge Research Computing Service (http://www.csd3.cam.ac.uk/), provided by Dell EMC and Intel using Tier-2 funding from the Engineering and Physical Sciences Research Council, and DiRAC funding from the Science and Technology Facilities Council (www.dirac.ac.uk).

# References

[1] CASTEP
[2] QuantumEspresso
[3] ASE
[3] pymatgen
[4] AiiDA
[5] Na-P paper
[6] Na-Sn paper
[7] spglib
[8] seekpath
[9] MongoDB
[10] matplotlib
[11] OptaDOS
