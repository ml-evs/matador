---
title: 'matador: a Python library for analysing, curating and performing high-throughput density-functional theory calculations'
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
 - name: Theory of Condensed Matter Group, Cavendish Laboratory, University of Cambridge, J. J. Thomson Avenue, Cambridge, CB3 0HE, U.K.
   index: 1
 - name: School of Metallurgy and Materials, University of Birmingham, Edgbaston, Birmingham, B15 2TT, U.K.
   index: 2

date: July 2020
bibliography: paper.bib
---

# Summary

The properties of materials depend heavily on their atomistic structure; knowledge of the possible stable atomic configurations that define a material is required to understand the performance of many technologically and ecologically relevant devices, such as those used for energy storage [@NaP; @CuP]. First-principles crystal structure prediction (CSP) is the art of finding these stable configurations using only quantum mechanics [@JM]. Density-functional theory (DFT) is a ubiquitous theoretical framework for finding approximate solutions to quantum mechanics; calculations using a modern DFT package are sufficiently robust and accurate that insight into real materials can be readily obtained. The computationally intensive work is performed by well-established, low-level software packages, such as CASTEP [@castep] or Quantum Espresso [@qe], which are able to make use of modern high-performance computers. In order to use these codes easily, reliably and reproducibly, many high-level libraries have been developed to create, curate and manipulate the calculations from these low-level workhorses; `matador` is one such framework.

# Statement of need

The purpose of `matador` is fourfold:
- to promote the use of local databases and high-throughput workflows to increase the reproducibility of the computational results, 
- to perform reliable analysis of the stability, structure and properties of materials derived from calculations, 
- to provide tools to create customisable, publication-quality plots of phase diagrams, spectral properties and electrochemistry,
- to make the above functionality available to those with limited programming experience.

# `matador`

`matador` is a Python 3.6+ library and set of command-line tools for performing and analysing high-throughput DFT calculations using the CASTEP [@castep] and Quantum Espresso [@qe] packages. It is well-tested and fully-documented at [ReadTheDocs](https://matador-db.readthedocs.io), and comes with several tutorials and examples. The package is available on PyPI under the name [`matador-db`](https://pypi.org/project/matador-db). As with many projects, `matador` is built on top of the scientific Python ecosystem of NumPy [@numpy; @numpybook], SciPy [@scipy] and matplotlib [@matplotlib].

`matador` has been developed with high-throughput CSP in mind and has found use in the application of CSP to energy storage materials [@NaP; @CuP]; in this use case, a single compositional phase diagram can consist of tens of thousands of structural relaxation calculations. This package is aimed at users of CASTEP or Quantum Espresso who are comfortable with the command-line, yet maybe lack the Python knowledge required to start from scratch with more sophisticated packages. There are many mature packages that provide overlapping functionality with `matador`, the most widespread of which being the Atomic Simulation Environment (ASE) [@ase] and pymatgen [@pymatgen]. A translation layer to and from the structure representation of both of these packages is provided, such that analysis can be reused and combined. 

![Li-Zn-P ternary phase diagram created with matador, plot generated with matplotlib [@matplotlib] and python-ternary [@ternary].\label{fig:hull}](./LiZnP_hull.pdf)

# Overview of functionality

There are two ways of working with `matador`, either from the command-line interface (CLI) or through the Python library directly, with some features that are unique to each. The functionality of `matador` can be broadly split into three categories: 

1. *Creation and curation of databases of the results of first-principles calculations.*

`matador` allows for the creation of [MongoDB](https://mongodb.com) databases of CASTEP (6.0+) geometry optimisations from the command-line, using `matador import`. Only calculations that are deemed "successful", and that have fully-specified inputs are stored, with errors displayed for the rest. The resulting database can be queried with `matador query`, either with Python or through the powerful CLI. The results can be filtered for structural "uniqueness" and written to one of several supported file types and exported for use in other frameworks, such as ASE or pymatgen. Prototyping of structures directly from the database is achieved using `matador swaps`, which uses the same interface as `matador query` to return structure files with "swapped" elements [@NaP].

2. *High-throughput calculations and automated workflows.*

The `run3` executable bundled with `matador` allows for high-throughput calculations to be performed with little setup and no programming knowledge. Specialised support for CASTEP and the post-processing tool OptaDOS [@optados; @optados2] is provided to perform high-throughput geometry optimisations, orbital-projected band structures and densities of states, phonon calculations and elastic properties, however `run3` can also be used to run generic MPI programs concurrently on a set of structures. Sensible defaults for these workflows are provided by leveraging the open-source SeeK-path [@seekpath] and spglib [@spglib] libraries. The bundled `dispersion` script and associated library functionality allows for the creation of publication-quality spectral and vibrational property plots, in a similar fashion to the `sumo` package [@sumo]. The `matador.compute` module behind `run3` also powers the `ilustrado` genetic algorithm code [@ilustrado].

3. *Stability and structural analysis (with an emphasis on battery materials).*

The construction of reliable compositional phase diagrams requires several independent calculations to be performed on different atomic configurations with a compatible set of external parameters. These can be generated from a database query using `matador hull`, which allows the user to filter between different sets of calculations, and, where relevant, `matador voltage` can provide the electrochemical properties of that same phase diagram. Structural fingerprints implemented include pair distribution functions, powder X-ray diffraction patterns, and periodic crystal bond graphs. As more calculations are performed, changes to phase diagrams stored in the local database can be tracked with `matador hulldiff`. Phase diagrams can also be constructed from multiple energy values per structure, for example to show the effects of finite temperature [@CuP], or in the specific case of ensemble-based exchange-correlation functionals like the Bayesian Error Estimate Functional (BEEF) [@beef]. An example of a ternary phase diagram is shown \autoref{fig:hull}.

# Acknowledgements

We acknowledge all the contributors, users and testers of this package, primarily Angela Harper, James Darby, Jordan Dorrell and Matthew Cliffe. M.E. would like to acknowledge the EPSRC Centre for Doctoral Training in Computational Methods for Materials Science for funding under grant number EP/L015552/1. A.J.M. acknowledges funding from EPSRC (EP/P003532/1). The authors acknowledge networking support via the EPSRC Collaborative Computational Projects, CCP9 (EP/M022595/1) and CCP-NC (EP/T026642/1). Much of the development and testing was performed on the Cambridge Service for Data Driven Discovery (CSD3) operated by the University of Cambridge Research Computing Service (http://www.csd3.cam.ac.uk/), provided by Dell EMC and Intel using Tier-2 funding from the Engineering and Physical Sciences Research Council, and DiRAC funding from the Science and Technology Facilities Council (www.dirac.ac.uk).
