.. _changelog:

Changelog
=========


New in release (0.10.1) [20/01/2023]
------------------------------------

- Support for continuation and restarting of magres workflows (#306)
- Fix for CASTEP scraper to scrape both constrained and unconstrained stress tensors when both are present (#312)


New in release (0.10.0) [26/10/2022]
------------------------------------

This release contains many accumulated changes over the course of a year.
The library is now primarily in maintenance mode unless new features are
contributed by others.

- Added a ray-tracing crystal visualisation module based on the `Fresnel library <https://github.com/glotzerlab/fresnel>`_.
- Refinements of almost all plotting types and styles, as well as making it
  easier to plot multiple comparative plots on the same figure (e.g., multiple
  phase diagrams computed with different functionals or from different sources).
- Added the :func:`matador.magres` module for referencing NMR calculations
  against experimental data.
- Optimised the CASTEP output reader for long concatenated geometry optimisations.
- Minimal integration with `OPTIMADE <https://optimade.org>`_ by being able to
  download an OPTIMADE structure from a URL as a matador object.
- Support for Python 3.10
- Bug fixes for scrapers, densities of states plotting and scattering simulations.


New in release (0.9.11) [03/06/2021]
------------------------------------

- Minor change: allow specification of arbitrary strings for CASTEP pseudopotential libraries (#156).
- Bug fix: ``standardize_cell`` script failing to use default symmetry tolerance value (#157).
- Bug fix: scraping of .cif files with single atoms and no symmetries (#173)
- Bug fix: scraping of Hubbard U values from .castep files, and associated bugs when performing relaxations with Hubbard U (#180)
- Dependency updates and Python 3.6 deprecation warning (#158, #181)

New in release (0.9.10) [23/02/2021]
------------------------------------

- Windows compatibility changes (#149)
- Dependency updates (#146, #148, #149)

New in release (0.9.9) [16/10/2020]
-----------------------------------

- Added support for CASTEP kpoint path ``BREAK`` directive (#107)
- Improvements to magres plotting and magres workflow (#112)
- Added ability to scrape electric field gradient values and compute quadrupolar quantities from NMR calculations (#115)
- Added ability to run all several examples under Binder (#106, #130).
- JOSS paper accepted! (#129)


New in release (0.9.8) [10/08/2020]
-----------------------------------
- Improvements to PDIS functionality (#94).

  - Rasterized scatter points for more efficient exporting and fewer graphical artifacts
  - Made underlying :func:`matador.plotting.spectral_plotting.dos_plot` and :func:`matador.plotting.spectral_plotting.dispersion_plot` more API friendly, and added example notebook.
  - Fixed bug in cell scraping for old ``BS_*`` style keywords.

- Improvements to magres functionality, including scraping of units (#90)
- Example notebooks that do not need external data/databases are now run as part of CI (#91).
- New workflow for NMR calculations and refactoring of old workflows (#96).

  - New workflow performs relaxation and high-quality SCF before NMR calculation.
  - Old workflows refactored and improved to enforce certain required parameters for e.g. checkpointing.
  - Enabled phonon workflow for CASTEP ``PHONON+EFIELD`` task.
  - Made file scrapers less dependent on file type.

- Updated CASTEP parameter list to 20.1 (#97).
- Tweaked spectral plotting defaults, including ``--colours`` flag to dispersion script (#98).


New in release (0.9.7) [29/07/2020]
-----------------------------------
- Bug fixes to problems introduced in 0.9.6.
- Cosmetic fixes to logging and misleading status reports in workflows.


New in release (0.9.6) [28/07/2020]
-----------------------------------
- Improvements to ASE and pymatgen interoperability (#80)
- Fixed bug in :class:`matador.hull.TemperatureDependentHull` which would crash when not provided a list of temperatures (#82).
- Added plotting functions for magres data, and improved its handling inside :class:`matador.crystal.Crystal` (#79).

New in release (0.9.5) [25/06/2020]
-----------------------------------
- This release is mostly to trigger Zenodo archiving.
- Updated README and tests for recent Python versions.


New in release (0.9.4) [08/06/2020]
-----------------------------------
- Fixed flag help strings for ``pxrd_calculator`` (#65)
- Changed default PDF broadening for 3x speedup (#65)
- Reverted ``cpu_count`` to use version that works correctly in most cases, by chance (#66).


New in release (0.9.3) [07/06/2020]
-----------------------------------

- Fixes for the CIF reader: now works with awkward linebreaks and alternative symmetry operation specifications (#61).
- Added several new flags to ``pxrd_calculator`` script (#60 and 61).
- Usability fixes for ``spectral_plotting`` in the case of projected dispersion curves (#59).


New in release (0.9.2) [01/06/2020]
-----------------------------------

- Optimised CIF reader considerably (#50)
- Updated PXRD calculator to allow for partial occupancy, monochromated beam angles and text export, and added ``pxrd_calculator`` script for convenience when handling CIF files.
- Added ability to choose which projectors are plotted with dispersion (#47)
- Various minor fixes and updates:

  - Updates to docs for CLI and configuration.
  - Allow nan-values to be reset inside :class:`matador.crystal.Crystal`.
  - Fixed display ordering of fingerprint-filtered cursors.


New in release (0.9.1) [20/05/2020]
-----------------------------------

- Fixed issue with local pip installs after 0.9 release
- Fixed issue with multi-node MPI tasks by switching to ``proc.communicate()`` after an initial polling stage (#37)
- Fixed issue where bands would be reordered multiple times in spectral plots (#40)
- Tweaked spectral plot defaults (#40)
- Replaced ``multiprocessing.cpu_count()`` calls with ``psutil.cpu_count(logical=False)`` to avoid double-counting hyperthreaded cores


New in release (0.9) [15/05/2020]
---------------------------------

- PyPI release! Can now install with ``pip install matador-db`` (unfortunately ``matador`` was taken, but they are sufficiently orthogonal that the package name ``matador`` is retained here.
- Much improved code structure and many additional classes that wrap the raw calculation dictionaries for e.g. :class:`matador.crystal.Crystal` and spectral classes.
- New module :mod:`matador.orm` containing useful models for data handling.

  - :class:`matador.orm.orm.DataContainer` as a base class that allows for easy
    access to underlying dictionaries.
  - :mod:`matador.orm.spectral` module that contains many useful classes for
    manipulating and plotting e.g. bandstructures, DOS and finite temperature
    properties.

- New features in :mod:`matador.hull` module:

  - Pseudo-ternary phase diagrams (building towards arbitrary n-dimensional phase diagrams).
  - :class:`matador.hull.EnsembleHull` class and submodule to support the Bayesian Error Estimate Functional (BEEF) and finite temperature phase diagrams.
  - Refactoring of hull calculation into light-weight :class:`matador.hull.PhaseDiagram` class.
  - Finite temperature hulls based on :class:`matador.hull.EnsembleHull` with
    :class:`matador.hull.TemperatureDependentHull`.

- Refactored old PDF `similarity` module into new module :mod:`matador.fingerprints`.

  - Added new fingerprint class, :class:`matador.fingerprints.PXRD`, with associated plots (thanks for James Darby for some initial code). Defaults calibrated with GSAS-II.
  - :class:`matador.fingerprints.PDF` sped up by an order of magnitude using `numba`.

- :class:`matador.workflows.castep.CastepSpectralWorkflow` extended to include latest projected dispersion curve developments from OptaDOS, with associated projected dispersion plots (see tutorial).

  - Updated dispersion script to automatically perform naive Gaussian smearing if OptaDOS output not detected.

- Abstracted and simplified :mod:`matador.compute` module to allow for extension to new codes via :mod:`matador.compute.calculators` submodule.

  - Should now be more robust and transferrable, with many HPC environments automatically detected.
  - Added ``--scratch_prefix`` to run3 to allow for temporary files to e.g. be written to faster filesystem with appropriate symlinks to work folder.

- All CASTEP 19 keywords supported, as well as `devel_code` blocks.
- Several new tests: coverage now around 75% when CASTEP is available.
- New tutorials:

  - :ref:`MongoDB setup<mongo>`
  - :ref:`Spectral calculations with run3<run3_spectral>`
  - Example notebooks


New in release (0.8b) [03/08/2018]
----------------------------------

- Wholesale changes, complete refactoring of most of the code.
- Released open source under the MIT license!
- Documentation now hosted on `readthedocs <matador-db.readthedocs.org>`_,
- Workflows: chaining up job steps with run3:

  - spectral and phonons (combined DOS, dispersion calculations) with automated kpoint paths.
  - bulk modulus calculations and EOS fitting.

- New tutorials:

  - :ref:`Geometry optimisations with run3<run3_geom>`

- Temperature-dependent convex hulls (thanks to Angela Harper).
- New per-used configuration that allows changing of plotting styles, colourschemes, database names etc.
- Improvements to compute module:

  - automatically handle walltime constraints for Slurm and PBS.
  - estimate memory usage with CASTEP and skip if exceeds machine capacity,

- All CASTEP 18 keywords supported.
- Better support for electronic structure data, OptaDOS, NMR calculations, CIF files, partial occupancy.


New in version (0.7b) [13/04/2017]
----------------------------------

-  Ternary voltage curves.
-  Similarity/uniqueness filtering with element-projected PDFs.
-  Updated compute engine for remote calculations (see ``compute.py`` and new script ``oddjob``).
-  Improved test suite and full pip compatiblity.
-  Many bugfixes and usability changes.

New in version (0.6b) [01/06/2017]
----------------------------------

-  Intercalation voltage curves, e.g. ``matador voltage -c Li:SnS2``.
-  Ternary phase diagrams with heatmaps for structure prediction sampling, gravimetric capacity and formation enthalpy ``matador hull -c ABC --sampmap --efmap --capmap``.
-  Substructural similarity interface with Can Kocer's code, as proposed by `Yang et al., PRB (2014) <http://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.054102>`_.
