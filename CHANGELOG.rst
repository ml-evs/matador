Changelog
=========

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
- Fixed issue with multinode MPI tasks by switching to ``proc.communicate()`` after an initial polling stage (#37)
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
------------------------------------------

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
