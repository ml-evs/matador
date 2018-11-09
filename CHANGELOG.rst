Changelog
=========

For planned developments, please see the `Boards tab <https://bitbucket.org/ml-evs/matador/addon/trello/trello-board>`_ on Bitbucket.

New in development version (rolling) (0.9a)
-------------------------------------------

- SpectralWorkflow extended to include latest projected dispersion curve developments from OptaDOS.
- Updated dispersion script to automatically perform naive Gaussian smearing if OptaDOS output not detected.
- New tutorials:

   + :ref:`MongoDB setup<mongo>`
   + :ref:`Spectral calculations with run3<run3_spectral>`

- Added ``--scratch_prefix`` to run3 to allow for temporary files to e.g. be written to faster filesystem with appropriate symlinks to work folder.

New in current release (0.8b) [03/08/2018]
------------------------------------------

- Wholesale changes, complete refactoring of most of the code.
- Released open source under the MIT license!
- Documentation now hosted on `readthedocs <matador-db.readthedocs.org>`_,
- Workflows: chaining up job steps with run3:

   + spectral and phonons (combined DOS, dispersion calculations) with automated kpoint paths.
   + bulk modulus calculations and EOS fitting.

- New tutorials:

   + :ref:`Geometry optimisations with run3<run3_geom>`

- Temperature-dependent convex hulls (by Angela Harper).
- New per-used configuration that allows changing of plotting styles, colourschemes, database names etc.
- Improvements to compute module:

   + automatically handle walltime constraints,
   + estimate memory usage with CASTEP and skip if exceeds machine capacity,

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
