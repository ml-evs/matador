Changelog
=========

New in current version (0.8b) [03/08/2018]
------------------------------------------

* Wholesale changes, complete refactoring of most of the code.
* Released open source under the MIT license!
* Documentation now hosted on `readthedocs <matador-db.readthedocs.org>`_,
* Workflows: chaining up job steps with run3:
    * spectral and phonons (combined DOS, dispersion calculations) with automated kpoint paths.
    * bulk modulus calculations and EOS fitting.
* Temperature-dependent convex hulls (by Angela Harper).
* New per-used configuration that allows changing of plotting styles, colourschemes, database names etc.
* Improvements to compute module:
    * automatically handle walltime constraints,
    * estimate memory usage with CASTEP and skip if exceeds machine capacity,
* All CASTEP 18 keywords supported.
* Better support for electronic structure data, OptaDOS, NMR calculations, CIF files, partial occupancy.


New in version (0.7b) [13/04/2017]
----------------------------------

*  Ternary voltage curves.
*  Similarity/uniqueness filtering with element-projected PDFs.
*  Updated compute engine for remote calculations (see ``compute.py`` and new script ``oddjob``).
*  Improved test suite and full pip compatiblity.
*  Many bugfixes and usability changes.

New in version (0.6b) [01/06/2017]
----------------------------------

*  Intercalation voltage curves, e.g. ``matador voltage -c Li:SnS2``.
*  Ternary phase diagrams with heatmaps for structure prediction
   sampling, gravimetric capacity and formation enthalpy
   ``matador hull -c ABC --sampmap --efmap --capmap``.
*  Substructural similarity interface with Can Kocer's code, as proposed
   by `Yang et al., PRB
   (2014) <http://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.054102>`__.
