.. index:: cli

.. _cli:

Command-line usage
==================

Many of ``matador``'s entry points will make use of any ``.matadorrc`` config file they find. Customisation options can be found under :ref:`getting_started`.

Main entry points
-----------------

``matador`` has 3 main CLI entry points that are added to your ``$PATH`` on installation:

- ``matador``: used to interact with databases.
- ``run3``: used to run calculations at high-throughput on the local machine.
- ``dispersion``: used to plot electronic or vibrational properties from local files containing the output of the relevant calculations.

Each of these has many options/sub-commands, which can be explored using ``--help``.

Script entry points
-------------------

There are also several convenience scripts bundled with ``matador``. These scripts are not as well-supported as the main API, but please raise an issue on GitHub if they do not behave as expected. These scripts are also added to your ``$PATH`` on installation and can be found in the ``scripts/`` folder of a development install. If you have a workflow that only requires simple script usage but is not yet supported, or if you have a script you would like to contribute, please raise an issue on GitHub.

Analysis scripts
~~~~~~~~~~~~~~~~

- ``compare_pdfs``: compute PDF overlaps and plot differences from ``.res`` file input.
- ``pxrd_calculator``: compute, plot and export PXRD patterns from ``.cif`` file input.
- ``dbtools``: simple script for dropping MongoDB collections created by ``matador``.
- ``oddjob``: a script for spawning jobs on compute nodes with no queuing and shared filesystems, e.g. ``oddjob 'run3 <seed>' -n 1 3 5 -p cpu-`` will spawn ``run3`` processes on compute nodes with names ``cpu-1``, ``cpu-3`` and ``cpu-5``.
- ``fryan``: a stripped down version of ``matador query`` that operates on folders of ``.res`` files.
- ``plot_convergence``: plots the output of a ``run3`` convergence run.
- ``standardize_cell``: for an input ``.res`` file, standardize the cell (optionally to primitive) with spglib and write out a new file.

File converter scripts
~~~~~~~~~~~~~~~~~~~~~~

These scripts all use ``matador``'s internal file readers/writers to operate on lists of files provided at the CLI: e.g. ``x3y *.x`` will read all ``.x`` files in  a folder and output a series of ``.y`` files. This is not an all-to-all process, these scripts have only been written when needed. If you need quick access to a file converter, either write your own script based on the below, or request the converter script in an issue on GitHub.

- ``atoms3shx``: use ASE's file reader on any arbitrary file, and output a ``.res`` file.
- ``castep3shx``, ``castep3cell``: convert a CASTEP output file into ``.res`` or ``.cell`` respectively.
- ``cell3shx``, ``cif3shx``, ``magres3shx``: convert a CASTEP ``.cell``, ``.cif``, ``.magres`` or Quantum Espresso output file into a ``.res`` file.
- ``shx3cell``, ``shx3cif``, ``shx3pwscf``: convert a ``.res`` file into either ``.cell``, ``.cif``, or Quantum Espresso ``.in`` input file.
