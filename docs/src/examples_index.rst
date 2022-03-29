.. index:: examples_index

.. _examples_index:

Tutorials and example notebooks
===============================

Real use-cases
--------------

- `Harper et al, Research Data Accompanying "Computational Discovery of Novel Copper Phosphide Conversion Anode for Lithium-Ion Batteries" <https://github.com/harpaf13/data.copper-phosphides>`_. |CuP Binder|

Example Jupyter notebooks
-------------------------

Interactive
~~~~~~~~~~~

These interactive examples can be run under Binder; some features that required system configuration (e.g. fonts when plotting) may not work correctly.

|Binder|

.. toctree::
    :maxdepth: 1

    notebooks/interactive/matador_plot_styles
    notebooks/interactive/voltage_from_res_files
    notebooks/interactive/pxrd_plotting
    notebooks/interactive/magres_plotting
    notebooks/interactive/pymatgen_and_ase_interface
    notebooks/interactive/visualisation

Non-interactive
~~~~~~~~~~~~~~~

These examples require external data, but can be used as example code.

.. toctree::
    :maxdepth: 1

    notebooks/non-interactive/hull
    notebooks/non-interactive/phonon_results
    notebooks/non-interactive/projected_spectral
    notebooks/non-interactive/projected_spectral_as_subplots

Tutorials
---------

These tutorials will guide you through common use cases of the command-line tools bundled with ``matador``.

.. toctree::
    :maxdepth: 2

    tutorials/run3
    tutorials/mongo

.. |Binder| image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/ml-evs/matador/master?filepath=examples/interactive/

.. |CuP Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/harpaf13/data.copper-phosphides/master?filepath=CuP_results.ipynb
