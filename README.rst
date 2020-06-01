=======
matador
=======

|PyPI Version| |GH Actions| |Coverage Status| |Documentation Status| |MIT License|

matador is an aggregator, manipulator and runner of first-principles calculations, written with a bent towards battery electrode materials.
The source can be found on `GitHub <https://github.com/ml-evs/matador>`_ and online documentation is hosted on `ReadTheDocs <https://docs.matador.science>`_.

Example Jupyter notebooks and tutorials can be found `online <https://docs.matador.science/en/latest/examples_index.html>`_ or in the ``examples/`` folder of the matador source code.

Written & maintained by `Matthew Evans <https://ml-evs.science>`_ (2016-).


.. image:: docs/src/img/lipzn.png
   :name: LiPZn
   :align: center

Installation
------------

In the simplest case (e.g. you already have Python 3.6+ set up), ``pip install matador-db`` is sufficient to get up and running, preferably in a fresh virtual environment. More detailed instructions can be found in the `Installation instructions <https://docs.matador.science/en/latest/install.html>`_.

Upgrading to the latest version should be as simple as ``pip install -U matador-db``.

Usage
------

``matador`` is primarily a Python *library* that can be used inside Python scripts/modules to create a custom workflow. There are, however, several command-line scripts bundled with ``matador`` itself. All of these scripts are listed under `CLI Usage <https://docs.matador.science/en/latest/cli.html>`_.

For basic command-line usage, please explore the help system for command. Common workflows can be found inside ``examples/`` and in the `online docs <http://docs.matador.science/en/latest/examples_index.html>`_.

Please consult the full `Python API documentation <http://docs.matador.science/en/latest/modules.html>`_ for programmatic usage.

Core functionality
-------------------

The API has many features that can be explored in the examples and API documentation. As a summary, ``matador`` can be used for:

- Scraping of CASTEP (and Quantum Espresso) input/output files into flexible Python dictionaries/models.
- The creation and curation of MongoDB collections of geometry optimisation, calculations, with a powerful querying CLI/API.
- Customisable, publication-ready plots for all models, e.g. phase diagrams, PDF, PXRD, voltage profiles, electronic/vibrational bandstructures etc.
- Automated high-throughput geometry optimisations, electronic and vibrational properties using CASTEP (and Quantum Espresso) with ``run3``. Tested on several supercomputers.
- Creation of phase diagrams and electrochemical voltage profiles from the results of DFT calculations.


.. |PyPI Version| image:: https://img.shields.io/pypi/v/matador-db?label=PyPI&logo=pypi
.. |GH Actions| image:: https://img.shields.io/github/workflow/status/ml-evs/matador/Run%20tests/master?label=master&logo=github
   :target: https://github.com/ml-evs/matador/actions?query=branch%3Amaster
.. |MIT License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://github.com/ml-evs/matador/blob/master/LICENSE
.. |Coverage Status| image:: https://img.shields.io/codecov/c/gh/ml-evs/matador/master?logo=codecov
  :target: https://codecov.io/gh/ml-evs/matador
.. |Documentation Status| image:: https://readthedocs.org/projects/matador-db/badge/?version=stable
   :target: https://matador-db.readthedocs.io/en/stable/?badge=stable
