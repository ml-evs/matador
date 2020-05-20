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

In the simplest case, ``pip install matador-db`` is sufficient to get up and running, preferably in a fresh virtual environment. More detailed instructions can be found in the `Installation instructions <https://docs.matador.science/en/latest/install.html>`_.

Usage
------

For basic command-line usage, please explore the help system for each sub-command. Common workflows can be found inside ``examples/`` and in the `online docs <http://docs.matador.science/en/latest/examples_index.html>`_.

Please consult the full `Python API documentation <http://docs.matador.science/en/latest/modules.html>`_ for programmatic usage.

.. code-block:: text

    usage: matador [-h] [--version]
                   {stats,query,import,hull,voltage,changes,hulldiff,swaps,refine}
                   ...

    MATerial and Atomic Database Of Refined structures.

    optional arguments:
      -h, --help            show this help message and exit
      --version             show program's version number and exit

    subcommands:
      valid sub-commands

      {stats,query,import,hull,voltage,changes,hulldiff,swaps,refine}
        stats               print some stats about the database.
        query               query and extract structures from the database
        import              import new structures in folder into database
        hull                create a convex hull from query results (currently
                            limited to binaries and ternaries)
        voltage             plot a voltage curve from query results (currently
                            limited to binaries and ternaries)
        changes             view database changelog or undo additions to database
                            (NB: not deletions!)
        hulldiff            diff two convex hulls with the --compare flag.
        swaps               perform atomic swaps on query results
        refine              update structures in the database according to
                            specified --task

Core functionality
-------------------

The core functionality can be summarised by the various sub-commands of the
command-line interface above. The API has many more features that can be explored
in the examples and API documentation.

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
.. |Documentation Status| image:: https://readthedocs.org/projects/matador-db/badge/?version=latest
   :target: https://matador-db.readthedocs.io/en/latest/?badge=latest
