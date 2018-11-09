=======
matador
=======

|Bitbucket Pipelines| |Coverage Status| |Documentation Status| |MIT License|

matador is an aggregator, manipulator and runner of first-principles
calculations, written with a bent towards battery electrode materials. 
The source can be found on `Bitbucket <https://bitbucket.org/ml-evs/matador>`_.

Written by `Matthew Evans <https://www.ml-evs.github.io>`_ (2016-). 


.. image:: docs/src/img/lipzn.png
   :name: LiPZn
   :align: center

Core functionality
-------------------

The core functionality can be summarised by the various sub-commands of the
command-line interface. The API has many more features that can be explored
in the examples and API documentation.


1. The scraping of CASTEP/Quantum Espresso output files into flexible
   Python dictionaries with a sensible schema via ``matador import``.
2. The transferal of the scraped dictionaries into a MongoDB database.
3. Powerful CLI querying of the database, with a focus on energy storage
   applications using ``matador query``.
4. Calculation and presentation of binary and ternary phase diagrams
   with ``matador hull``.
5. 0K voltage curves for binary and ternary systems, as well as arbitrary intercalation electrodes, using
   ``matador voltage``.
6. Atomic species swapping and polishing from previous calculations using 
   ``matador swaps``.
7. Automated high-throughput geometry optimisations, electronic and vibrational properties, 
   plus convergence tests, using CASTEP or Quantum Espresso, with the ``run3`` entrypoint.
   Tested on several supercomputers.
8. Refinement of structural and chemical data, and symmetries powered by ``spglib``, via
   ``matador refine``.

Usage
------

For basic command-line usage, please explore the help system for each sub-command. Common workflows can be found under ``examples/`` and :ref:`tutorials_index`.

The full Python API documentation can be found under :ref:`modules`.

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

License
--------

matador is available under the `MIT License <https://bitbucket.org/ml-evs/matador/src/master/LICENSE>`_.

.. |Bitbucket Pipelines| image:: https://img.shields.io/bitbucket/pipelines/ml-evs/matador/master.svg
   :target: https://bitbucket.org/ml-evs/matador/addon/pipelines/home
.. |MIT License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://bitbucket.org/ml-evs/matador/src/master/LICENSE
.. |Coverage Status| image:: https://codecov.io/bb/ml-evs/matador/branch/master/graph/badge.svg
  :target: https://codecov.io/bb/ml-evs/matador
.. |Documentation Status| image:: https://readthedocs.org/projects/matador-db/badge/?version=latest
   :target: https://matador-db.readthedocs.io/en/latest/?badge=latest
