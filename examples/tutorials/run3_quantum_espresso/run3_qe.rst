.. index:: run3_qe

.. highlight:: bash

.. _run3_elastic:

Example 5: Geometry optimisations with Quantum Espresso and run3
================================================================

This tutorial uses the files found in ``examples/run3_quantum_espresso/vc-relax`` to relax
a folder of res files in a similar way to the CASTEP tutorial. Quantum Espresso (QE) uses a single
input file (as opposed to CASTEP's two), and these must be prepared beforehand (rather than
generated per structure with run3).

First things first, set up your directory so that you have:

* a load of res files
* a template QE input file containing DFT parameters

First, we make the QE input files using the ``shx3pwscf`` script:
``shx3pwscf *.res --template vcr.template --kpoint_spacing 0.03``

Then, to run the optimisations with run3 (either interactively, or at the bottom of a job script), run3 must be called on all the ``*.in`` files:
``run3 --redirect "$seed.out" -nc 4 --mode generic --executable 'pw.x -i $seed.in' *.in``
The ``$seed`` variables will be expanded by run3 to take the values of the ``.res`` file names. In this case, our only res file is called NaP.res, so the only calculation will be called as ``mpirun -n 4 pw.x -i NaP.in > NaP.out`` (or equivalent).
