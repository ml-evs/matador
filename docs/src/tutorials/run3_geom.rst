.. index:: run3_geom

.. _run3_geom:

Example 1: High-throughput geometry optimisations with CASTEP
-------------------------------------------------------------

.. _ex1:


Example 1.1: Using run3 locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


In this example, we will suppose that you want to perform a geometry optimisation on several different polymorphs of TiO\ :sub:`2` from the ICSD. The files for this example can be found in ``examples/run3_tutorial``, `here <https://github.com/ml-evs/matador/blob/develop/examples/run3_tutorial>`_.

Setting up the files
^^^^^^^^^^^^^^^^^^^^

By default, run3 expects the following files in the job folder:

* one ``.res``/SHELX file `per structure` to optimise
* one ``$seed.cell`` file that contains the CASTEP CELL keywords common to every structure, e.g. pseudopotentials, kpoint spacing, etc.
* one ``$seed.param`` file that contains the CASTEP PARAM keywords common to every calculation, e.g. ``cut_off_energy``, ``xc_functional``, ``task``, ``geom_max_iter``, etc.
* any required pseudopotential files.

.. tip:: If you have a database set up with structures from the OQMD these could be obtained via ``matador query --db oqmd_1.1 -f TiO2 --icsd --res``.

.. tip:: Alternatively, you can turn many file types into ``.res`` using the various ``<format>3shx`` scripts (shx standing for SHELX), e.g. ``cell3shx *.cell``.

The job folder should look something like this::

    $ ls
    O2Ti-OQMD_112497-CollCode171670.res  O2Ti-OQMD_2575-CollCode9852.res     O2Ti-OQMD_7500-CollCode41493.res
    O2Ti-OQMD_117323-CollCode182578.res  O2Ti-OQMD_3070-CollCode15328.res    O2Ti-OQMD_84685-CollCode97008.res
    O2Ti-OQMD_13527-CollCode75179.res    O2Ti-OQMD_31247-CollCode657748.res  O2Ti-OQMD_97161-CollCode154036.res
    O2Ti-OQMD_19782-CollCode154035.res   O2Ti-OQMD_5979-CollCode31122.res    TiO2.cell
    O2Ti-OQMD_2475-CollCode9161.res      O2Ti-OQMD_7408-CollCode41056.res    TiO2.param

with ``.param`` file containing::

    $ cat TiO2.param
    task                 : geometryoptimization
    xc_functional        : LDA 
    cut_off_energy       : 300.0 eV
    geom_force_tol       : 0.1
    spin_polarized       : false
    fix_occupancy        : false
    max_scf_cycles       : 100
    opt_strategy         : speed
    page_wvfns           : 0
    perc_extra_bands     : 40
    num_dump_cycles      : 0
    backup_interval      : 0
    geom_method          : LBFGS
    geom_max_iter        : 300
    mix_history_length   : 20
    finite_basis_corr    : 0
    fixed_npw            : false
    write_cell_structure : true
    write_checkpoint     : none
    write_bib            : false
    bs_write_eigenvalues : false
    calculate_stress     : true

and ``.cell`` file containing::

    $ cat TiO2.cell
    kpoint_mp_spacing: 0.07
    
    %block species_pot
    QC5
    %endblock species_pot

    symmetry_generate
    symmetry_tol: 0.01
    snap_to_symmetry

.. highlight:: bash


Calling run3
^^^^^^^^^^^^

Once these files are in place, we can begin the geometry optimisations. To run the current host machine, simply call::

    $ run3 TiO2

This will start a single node CASTEP job on the current machine, using all available cores. If you are on a local cluster without a queuing system, and wish to run on several nodes at once (say ``node3``, ``node6`` and ``node8``), the oddjob script can be used as follows::

    $ oddjob 'run3 TiO2' -n 3 6 8

This will start 3 single node CASTEP jobs on the desired nodes. If instead your nodes are called ``cpu00010912``, ``cpu323232`` and ``cpu123123``, the ``--prefix`` flag is needed::

    $ oddjob 'run3 TiO2' --prefix cpu -n 00010912 323232 123123


.. tip:: On a supercomputer with a queuing system, e.g. PBS or slurm, run3 must be called in your submission script. Array jobs are typically an effective way of spreading out over multiple nodes. An example of this kind can be found in `example 1.2 <ex.1.2_>`__.


Monitoring your calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you look at the job folder as run3, er... runs, you will see several files and folders being created. Firstly, 3 ``.txt`` files will be made:

* ``jobs.txt``: this file contains a list of jobs that, at some point, __started__ running.
* ``finished_cleanly.txt``: this file lists jobs that completed without error.
* ``failures.txt``: this file lists jobs that produced an error.

Every structure in progress will have a ``<structure_name>.lock`` file to prevent clashes with other nodes.

Several folders will also be created:

* ``logs/``: log file per structure containing a complete history of the run.
* ``input/``: a backup of the starting configuration as a ``.res`` file.
* ``completed/``: all successful calculations will end up here, usually as a ``.res`` file with the final configuration, a concatenated ``.castep`` file containing every job step, and if requested (via ``write_cell_structure: true``), CASTEP's ``-out.cell`` file.
* ``bad_castep/``: all failed calculations end up here, including all auxiliary files.
* ``<node_name>/``: a folder is created per hostname (e.g. when running on multiple nodes) that contains the interim calculations. On failures/timeouts, all files in here are moved back to the main job folder.

Eventually, all jobs will hopefully be moved to ``completed/``, then you are done!


Example 1.1.1: High-throughput geometry optimisations with CASTEP with per-structure parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are a few occasions where you might need a custom ``.param`` file for each structure, for example, if using the implicit nanotube ``%devel_code`` in CASTEP.

These calculations are performed in exactly the same was as above, except a ``<structure_name>.param`` file must be made containing the required DFT parameters AND the nanotube parameters. In this case, run3 must now be called as::

    $ run3 --custom_params TiO2

.. tip:: If you have a .res file that contains a PyAIRSS "REM NT_PROPS" line, this will be ignored.


Example 1.2: High-throughput geometry optimisations with CASTEP on a supercomputer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each HPC facility has its own quirks, so in this example we will try to be as explicit as possible. The set up of the job is exactly the same as in `example 1 <ex1_>`__, but we now must add run3 to our job submission script. The following examples are for the SLURM queuing system on the BlueBear machine at the University of Birmingham and PBS on ARCHER (Tier-1), but run3 has also been tested on CSD3 (Tier-2), HPC Midlands-Plus (Tier-2), Thomas (Tier-2) and several local group-scale clusters.

Example 1.2.1: SLURM on BlueBear
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this job, we will submit a run3 job that performs CASTEP calculations across 2 24-core nodes per structure. Let us presume we have many thousand structures to run. The submission script looks as follows::

    $ cat run3.sub
    #!/bin/bash -l
    
    ###### MACHINE/USER-SPECIFIC OPTIONS ######
    
    #SBATCH --ntasks 48
    #SBATCH --nodes 2-2
    #SBATCH --time 24:00:00
    #SBATCH --qos <REDACTED> 
    ##SBATCH --qos bbshort
    #SBATCH --mail-type ALL
    #SBATCH --account=<REDACTED>

    module purge
    export PATH="$HOME/bin/CASTEP-17.21:$HOME/.conda/bin"
    module load bluebear
    module load mpi/impi/2017.1.132-iccifort-2017.1.132
    unset I_MPI_PMI_LIBRARY
    
    # RUN3 COMMANDS
    # (assuming installation guide followed at
    #  https://matador-db.readthedocs.io/en/latest/install.html)

    source activate matador
    run3 -nc 48 -v 4 --executable castep.mpi --ignore_jobs_file TiO2

Let's unpick a few of the flags used to call run3 here:

* ``-nc/--ncores``: the number of cores to use per structure, per calculation. It is often worth specifying this if more than one node is being used, as the correctness of run3's core counter is queue/machine-specific.
* ``-v 4``: sets the verbosity in the log file to the highest level.
* ``--ignore_jobs_file``: by default run3 will for both ``<structure>.lock`` files and entries in ``jobs.txt`` before running a new structure. It is often worth disabling the ``jobs.txt`` check if it is not expected that all structures complete in one job submission (see below).
  
try to call an executable called simply ``castep``. On many machines, CASTEP is installed as ``castep.mpi``.

Now to submit this script as a 200-node array job (i.e. running a maximum of 100 structures concurrently, depending on the queue), we call the following::

    $ sbatch --array=1-100 run3.job

It may be that this job is not large enough to optimise all structures within the walltime limit. In this case, it can be resubmitted using the same command. Jobs that were running when the walltime ran out should automatically be pushed back into the job folder so that they will be available to the next run3 call. In the event that this does not happen (for example MPI kept control of the Python thread for too long so the queuieng system interrupted run3's clean up), ``<hostname>`` folder
will be left hanging around in the main jobs folder. Jobs must then be manually made restartable by deleting ``<structure>.lock`` (and removing ``<structure>>`` from ``jobs.txt`` if not using ``--ignore_jobs_file``). It may also be that the intermediate CASTEP calculation was not copied over from the ``<hostname>`` folder: in this case, the CASTEP files can be updated by running::
    
    $ cp -u node*/*.castep .

from inside the root job folder.

Example 1.2.2: PBS on ARCHER
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ex.1.2:

Instructions are almost identical to the above, but the array job script looks a little different, for the same 100 copies of 2 node jobs (this time 24 cores per node)::

    $ cat run3.job
    #!/bin/bash --login
    # PBS job options (name, compute nodes, job time) # PBS -N is the job name (e.g. Example_MixedMode_Job)
    #PBS -N my_run3_job
    # PBS -l select is the number of nodes requested (e.g. 128 node=3072 cores)
    #PBS -l select=2
    # PBS -l walltime, maximum walltime allowed (e.g. 6 hours)
    #PBS -l walltime=24:00:00
    # Replace [budget code] below with your project code (e.g. t01)
    #PBS -A <REDACTED>
    #PBS -m abe
    #PBS -M <REDACTED>
    #PBS -J 1-100
    #PBS -r y

    # Make sure any symbolic links are resolved to absolute path
    export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

    # Change to the direcotry that the job was submitted from
    # (remember this should be on the /work filesystem)
    cd $PBS_O_WORKDIR

    source $HOME/.bashrc
    module load anaconda-compute/python3
    source activate $HOME/work/.conda/matador

    run3 --archer -v 4 -nc 48 KSnP

Notice here we have specified ``--archer``: again, run3 should be able to detect that ``mpirun`` is missing and thus try ``aprun``, but it can be worth specifying just in case. With PBS, the whole array can be submitted with just::

    $ qsub run3.job
