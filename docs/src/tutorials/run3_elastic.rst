.. index:: run3_elastic

.. highlight:: bash

.. _run3_elastic:


Example 4: Bulk moduli with CASTEP and run3
-------------------------------------------

In this tutorial we will calculate the bulk moduli of diamond, silicon and lithium, using CASTEP and run3. 

The general process for calculating the bulk modulus is as follows:

1. Relax a structure to its equilibrium volume :math:`V_0` (ideally relaxing the positions
   too).
2. Run total energy calculations at a series of perturbed volumes :math:`\alpha V_0`
   for :math:`\alpha \in [0.95, 1.05]`.
3. Fit a particular form of equation of state to the resulting :math:`E(V)` curve and extract the bulk modulus.

Note: there is a supported set of scripts for calculating all the elastic constants in CASTEP, which
can be found on `GitHub <https://github.com/andreww/elastic-constants>`_ with an additional tutorial on the `CASTEP website <http://www.castep.org/Tutorials/ElasticConstants>`_.

As with the other tutorials, run3 expects to find a series of structures as ``.res`` files and single ``$seed.cell`` and ``$seed.param`` files. The files for this example can be found in ``examples/run3_elastic``. The ``.cell`` file in this tutorial is basically identical to that in the geometry optimisation tutorial, but the ``.param`` file this time contains the line ``task: bulk_modulus``. This is *not* a valid CASTEP task (as of 2019), but instead run3 will capture this task and use it to spawn a series of single point jobs.::

    $ cat bulk_mod.cell
    kpoints_mp_spacing: 0.05
    snap_to_symmetry
    symmetry_tol: 0.01
    symmetry_generate
    %block species_pot
    QC5
    %endblock species_pot

    $ cat bulk_mod.param
    task: bulk_modulus
    cut_off_energy: 300 eV
    write_bib: false
    write_checkpoint: none
    geom_max_iter: 100
    xc_functional: LDA


Once the files are in place, calling run3 as ``run3 bulk_mod`` will begin to relax
the first structure before deforming it. You can track the calculations
progress in the log files found in ``logs/``.

After the calculations have completed, your completed folder should look
something like this::

   $ ls completed
    C-OQMD_675640-CollCode28857.bands             Li-OQMD_30676-CollCode642106_bulk_mod.castep
    C-OQMD_675640-CollCode28857.castep            Li-OQMD_30676-CollCode642106_bulk_mod.cell
    C-OQMD_675640-CollCode28857.cell              Li-OQMD_30676-CollCode642106_bulk_mod.cst_esp
    C-OQMD_675640-CollCode28857.cst_esp           Li-OQMD_30676-CollCode642106_bulk_mod.param
    C-OQMD_675640-CollCode28857.geom              Li-OQMD_30676-CollCode642106_bulk_mod.png
    C-OQMD_675640-CollCode28857.param             Li-OQMD_30676-CollCode642106_bulk_mod.res
    C-OQMD_675640-CollCode28857.res               Li-OQMD_30676-CollCode642106_bulk_mod.results
    C-OQMD_675640-CollCode28857_bulk_mod.bands    Si-OQMD_5714-CollCode29287.bands
    C-OQMD_675640-CollCode28857_bulk_mod.castep   Si-OQMD_5714-CollCode29287.castep
    C-OQMD_675640-CollCode28857_bulk_mod.cell     Si-OQMD_5714-CollCode29287.cell
    C-OQMD_675640-CollCode28857_bulk_mod.cst_esp  Si-OQMD_5714-CollCode29287.cst_esp
    C-OQMD_675640-CollCode28857_bulk_mod.param    Si-OQMD_5714-CollCode29287.geom
    C-OQMD_675640-CollCode28857_bulk_mod.png      Si-OQMD_5714-CollCode29287.param
    C-OQMD_675640-CollCode28857_bulk_mod.res      Si-OQMD_5714-CollCode29287.res
    C-OQMD_675640-CollCode28857_bulk_mod.results  Si-OQMD_5714-CollCode29287_bulk_mod.bands
    Li-OQMD_30676-CollCode642106.bands            Si-OQMD_5714-CollCode29287_bulk_mod.castep
    Li-OQMD_30676-CollCode642106.castep           Si-OQMD_5714-CollCode29287_bulk_mod.cell
    Li-OQMD_30676-CollCode642106.cell             Si-OQMD_5714-CollCode29287_bulk_mod.cst_esp
    Li-OQMD_30676-CollCode642106.cst_esp          Si-OQMD_5714-CollCode29287_bulk_mod.param
    Li-OQMD_30676-CollCode642106.geom             Si-OQMD_5714-CollCode29287_bulk_mod.png
    Li-OQMD_30676-CollCode642106.param            Si-OQMD_5714-CollCode29287_bulk_mod.res
    Li-OQMD_30676-CollCode642106.res              Si-OQMD_5714-CollCode29287_bulk_mod.results
    Li-OQMD_30676-CollCode642106_bulk_mod.bands

You can see the results of the bulk modulus fits for each structure in
``completed/*_bulk_mod.results``, along with plots of the fits saved as pngs.

Compare these LDA values from the three different fits to the experimental values from Wikipedia:

+----------+--------------------------+------------------------------+
| Material | Expt. Bulk modulus (GPa) | Predicted Bulk Modulus (GPa) |
+==========+==========================+==============================+
| Lithium  | 11                       | 17.6, 17.7, 17.8             |
+----------+--------------------------+------------------------------+
| Silicon  | 98                       | 97.5, 98.0, 97.4             |
+----------+--------------------------+------------------------------+
| Diamond  | 442                      | 503, 496, 496                |
+----------+--------------------------+------------------------------+
