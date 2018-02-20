# run3 for geometry optimisations with quantum espresso     
                                                           
## Requirements:                                             
   * folder full of res files                              
   * a template QE input file containing DFT parameters    
   * a desired kpoint spacing                              

First, we must make the QE input files
`shx3pwscf *.res --template vc-relax.template --kpoint_spacing 0.03`

Then, to run the optimisations with run3 (either interactively, or at the bottom of a job script) run3 must be called on all the *.in files.
`run3 -nc 4 --executable 'pw.x -i $seed.in > $seed.out' *.in`
in this case, our res file is called NaP.res, so the only calculation will be `mpirun -n 4 pw.x -i NaP.in > NaP.out`.
