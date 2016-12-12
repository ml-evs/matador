usage: matador [-h] [-v]
               
               {stats,query,import,rebuild,pdffit,hull,voltage,swaps,polish,refine}
               ...

MATerial and Atomic Database Of Refined structures.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  valid sub-commands

  {stats,query,import,rebuild,pdffit,hull,voltage,swaps,polish,refine}
    stats               print some stats about the database.
    query               query and extract structures from the database
    import              import new structures in folder into database
    rebuild             rebuild whole database.
    pdffit              provide experimental .gr file and fit to calculated
                        PDF of structures in query
    hull                create a convex hull from query results (currently
                        limited to binaries)
    voltage             plot a voltage curve from query results (currently
                        limited to binaries)
    swaps               perform atomic swaps on query results
    polish              re-relax a series of structures with new parameters.
    refine              refine database structures

Written by Matthew Evans (2016).

API documentation: www.tcm.phy.cam.ac.uk/~me388/matador_docs/