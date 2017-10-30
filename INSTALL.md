# Installation

Matador is pinned to various versions of its dependencies; the tl;dr way to install is to create a new conda environment, manually install some of the heavier libraries (e.g. numpy and Scipy) directly with conda, then run `pip install .` from inside the top-level matador directory.

## Installation on ARCHER (30th October 2017)

These steps assume you have your `work` directory symlinked inside `$HOME`.

1. `module load anaconda-compute/2.2.0-python3`
2. `conda create -p $HOME/work/.conda/matador python=3.6 numpy scipy`
3. `source activate $HOME/work/matador`
4. `git clone git@bitbucket.org:me388/matador.git`
5. `cd matador`
6. `pip install .`
