Installation
============


Installing with conda/pip (OUT OF DATE)
---------------------------------------

The tl;dr way to install matador, on e.g. a computing cluster, is as follows:

-  create a new conda environment (``conda create matador python=3.6 numpy scipy``)
-  manually install some of the heavier libraries (e.g. numpy and Scipy)
   through conda with ``conda install --yes --file requirements.txt``.
-  run ``pip install .`` from inside the top-level matador directory.
-  this will provide a basic matador installation, if you wish to use
   e.g. bandstructure/plotting/jupyter notebook visualisation
   functionality, use ``pip install .[plotting,bandstructure,viz]``.

Esoteric installation on ARCHER (30th October 2017)
---------------------------------------------------

These steps assume you have your ``work`` directory symlinked inside
``$HOME``.

1. ``module load anaconda-compute/2.2.0-python3``
2. ``conda create -p $HOME/work/.conda/matador-env python=3.6 numpy scipy``
3. ``source activate $HOME/work/.conda/matador-env``
4. ``git clone git@bitbucket.org:ml-evs/matador.git``
5. ``cd matador``
6. Due to some weird permissions on ARCHER, we must move the git folder
   so that ``pip`` works...
7. ``mv .git $HOME/matador_git_stash; pip install . ; mv $HOME/matador_git_stash .git``
