.. _install:

Installation
============


Installing with conda/pip
-------------------------

The tl;dr way to install matador, on e.g. a computing cluster, is as follows:

1. Clone the matador source onto your local machine ``git clone https://bitbucket.org/ml-evs/matador.git``.

Optional (but recommended) steps:

2. Install conda, if you have not already (https://conda.io/miniconda.html) (choose the version suitable for you, the installers Python version does not matter that much, but it may as well be 3+!).
3. Create a new conda environment to install matador into (``conda create -n matador python=3.6``)
4. Install some of the heavier requirements (e.g. NumPy and SciPy) through conda with ``conda install --yes --file requirements/requirements.txt``. You may also want to install some of the optional dependencies in this manner, e.g. ``conda install --yes --file requirements/plotting_requirements.txt``.

Required steps:

4. Run ``pip install .`` from inside the top-level matador directory, or ``pip install -e .`` for an editable developer install.
5. You now have a basic matador API installation, if you wish to use all matador features, install extra dependencies from the other requirements files inside ``requirements/`` using either conda (preferably) or pip. e.g. for plotting, running your own database and Jupyter notebook visualisation functionality, use ``pip install .[plotting,db,viz]``. If you wish to just install everything use ``pip install .[all]``.
6. To use matador, you will need to work inside the conda environment from step 2, by running ``source activate matador``.


More esoteric installation on ARCHER (30th October 2017)
--------------------------------------------------------

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
