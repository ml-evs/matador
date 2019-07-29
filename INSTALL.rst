.. _install:

Installation
============


Installing with conda/pip
-------------------------

The tl;dr way to install matador, on e.g. a computing cluster, is as follows:

1. Clone the matador source onto your local machine ``git clone https://bitbucket.org/ml-evs/matador.git``.

Optional (but recommended) steps:

2. `Install conda <https://conda.io/miniconda.html>`_, if you have not already. There may be a package available already if you are using a supercomputer (e.g. `anaconda-compute/2.2.0-python3` on ARCHER 30/10/2017).
3. Create a new conda environment to install matador into (``conda create -n matador python=3.6``) and activate it with (``conda activate matador``).
4. Install some of the heavier requirements (e.g. NumPy and SciPy) through conda with ``conda install --yes --file requirements/requirements.txt``. You may also want to install some of the optional dependencies in this manner, e.g. ``conda install --yes --file requirements/plotting_requirements.txt``.

Required steps:

5. Run ``pip install .`` from inside the top-level matador directory, or ``pip install -e .`` for an editable developer install. 
Note: On some machines (e.g. ARCHER/Thomas) you may receive permissions errors at this stage; if so, try moving matador's `.git` and install again (``mv .git $HOME/matador_git_stash; pip install . ; mv $HOME/matador_git_stash .git``)
6. You now have a basic matador API installation, if you wish to use all matador features, install extra dependencies from the other requirements files inside ``requirements/`` using either conda (preferably) or pip. e.g. for plotting, running your own database and Jupyter notebook visualisation functionality, use ``pip install .[plotting,db,viz]``. If you wish to just install everything use ``pip install .[all]``.
7. To use matador, you will need to activate the conda environment from step 2, by running ``conda activate matador``. You will also need this in e.g. any job scripts. You can test your installation using ``python -m unittest discover`` or simply ``py.test`` (if you have it installed). By default this will look for an executable called `castep` to run CASTEP tests, which are probably the most useful.
