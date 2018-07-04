Installation
============


Installing with conda/pip
-------------------------

The tl;dr way to install matador, on e.g. a computing cluster, is as follows:

Optional conda steps:

1. install anaconda, if you have not already (https://conda.io/miniconda.html).
2. create a new conda environment (``conda create -n matador python=3.6 numpy scipy``)
3. manually install some of the heavier libraries (e.g. numpy and Scipy) through conda with ``conda install --yes --file requirements/requirements.txt``.

Required steps:

4. run ``pip install .`` from inside the top-level matador directory, or ``pip install -e .`` for an editable developer install.
5. you now have a basic matador installation, if you wish to use all matador features, install extra dependencies from the other requirements files inside ``requirements/``. e.g. for plotting, running your own database and Jupyter notebook visualisation functionality, use ``pip install .[plotting,db,viz]``.


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
