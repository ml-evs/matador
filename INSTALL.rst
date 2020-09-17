.. _install:

Installation
============

If you have any issues with installation, feel free to raise an issue on GitHub outlining your approach and any errors you received.


Simple installation with pip
----------------------------

The matador package can be found on PyPI under the name `matador-db <https://pypi.org/project/matador-db>`_ and installed with
``pip install matador-db``, preferably in a fresh virtual environment (see conda instructions below). Extra dependencies may be installed with e.g. ``pip install matador-db[all]``.

Development installation with conda/pip
---------------------------------------

The tl;dr way to install matador, on e.g. a computing cluster, is as follows:

1. Clone the matador source onto your local machine ``git clone https://github.com/ml-evs/matador.git``.

Optional (but recommended) steps:

2. `Install conda <https://conda.io/miniconda.html>`_, if you have not already. There may be a package available already if you are using a supercomputer (e.g. `anaconda-compute/2.2.0-python3` on ARCHER 30/10/2017).
3. Create a new conda environment to install matador into (``conda create -n matador python=3.7``) and activate it with (``conda activate matador``).
4. Install some of the heavier requirements (e.g. NumPy and SciPy) through conda with ``conda install --yes --file requirements/requirements.txt``.

Required steps:

5. Run ``pip install .`` from inside the top-level matador directory, or ``pip install -e .`` for an editable developer install.
6. You now have a basic matador API installation, if you wish to use all matador features, install extra dependencies from the other requirements files inside ``requirements/`` using either conda or pip. If you wish to just install everything use ``pip install .[all]``.
7. To use matador, you will need to activate the conda environment from step 2, by running ``conda activate matador``. You will also need this in e.g. any job scripts. After installing the test dependencies with ``pip install .[test]``, you can test your installation using ``python -m unittest discover -v -b`` or simply ``py.test``. By default this will look for an MPI-enabled executable called ``castep`` on your ``$PATH`` to run CASTEP tests.

Troubleshooting
---------------

Below are some problems encountered on various machines that may be helpful:

1. (10/09/2019) When installing with ``conda``, if you receive the following error (or
   similar): ``/home/#####/.local/conda/envs/matador/compiler_compat/ld: build/temp.linux-x86_64-3.6/psutil/_psutil_common.o: unable to initialize decompress status for section .debug_info``, then you are using a modern compiler that breaks ``conda``'s attempts to be backwards compatible (in this case it was GCC 9). The simple fix is to rename/remove the copy of ``ld`` inside your conda environment (path in the message above) such that your system ``ld`` is used.
2. (10/10/2017) On some machines (e.g. ARCHER/Thomas) you may receive permissions errors at step 5; if so, try moving matador's `.git` and install again (``mv .git $HOME/matador_git_stash; pip install . ; mv $HOME/matador_git_stash .git``).
3. Some dependencies may not have compiled packages (wheels) for your distribution on PyPI, and may have compilation errors. In this case, you could try finding the relevant package on conda instead.
