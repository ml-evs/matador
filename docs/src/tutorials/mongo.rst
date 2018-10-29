.. index:: mongo

.. highlight:: bash

.. _mongo:

Setting up your own database
============================

``matador`` uses MongoDB to store and operate on results from first principles calculations. In this tutorial, we will set up a MongoDB database containing some dummy test data.

Installing and configuring MongoDB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instructions on how to install MongoDB community edition can be found at the `MongoDB website <https://docs.mongodb.com/manual/administration/install-community/>`_. If you are using Linux, you may find that your OS provides a MongoDB package that you can install.

Once installed, you can run your MongoDB server with the ``mongod`` command. An example configuration file can be found in ``matador/config/example_mongod.conf``, you may wish to change the port and default data directory (where the raw database is written to disk). This command will run persistently so you may wish to run it in a persistent terminal (using e.g. GNU screen or tmux), or as an OS-wide service (e.g. using systemctl or your OS's equivalent).

.. warning::

   It is up to the user to ensure that their MongoDB server is secure! This guide will assume access from ``localhost`` only, and thus will not set a password. If you do not wish external parties to have access to your database, make sure your server is running on secure/closed port, or set up authorisation as outlined on the `MongoDB website <https://docs.mongodb.com/manual/tutorial/enable-authentication/>`_.

.. tip::

   You should ensure the path found under the ``storage.dbPath`` option in your configuration file exists before running ``mongod``, and that you have access to it!

With a MongoDB server running, now we can import some structures. Navigate to ``examples/mongo/castep_files`` and run ``matador import``. If this is your first time using ``matador``, you should be promoted to make a configuration file. Follow the instructions and enter the values you used to set up your database. You will be given the option to test the connection to your database at the end.

If all is successful, you should see that 3 structures are added to your database, along with 4 failures. You can have a look at the spatula output files to see what happened. To finish this tutorial, use ``matador query`` to look at the structures you have just imported, before moving on to the full query tutorial at :ref:`query_tutorial`.
