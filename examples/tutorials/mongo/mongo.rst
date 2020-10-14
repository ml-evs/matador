.. index:: mongo

.. highlight:: bash

.. _mongo:

Setting up your own database
============================

``matador`` uses MongoDB to store and operate on results from first principles calculations. In this tutorial, we will set up a MongoDB database containing some dummy test data.

Installing and configuring MongoDB
----------------------------------

Instructions on how to install MongoDB community edition can be found at the `MongoDB website <https://docs.mongodb.com/manual/administration/install-community/>`_. If you are using Linux, you may find that your OS provides a MongoDB package that you can install.

Once installed, you can run your MongoDB server with the ``mongod`` command. An example configuration file can be found in ``matador/config/example_mongod.conf``, you may wish to change the port and default data directory (where the raw database is written to disk). This command will run persistently so you may wish to run it in a persistent terminal (using e.g. GNU screen or tmux), or as an OS-wide service (e.g. using systemctl or your OS's equivalent).

.. warning::

   It is up to the user to ensure that their MongoDB server is secure! This guide will assume access from ``localhost`` only, and thus will not set a password. If you do not wish external parties to have access to your database, make sure your server is running on secure/closed port, or set up authorisation as outlined on the `MongoDB website <https://docs.mongodb.com/manual/tutorial/enable-authentication/>`_.

.. tip::

   You should ensure the path found under the ``storage.dbPath`` option in your configuration file exists before running ``mongod``, and that you have access to it!

Importing structures
--------------------

With a MongoDB server running, you can now import some structures. Navigate to ``examples/mongo/castep_files`` and run ``matador import``. If this is your first time using ``matador``, you should be prompted to make a configuration file. Follow the instructions and enter the values you used to set up your database. You will be given the option to test the connection to your database at the end.

If all is successful, you should see that 3 structures are added to your database, along with 4 failures. You can have a look at the spatula output files to see what happened. To finish this tutorial, use ``matador query`` to look at the structures you have just imported, before moving onto the command-line tutorial at :ref:`cli`.

Caring for your database
------------------------

Collections
^^^^^^^^^^^

MongoDB databases are made out of "collections" of data. To create a new collection, simply specify ``--db <collection_name>`` when performing an import. You may wish to create multiple collections of data, or one monolithic database.
   
The default collection can be set in your ``.matadorrc`` with the ``mongo.default_collection`` option. If you want to ensure that this database can only be imported to from a particular path, set the ``mongo.default_collection_file_path`` option.

A summary of the contents of existing collections can be generated using ``matador stats --db <collection_name>``. Similarly, a collection can be deleted using ``matador stats --db <collection_name> --delete``. If you wish to get a list of all existing collections, use ``matador stats --list``.

.. tip:: 
   In our group we combine these two approaches: a watched folder of "finalised" data is automatically scraped every day into a single monolithic collection (a cron job calls ``matador import``), interim calculations are then stored in collections per research project.

.. tip:: 
   Each collection also saves its own "changelog". You can view this changelog and undo changes with the ``matador changes`` interface.

Backing up
^^^^^^^^^^

You should always back up the data used to create a database, but if you wish, you can also directly backup the database itself using the MongoDB tool ``mongodump`` (`documentation <https://docs.mongodb.com/manual/reference/program/mongodump/>`_). Similarly, restoration can be performed using ``mongorestore`` (`documentation <https://docs.mongodb.com/manual/reference/program/mongorestore/>`_). You may also wish to read the general back up and restore tutorial on the `MongoDB website <https://docs.mongodb.com/manual/tutorial/backup-and-restore-tools/>`_.
