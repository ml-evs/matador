#!/bin/bash
cd docs
rm matador*.rst modules.rst
cd ../
sphinx-apidoc --separate -o docs matador matador/tests/* setup.py
cd docs
make html
