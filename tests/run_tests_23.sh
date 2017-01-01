#!/bin/bash
echo "Testing python2..."
python2 --version
python2 -m unittest discover . "*_test.py"
echo "Testing python3..."
python3 --version
python3 -m unittest discover . "*_test.py"
