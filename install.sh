#!/bin/bash

virtualenv venv
source venv/bin/activate

pip install --upgrade pip
pip install numpy
pip install matplotlib 
pip install cython
pip install mdtraj
