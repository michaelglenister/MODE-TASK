#!/bin/bash

cd ..
export DISPLAY=0.0
./internal_pca.py --trj tests/pca_test_trj.xtc --top tests/complex.pdb
