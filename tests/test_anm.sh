#!/bin/bash
cd ..

g++ -I cpp/src/ ANM.cpp -o ANM

./ANM --pdb output/3VBSPent4_SCA.pdb --cutoff 2
