#!/bin/bash
cd ..

g++ -I cpp/src/ ANM.cpp -o ANM

./ANM --pdbFile tests/3VBSPent4_SCA.pdb
