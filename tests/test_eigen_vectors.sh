#!/bin/bash
cd ..

g++ -I cpp/src/ getEigenVectors.cpp -o getEigenVectors

./getEigenVectors --pdb output/VT_values.txt
