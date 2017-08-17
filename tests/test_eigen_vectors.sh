#!/bin/bash
cd ..

g++ -I cpp/src/ getEigenVectors.cpp -o getEigenVectors

./getEigenVectors --vt_values output/VT_values.txt --total_residues 627 --first_mode 617
