#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./commonResidues.py --conf1 tests/${prot}.pdb --conf2 output/ComplexCG.pdb
