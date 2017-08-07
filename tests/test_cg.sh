#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./coarseGrain.py --pdbFile tests/${prot}.pdb --cg 4
