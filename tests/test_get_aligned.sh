#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./getAlignedCG.py --pdbAligned tests/${prot}.pdb --pdbCG output/ComplexCG.pdb
