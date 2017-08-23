#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./conformationMode.py --pdbConfAligned tests/Apart_PentamerAligned.pdb --pdbProtAligned output/ComplexCG.pdb --pdbANM output/ComplexCG.pdb --vtMatrix output/VT_values.txt
