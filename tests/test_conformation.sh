#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./conformationMode.py --pdbAligned tests/${prot}.pdb --pdbProtomerAligned output/${prot:0:4}_SCA.pdb --pdbSca output/${prot}4_SCA.pdb --vtProtomer output/VT_values.txt
