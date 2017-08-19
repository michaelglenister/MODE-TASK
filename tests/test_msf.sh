#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./meanSquareFluctuations.py --commonResidues output/common_residues --pdbProtomer output/${prot}4_SCA.pdb --wMatrix output/W_values.txt --vtMatrix output/VT_values.txt
