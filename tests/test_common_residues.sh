#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./commonResidues.py --fullCapsid tests/${prot}.pdb --protomer output/${prot}4_SCA.pdb
