#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./getAligned.py --pdbAligned tests/${prot}.pdb --pdbSca output/${prot}4_SCA.pdb
