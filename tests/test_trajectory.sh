#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

./trajectoryPentamer.py --pdb output/${prot}4_SCA.pdb --modeFile output/ProtomerMode.txt
