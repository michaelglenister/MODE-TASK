#!/bin/bash

file="input_pdb"
while IFS= read line
do
	prot=${line}
done <${file}

cd ..

g++ -I cpp/src/ ANM.cpp -o ANM

./ANM --pdb output/${prot}4_SCA.pdb --cutoff 24
