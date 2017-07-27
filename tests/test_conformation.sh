#!/bin/bash
cd ..

./conformationMode.py --pdbAligned tests/3VBSPent.pdb --pdbProtomerAligned output/3VBS_SCA.pdb --pdbSca output/3VBSPent4_SCA.pdb --vtProtomer output/VT_values.txt
