Normal Mode Analysis
====================================

Coarse grain
-------------------------------

Takes a protomer structure and coarse grains to select a set amount of C-Beta

**Command:** ::
	
	coarseGrain.py <options> --pdbFile <pdb file> --cg <int> --startingAtom <int> --protomerAtoms <int>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *        	 | File       |``--pdbFile``       | PDB input file              |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Course grain level     | Integer    |``--cg``            | Default: 4                  |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Starting atom          | Integer    |``--startingAtom``  | Residue number of the    	 |
|                        |            |                    | starting atom.              |
|                        |            |                    | Default: 1                  |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PDB file               | Coarse grained protein      |
|                        |                             |
+------------------------+-----------------------------+

ANM
-------------------------------

**Compile:** ::

    g++ -I cpp/src/ ANM.cpp -o ANM

**Command:** ::

	ANM <options> --pdb <pdb file> --cutoff <int>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *             | File       |``--pdb``           |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Cutoff                 | Integer    |``--cutoff``        | Default: 24                 |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| W matrix               |                             |
|                        |                             |
+------------------------+-----------------------------+
| VT matrix              |                             |
|                        |                             |
+------------------------+-----------------------------+
| U matrix               |                             |
|                        |                             |
+------------------------+-----------------------------+

Get eigen vectors
-------------------------------

**Compile:** ::

	g++ -I cpp/src/ getEigenVectors.cpp -o getEigenVectors

**Command:** ::

	getEigenVectors <options> 

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| VT matrix file *    	 | File       |``--vt_values``     | VT values from ANM script   |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
|                        |                             |
|                        |                             |
+------------------------+-----------------------------+

Common residues
-------------------------------

Takes two pdb models and determines the common residues

**Command:** ::

	commonResidues.py --fullCapsid <pdb file> --protomer <pdb file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Full capsid file *     | File       |``--fullCapsid``    |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Protomer file *        | File       |``--protomer``	   |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PDB file               |                             |
|                        |                             |
+------------------------+-----------------------------+

Mean square fluctuation
-------------------------------

Calculates and Returns Diagonals of Correlated Matrix for a given set of modes

Lets say that the user has performed NMA on two coarse grained models of the same protein and now wants to compare
to see if the additional coarse graining decreased the accuracy. If we obtain the same mean square fluctuations for
each residue then in each model then we can say that the results are comparable regardless of the coarse graining
level. But obviously must compare only the residues that are common in each model. hence we specify commonResidues

**Command:** ::

	meanSquareFluctuations.py --commonResidues <file> --pdbProtomer --totalModes --firstMode --lastMode --firstResidue --lastResidue --wMatrix --vtMatrix

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Common residue file *  | File       |``--commonResidues``|                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Protomer PDB file *    | File       |``--pdbProtomer``   |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| First mode             | Integer    |``--firstMode``	   | When unassigned the script  |
|                        |            |                    | will use the last six modes |
+------------------------+------------+--------------------+-----------------------------+
| Last mode              | Integer    |``--lastMode``	   | When unassigned the script  |
|                        |            |                    | will use the last six modes |
+------------------------+------------+--------------------+-----------------------------+
| W matrix file *        | File       |``--wMatrix``	   | W values from ANM script    |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``	   | VT values from ANM script   |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+ 

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
|                        |                             |
|                        |                             |
+------------------------+-----------------------------+


Conformation mode
-------------------------------

Identifies Modes responsible for conformational change for a molecule wth 15 copies of each atom

**Command:** ::

	conformationMode.py 

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Aligned PDB file *     | File       |``--pdbAligned``    |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Aligned protomer PDB   | File       |``--pdbProtAligned``|                             |
| file *                 |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| PDB *                  | File       |``--pdbSca``        |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtProtomer``    | VT values from ANM script   |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
|                        |                             |
|                        |                             |
+------------------------+-----------------------------+

Get aligned
-------------------------------

Creates a PDB for a multiple protomer structure, containing co-ords of an aligned PDB structure

**Command:** ::

	getAligned.py 

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Aligned PDB file *     | File       |``--pdbAligned``    |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| PDB *                  | File       |``--pdbSca``        |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
|                        |                             |
|                        |                             |
+------------------------+-----------------------------+

Trajectory pentamer
-------------------------------

Makes a trajectory of 100 PDB files

**Command:** ::

	trajectoryPentamer.py 

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Coarse grained PDB     | File       |``--pdb``           |                             |
| file *                 |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
|                        | File       |``--modeF``         |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
|                        | File       |``--modeL``         |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Mode file *            | File       |``--modeFile``      | File containing eigen       |
|                        |            |                    | vectors                     |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
|                        |                             |
|                        |                             |
+------------------------+-----------------------------+


