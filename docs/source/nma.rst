NMA Scripts
====================================

Coarse grain
-------------------------------

Takes a protomer structure and coarse grains to select a set amount of C-Beta

**Command:** ::
	
	coarseGrain.py <options> --pdbFile <pdb file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *        	 | File       |``--pdb``           | PDB input file              |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Course grain level     | Integer    |``--cg``            | Default: 4                  |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Starting atom          | Integer    |``--startingAtom``  | Residue number of the    	 |
|                        |            |                    | starting atom.              |
|                        |            |                    | Default: 1                  |
+------------------------+------------+--------------------+-----------------------------+
| Output file            | File       |``--output``        | Specify a name for the PDB	 |
|                        |            |                    | output file.                |
|                        |            |                    | Default: ComplexCG.pdb      |
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

	ANM <options> --pdb <pdb file>

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

	g++ -I cpp/input/ getEigenVectors.cpp -o getEigenVectors

**Command:** ::

	getEigenVectors <options> --vt <text file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| VT matrix file *    	 | File       |``--vt``            | VT values from ANM script   |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Mode index *           | Integer    |``--mode``          | Specify the index of the    |
|                        |            |                    | mode you wish to target     |
+------------------------+------------+--------------------+-----------------------------+
| Direction              | Boolean    |``--direction``     | Direction of overlap        |
|                        | integer    |                    | correction                  |
|                        | (1 or -1)  |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| Eigen vectors file     | Text file containing a      |
|                        | list of eigen vectors       |
+------------------------+-----------------------------+

Common residues
-------------------------------

Extracts the common residues between two pdb models such as two protein conformations

**Command:** ::

	commonResidues.py <options> --conf1 <pdb file> --conf2 <pdb file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Conformation PDB *     | File       |``--conf1``         | Path to first PDB file      |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Conformation PDB *     | File       |``--conf2``         | Path to second PDB file     |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| Common residues file   | Text file containing common |
|                        | residues, this file is used |
|                        | as input to the mean square |
|                        | fluctuations script         |
+------------------------+-----------------------------+

Mean square fluctuation
-------------------------------

Calculates and Returns Diagonals of Correlated Matrix for a given set of modes

Lets say that the user has performed NMA on two coarse grained models of the same protein and now wants to compare
to see if the additional coarse graining decreased the accuracy. If we obtain the same mean square fluctuations for
each residue then in each model then we can say that the results are comparable regardless of the coarse graining
level. But obviously must compare only the residues that are common in each model. hence we specify commonResidues

**Command:** ::

	meanSquareFluctuations.py <options> --pdb <PDB file> --wMatrix <text file> --vtMatrix <text file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *             | File       |``--pdb``           |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Conformation PDB       | File       |``--pdfConf2``      | When assigned, caclulates   |
|                        |            |                    | mean square fluctautions    |
|                        |            |                    | based on in common residues |
|                        |            |                    | between the two proteins    |
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
| Beta values file       | Text file listing beta      |
|                        | values for common residues  |
+------------------------+-----------------------------+


Conformation mode
-------------------------------

Identifies Modes responsible for conformational change for a molecule wth 15 copies of each atom

**Command:** ::

	conformationMode.py <options> --pdbConfAligned <PDB file> --pdbProtAligned <PDB file> --pdbANM <PDB file> --vtProtomer <text file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Aligned PDB file *     | File       |``--pdbConfAligned``|                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Aligned protomer PDB   | File       |``--pdbProtAligned``|                             |
| file *                 |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| PDB *                  | File       |``--pdbANM``        |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``      | VT values from ANM script   |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Output file            | File       |``--output``        | Specify a name for the PDB	 |
|                        |            |                    | output file. Default:       |
|                        |            |                    | ModesOfConfChange.pdb       |
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

Generates a trajectory with arrows that can be viewed in the tool VMD

**Command:** ::

	visualiseVector.py <options> --pdb <PDB file> --vectorFile <text file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Coarse grained PDB     | File       |``--pdb``           |                             |
| file *                 |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
|                        | File       |``--mode``          |                             |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Vector file *          | File       |``--vectorFile``    | File containing eigen       |
|                        |            |                    | vectors                     |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

Outputs are generated in output/VISUALISE directory by default.

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PDB file               |                             |
|                        |                             |
+------------------------+-----------------------------+
| Arrows file            | Text file to draw arrows in |
|                        | the VMD visualizer          |
+------------------------+-----------------------------+
