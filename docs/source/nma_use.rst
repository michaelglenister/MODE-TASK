NMA Scripts
====================================

Coarse grain
-------------------------------

Takes a protomer structure and coarse grains to select a set amount of CB atoms
**Command:** ::
	
	coarseGrain.py <options> --pdbFile <pdb file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *        	 | File       |``--pdb``           | PDB structure to coarse     |
|                        |            |                    | grain                       |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Course grain level     | Integer    |``--cg``            | Level by which to coarse    |
|                        |            |                    | grain the protein. Lower    |
|                        |            |                    | is less coarse grained.     |
|                        |            |                    | Default: 4                  |
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
| PDB file *             | File       |``--pdb``           | PDB input file              |
|                        |            |                    |                             |
+------------------------+------------+--------------------+-----------------------------+
| Cutoff                 | Integer    |``--cutoff``        | Cuttoff radius in Å.        |
|                        |            |                    | Default: 15                 |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| W matrix               | Text file of :math:`3N`     |
|                        | eigenvalues                 |
+------------------------+-----------------------------+
| VT matrix              | Text file of :math:`3N`\ x\ |
|                        | :math:`3N` eigenvectors.    |
|                        | Printed in rows             |
+------------------------+-----------------------------+
| U matrix               | Text file of :math:`3N`\ x\ |
|                        | :math:`3N` eigenvectors.    |
|                        | Printed in columns          |
+------------------------+-----------------------------+

Extract eigenvectors
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

Mean square fluctuation
-------------------------------

Calculates and returns the diagonals of the correlation matrix for a given set of modes.

Let's say that the user has performed NMA on two coarse grained models of the same protein, and now wants to compare
if the additional coarse graining decreased the accuracy. If we obtain the same mean square fluctuations for
each residue, then in each model we can say that the results are comparable regardless of the coarse graining
level. Obviously, we must compare only the residues that are common in each model. Hence we specify common residues

**Command:** ::

	meanSquareFluctuations.py <options> --pdb <PDB file> --wMatrix <text file> --vtMatrix <text file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| PDB file *             | File       |``--pdb``           | PDB input file              |
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

Identifies modes responsible for the conformational change of a molecule.

**Command:** ::

	conformationMode.py <options> --pdbConfAligned <PDB file> --pdbProtAligned <PDB file> --pdbANM <PDB file> --vtProtomer <text file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Aligned PDB file *     | File       |``--pdbConfAligned``| PDB file of the             |
|                        |            |                    | conformational change       |
+------------------------+------------+--------------------+-----------------------------+
| PDB *                  | File       |``--pdbANM``        | PDB file that was useed to  |
|                        |            |                    | run ANM                     |
+------------------------+------------+--------------------+-----------------------------+
| VT matrix file *       | File       |``--vtMatrix``      | Eigenavalues obtained from  |
|                        |            |                    | ANM script                  |
+------------------------+------------+--------------------+-----------------------------+
| Output file            | File       |``--output``        | Specify a name for the PDB	 |
|                        |            |                    | output file. Default:       |
|                        |            |                    | ModesOfConfChange.pdb       |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| Conformation file      | Text file with the overlap  |
|                        | and correlation of each     |
|                        | mode                        |
+------------------------+-----------------------------+


Mode visualisation
-------------------------------

Generates a trajectory with arrows that can be viewed in the tool VMD

**Command:** ::

	visualiseVector.py <options> --pdb <PDB file> --vectorFile <text file>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
| Coarse grained PDB     | File       |``--pdb``           | Coarse grained PDB input    |
| file *                 |            |                    | file                        |
+------------------------+------------+--------------------+-----------------------------+
| Mode index value       | Ingeter    |``--mode``          | Value specifying the index  |
|                        |            |                    | of the mode                 |
+------------------------+------------+--------------------+-----------------------------+
| Vector file *          | File       |``--vectorFile``    | File containing eigen       |
|                        |            |                    | vectors                     |
+------------------------+------------+--------------------+-----------------------------+

**Outputs:**

Outputs are generated in output/VISUALISE directory by default.

+------------------------+-----------------------------+
| Output                 | Description                 |
+========================+=============================+
| PDB file               | Output PDM to be opened in  |
|                        | VMD                         |
+------------------------+-----------------------------+
| Arrows file            | Tcl script that can be      |
|                        | copied into the VMD TK      |
|                        | console                     |
+------------------------+-----------------------------+
