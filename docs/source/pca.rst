Principle Component Analysis
====================================

PCA on cartesian cordinates
-----------------------------

This programe performs the PCA (Principal Component Analysis) on a MD trajectory.

**1. SVD PCA** ::

To perform the singular value decomposition (SVD) PCA on C-alpha atoms of a MD trajectory (pca_test_trj.xtc)

	``pca.py -t pca_test_trj.xtc -p complex.pdb -at CA -pt svd``	

**2. Kernel PCA** ::

To perform the Kernel PCA with linear kernel on cartesian cordinates of C-alpha atoms of a MD trajectory (pca_test_trj.xtc)

	``pca.py -t pca_test_trj.xtc -p complex.pdb -at CA -pt kpca -kt linear``

**3. IncrementalPCA** ::

IPCA (IncrementalPCA) is a varient of normal PCA, which uses low-rank approximation of the input MD trajectory. It uses the amount of memory to store the input trajectory which is independent of trajectory size. IPCA is very usefull in case the size of trajectory size bigger than availaible computer memory.

	  ``pca.py -t pca_test_trj.xtc -p complex.pdb -at CA -pt ipca``

**4. Eigenvalue decomposition (EVD) PCA** ::

To perform the PCA by eigenvalue decomposition

	``pca.py -t pca_test_trj.xtc -p complex.pdb -at CA -pt evd``

**Detailed usage:** ::

Run the following command to see the detailed usage and other options:

	``pca.py -h``


PCA on internal cordinates
-----------------------------

This programe performs the PCA on internal cordinates of a MD trajectory. User can select different types of internal cordinates such as:*pairwise distance between atoms*, *1-3 angle between backbone atoms*, *psi angle*, and *phi angle*. 

**1.PCA on pairwise distance between C-alpha atoms:** ::

To perform the PCA on pairwise distance between C-alpha atoms of pca_test_trc.xtc trajectory

	``internal_pca.py -t pca_test_trj.xtc -p complex.pdb -at CA -ct distance``	

**Detailed usage:** ::

Run the following command to see the detailed usage and other options:

	``internal_pca.py -h``

MDS  on MD trajectory
-------------------------------------------------

It perform the metric and non metric MDS (Multi-dimentional scaling) on a given MD trajectory. 
MDS is a tool to visualize the similarity or dissimilarity in a dataset. Two types of dissimilarity measures can be used in the case of a MD trajectory. First is Euclidean distance between internal cordinates of a protein structure, second is pairwise RMSD between a set of atoms over the frames of a MD trajectory.

**1. MDS on pairwise RMSD:** :: 

To perform the MDS on pairwise RMSD between C-alpha atoms of pca_test_trc.xtc trajectory
	
	``mds.py -t pca_test_trj.xtc -p complex.pdb -dt rmsd -ag CA``

**2. MDS on internal cordinates:** :: 

To perform the MDS on pairwise distance between C-alpha atoms of pca_test_trc.xtc trajectory. 

	``mds.py -t pca_test_trj.xtc -p complex.pdb -dt euc -ag CA``

**Detailed usage:** ::

Run the following command to see the detailed usage and other options:

	``mds.py -h``

t-SNE on MD trajectory
--------------------------------------------------------------------

t-SNE (t-distributed Stochastic Neighbor Embedding) is a tool for dimentionality reduction. 

**Detailed usage:** ::

Run the following command to see the detailed usage and other options:

	``tsne.py -h``

*Page created by: Bilal Nizami*
