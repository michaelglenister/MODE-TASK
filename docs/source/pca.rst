Principle Component Analysis
====================================

PCA on cartesian cordinates
-----------------------------

This programe performs the PCA (Principal Component Analysis) on a MD trajectory	

Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >

Options:
  -h, --help            show this help message and exit
  -t TRJ, --trj=TRJ     file name of the MD trajectory
  -p TOPOLOGY, --top=TOPOLOGY
                        topology file
  -a ATM_GRP, --ag=ATM_GRP
                        group of atom for PCA. Default is C alpha atoms. Other
                        options are :all= all atoms, backbone = backbone
                        atoms, CA= C alpha atoms, protein= protein's atoms
  -r REFERENCE, --ref=REFERENCE
                        reference structure for RMSD
  -m PCA_TYPE, --pca_type=PCA_TYPE
                        PCA method. Default is normal PCA. Options are:
                        KernelPCA, normal, ipca. If normal is selected,
                        additional arguments can be passed by flag -svd. If
                        KernelPCA is selected kernel type can also be defined
                        by flag -k
  -k KERNEL_TYPE, --kernel_type=KERNEL_TYPE
                        Type of kernel for KernalPCA. default is linear.
                        Options are :linear, poly, rbf, sigmoid, cosine,
                        precomputed
  -s SVD_SOLVER, --svd_solver=SVD_SOLVER
                        Type of svd_solver for normal PCA. Default is auto.
                        Options are :auto, full, arpack, randomized


PCA on internal cordinates
-----------------------------

This programe performs the PCA on internal cordinates of a MD trajectory	

Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >

Options:
  -h, --help            show this help message and exit
  -t TRJ, --trj=TRJ     file name of the MD trajectory
  -p TOPOLOGY, --top=TOPOLOGY
                        topology file
  -a ATM_GRP, --ag=ATM_GRP
                        group of atom for PCA. Default is C alpha atoms. Other
                        options are :all= all atoms, backbone = backbone
                        atoms, CA= C alpha atoms, protein= protein's atoms
  -c CORDINATE_TYPE, --ct=CORDINATE_TYPE
                        Internal cordinate type. Options are: distance,
                        angles, dihedral

