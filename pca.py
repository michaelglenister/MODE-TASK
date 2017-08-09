#!/usr/bin/python
#filename: pca.py
import os, sys
import time, re, math
import argparse
from time import sleep, gmtime, strftime
from datetime import datetime
import mdtraj as md
import numpy as np
from sklearn.decomposition import PCA, KernelPCA, IncrementalPCA
from sklearn.metrics import euclidean_distances
from sklearn import preprocessing
from write_plot import write_plots, write_pcs
from traj_info import trajectory_info, get_cosine, get_kmo, print_kmo
from welcome_msg import welcome_msg
import scipy.integrate
def main():
	
	return;


#==============================================================================#
#											PCA MD 
#
#			This programe performs the PCA on a MD trajectory
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#==============================================================================#

##===============================================================================
##								 Welcome message
##===============================================================================
title='PCA MD'
welcome_msg(title)


#==============================================================================
#                            Setting the options
#==============================================================================

def set_option():
	parser = argparse.ArgumentParser( usage='%(prog)s -t <MD trajectory> -p <topology file>')
	#"Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >"
	
	parser.add_argument("-t", "--trj", dest="trj", help="file name of the MD trajectory", action="store")
	parser.add_argument("-p", "--top", dest="topology", help="topology file")      
	parser.add_argument("-ag", "--ag", dest="atm_grp", help="group of atom for PCA. Default is C alpha atoms. Other options are :"				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")	
	parser.add_argument("-r", "--ref", dest="reference", help="reference structure for RMSD") 
	parser.add_argument("-pt", "--pca_type", dest="pca_type", help="PCA method. Default is svd (Single Value Decomposition) PCA. Options are:\
					evd, kpca, svd, ipca. If svd is selected, additional arguments can be passed by flag -st. If KernelPCA is selected kernel type can also be defined by flag -kt") 	
	parser.add_argument("-nc", "--comp", type=int, dest="comp", help="Number of components to keep in a PCA object. If not set, by default all the components will be kept.")	
	parser.add_argument("-kt", "--kernel_type", dest="kernel_type", help="Type of kernel for KernalPCA. default is linear. Options are :"
					"linear, poly, rbf, sigmoid, cosine, precomputed") 
	parser.add_argument("-st", "--svd_solver", dest="svd_solver", help="Type of svd_solver for SVD (Single Value Decomposition) PCA. Default is auto. Options are :"				  "auto, full, arpack, randomized") 
	args = parser.parse_args()	
	
	
	#====================================================================
	# if no arguments are passed
	#====================================================================
	if args.trj is None: 
		print 'ERROR: Missing trajectory argument.... :(  \nPlease see the help by running \n\nsystem_setup.py -h\n\n '
		parser.print_help()
		sys.exit(1)
	
	if args.topology is None:
		print 'ERROR: Missing topology.... :( \nPlease see the help by running \n\nsystem_setup.py -h\n\n '
		parser.print_help()
		sys.exit(1)
	
	if not os.path.exists(args.trj ):
		print('\nERROR: {0} not found....:(  Please check the path\n' .format(args.trj ))
		parser.print_help()
		sys.exit(1)
	
	if not os.path.exists(args.topology):
		print('\nERROR: {0} not found....:(  Please check the path\n' .format(args.topology ))
		parser.print_help()
		sys.exit(1)
	if args.pca_type not in  ('svd', 'evd', 'kpca', 'ipca', None):
		print 'ERROR: no such option as', args.pca_type, 'for flag -pt \nPlease see the help by running \n pca.py -h..\n\n '
		sys.exit(1)
	if args.kernel_type not in  ('linear', 'poly', 'sigmoid', 'cosine', 'precomputed', 'rbf',None):
		print 'ERROR: no such option as', args.kernel_type, 'for flag -kt \nPlease see the help by running \n pca.py -h..\n\n '
		sys.exit(1)
		
	if args.kernel_type != None and args.pca_type != 'kpca':
		print 'WARNING: -kt', args.kernel_type, 'is meaningless with -pt', args.pca_type, '. Flag -kt is being ignored!'
		
	if args.svd_solver not in  ('auto', 'full', 'arpack', 'randomized', None):
		print 'ERROR: no such option as', args.svd_solver, 'for flag -st \nPlease see the help by running \n pca.py -h.\n\n'
		sys.exit(1)
	
	if args.svd_solver != None and args.pca_type != 'svd':
		print 'WARNING: -st', args.svd_solver, 'is meaningless with -pt', args.pca_type, '. Flag -st is being ignored!'
	return args
	
args = set_option()
atm_name = args.atm_grp
#=======================================
# assign the passed arguments and read the trajectory 
#=======================================
traj = args.trj
topology = args.topology
ref = args.reference
ptype=args.pca_type
comp = args.comp
n_eivec=5
## read the reference structure
if ref:
	try:
		ref = md.load(args.reference)
	except:
			raise IOError('Could not open reference structure {0} for reading. \n' .format(args.reference))
	
try:
	pca_traj = md.load(traj, top=topology)
except:
	raise IOError('Could not open trajectory {0} for reading. \n' .format(trj))
	
top = pca_traj.topology


#==============================================
#
# Setting the default options
#
#===============================================

if args.atm_grp == None:
	print 'No atom has been selected. PCA will be performed on C alpha atoms '
	atm_name = 'CA'  # set to default C-alpha atoms
if args.reference == None:
	print "No reference structure given, RMSD will be computed to the first frame in the trajectory"
	ref = pca_traj # set reference to current trajectory
if args.pca_type == None:
	ptype = 'svd'

if args.svd_solver == None:
	svd='auto'
	
#==========================================================================
#
#		selecting the atoms 
#
# User passes the arguements to select the subset of atoms for PCA
#===========================================================================


def get_trajectory():
	'get the part of system for PCA based on users input of atom group'
	if atm_name == 'CA':
		sele_grp=top.select("name CA")	

	if atm_name == 'backbone':
		sele_grp=top.select("backbone")
		
	if atm_name == 'all':
		sele_grp=top.select("all")
	return sele_grp;
	
	if atm_name == 'protein':
		sele_grp=top.select("protein")
	return sele_grp;

sele_grp = get_trajectory()

# Print trajectory information
trajectory_info(pca_traj, traj, atm_name, sele_grp)

# print KMO 
print_kmo(pca_traj, traj, atm_name, sele_grp)

#===============================================================
#
#  RMSD in reference with first frame
#
#===============================================================

def get_rmsd():
	rmsd = md.rmsd(pca_traj, ref, 0, atom_indices=sele_grp)
	print "RMSD written to rmsd.agr \n "
	## write the RMSD file
	np.savetxt('rmsd.agr', rmsd)
	rf = open('rmsd.agr', 'r')
	rf_cont = rf.read()
	rf.close()
	
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends = '@    title "RMSD"\n\
	@    xaxis  label "Time"\n\
	@    yaxis  label "RMSD"\n\
	@	TYPE xy\n'
	
	pf = open('rmsd.agr', 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+rf_cont)
	pf.close()
	return;
	
get_rmsd()

#===============================================================
#
# PCA using sci-kit learn library
#===============================================================

def svd_pca(svd):
	'single value decomposition based PCA'
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped, axis=0, with_std=False) # center to the mean
	
	pca_sele_traj = PCA(n_components=comp)
	pca_sele_traj.fit(sele_traj_reshaped_scaled)
	pca_sele_traj_reduced = pca_sele_traj.transform(sele_traj_reshaped_scaled)
	#pca_sele_traj_reduced=pca_sele_traj_reduced[:,0:n_pc]
	np.savetxt('pca_sele_traj_reduced.txt', pca_sele_traj_reduced)
	print "Trace of the covariance matrix is: ", np.trace(pca_sele_traj.get_covariance())
	print "Wrote covariance matrix..."
	np.savetxt('cov.dat', pca_sele_traj.get_covariance())
	
	# write the plots 
	write_plots('pca_projection', pca_sele_traj_reduced)
	#write the pcs variance
	#print type(pca_sele_traj.explained_variance_ratio_)
	write_pcs('pca_variance', pca_sele_traj)
	
	pc1_cos=get_cosine(pca_sele_traj_reduced, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(pca_sele_traj_reduced, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(pca_sele_traj_reduced, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(pca_sele_traj_reduced, 3)
	print 'cosine content of 4th PC=', pc4_cos
	
	return;


#==============================================================
#
# Kernel PCA
#
# ==============================================================
def my_kernelPCA(kernel):
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped, axis=0, with_std=False) # center to the mean

	kpca = KernelPCA(kernel = kernel, fit_inverse_transform=True, gamma=10)
	kpca.fit(sele_traj_reshaped_scaled)
	#print "Trace of the covariance matrix is: ", np.trace(kpca.get_covariance())
	kpca_reduced = kpca.transform(sele_traj_reshaped_scaled)
	
	#write plots
	write_plots('kpca_projection', kpca_reduced)
	
	#write variance
	np.savetxt('kpca_variance', kpca.lambdas_)
	
	pc1_cos=get_cosine(kpca_reduced, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(kpca_reduced, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(kpca_reduced, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(kpca_reduced, 3)
	print 'cosine content of 4th PC=', pc4_cos
	return;


#=============================================================
#
# Incremental PCA
#
#=============================================================

def incremental_pca():
	' normal PCA is very memory intesive. It can be problemetic for large dataset, \
	since dataset is stored in memory. Incremental principal component analysis (IPCA) is \
	typically used for such cases. '
	
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped, axis=0, with_std=False) # center to the mean

	ipca = IncrementalPCA()
	ipca = ipca.fit(sele_traj_reshaped_scaled)
	ipca_reduced=ipca.transform(sele_traj_reshaped_scaled)
	
	#write plots
	write_plots('ipca_projection', ipca_reduced)
	
	#write variance
	#np.savetxt('ipca_variance', kpca.lambdas_)
	pc1_cos=get_cosine(ipca_reduced, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(ipca_reduced, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(ipca_reduced, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(ipca_reduced, 3)
	print 'cosine content of 4th PC=', pc4_cos

	return;

#===============================================================
#
#  Eigenvalue decomposition based PCA
#
#=================================================================
def my_pca():
	'eigenvales decomposition PCA'
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped, axis=0, with_std=False) # center to the mean
	arr = sele_traj_reshaped_scaled
	
	
	#===============================================
	# covariance matrix 
	cov_mat = np.cov(arr, rowvar=False)
	trj_eval, trj_evec=np.linalg.eig(cov_mat)
	
	print "Trace of cov matrix is ",  np.trace(cov_mat)
	
	#=============================
	# sanity check of calculated eigenvector and eigen values 
	# it must be cov matrix * eigen vector = eigen vector * eigen value
	
	for i in range(len(trj_eval)):
		eigv = trj_evec.real[:,i].reshape(1,len(trj_evec[:,0]),).T
		np.testing.assert_array_almost_equal(cov_mat.dot(eigv), trj_eval[i]*eigv, decimal=3, err_msg='', verbose=True)

#=============================================
	# sort the eigenvalues and eigenvector
	sort_idx = trj_eval.argsort()[::-1]
	trj_eval = trj_eval[sort_idx]
	trj_evec = trj_evec[sort_idx]
	
	# join first two eigenvector into a single matrix
	#eivec_1 = trj_evec.real[:,0].reshape(len(trj_evec[:,0]),1)
	#eivec_2 = trj_evec.real[:,1].reshape(len(trj_evec[:,1]),1)
	#pca = np.concatenate((eivec_1, eivec_2, eivec_3, eivec_4, eivec_5), axis=1)	

	tot_var = np.sum(trj_eval.real)
	variation = []
	cum = []
	j = 0
	eigv = []
	#i=0
	n_comp=100
	pca = trj_evec.real[:,0:n_comp]    ## keep first 100 eigenvectors
	for i in trj_eval.real[0:n_comp]:
		eigv.append(i)
		variation.append(i/tot_var)
		j +=1
	#write_plots('variation', variation)
	#========================================================
	# transform the input data into choosen pc
	arr_transformed = pca.T.dot(arr.T)
	#arr_transformed = np.concatenate((arr_transformed[0,:].reshape(len(arr_transformed[0,:]),1), arr_transformed[1,:].reshape(len(arr_transformed[1,:]),1)), axis=1)
	#arr_transformed[:,0]
	write_plots('pca_projection', arr_transformed)
	
	pc1_cos=get_cosine(arr_transformed, 0)
	print 'cosine content of first PC=',pc1_cos
	pc2_cos=get_cosine(arr_transformed, 1)
	print 'cosine content of second PC=', pc2_cos
	pc3_cos=get_cosine(arr_transformed, 2)
	print 'cosine content of 3rd PC=',pc3_cos
	pc4_cos=get_cosine(arr_transformed, 3)
	print 'cosine content of 4th PC=', pc4_cos

	return;



if ptype == 'kpca':
	kernel = ''
	kernel = args.kernel_type
	if args.kernel_type:
		print "Performing Kernel PCA with", kernel, 'kernel'
		my_kernelPCA(kernel)
		print "\nFINISHED. !"
	else:
		print "Performing Kernel PCA with default linear kernel"
		my_kernelPCA('linear')
		print "\nFINISHED. !"

if ptype == 'svd':
	svd=''
	svd = args.svd_solver
	if svd:
		print "Performing SVD (Single Value Decomposition) PCA with ",svd,"svd_solver"
		svd_pca(svd)
		print "\nFINISHED. !"
	else:
		print "Performing SVD (Single Value Decomposition) PCA with 'auto' svd_solver"
		svd_pca(svd)
		print "\nFINISHED. !"
if ptype == 'ipca':
	print "Performing Incremental_pca (IPCA)"
	incremental_pca()
	print "\nFINISHED. !"

if ptype == 'evd':
	my_pca()
	print "\nFINISHED. !"

if __name__=="__main__":
	main()