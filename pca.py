#!/usr/bin/python
#filename: pca.py
import os, sys
import shlex, subprocess, time, re, argparse, traceback, math, matplotlib
import argparse
from time import sleep, gmtime, strftime
from datetime import datetime
import mdtraj as md
import numpy as np
from matplotlib import cm
matplotlib.use('Agg')
from sklearn.decomposition import PCA, KernelPCA, IncrementalPCA
from sklearn.metrics import euclidean_distances
from sklearn.manifold import MDS
from sklearn import preprocessing
from itertools import combinations
from write_plot import write_plots, write_pcs


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
print '\n\n'
print '|=======================================================|'
print '|\t\t\t\t\t\t\t|'
print '|\t :-) >>-----------> PCA MD <-----------<< (-:	|'
print '|\t\t\t\t\t\t\t|'
print '|\t\t\t\t\t\t\t|'
print '|   This programe performs the PCA (Principal Component\t| \n|             Analysis) on a MD trajectory\t\t|'
print '|\t\t\t\t\t\t\t|\n', '|\tAuthors:  Bilal Nizami\t\t\t\t|\n','|\tResearch Unit in Bioinformatics (RUBi)\t\t|\n', '|\tRhodes University, 2017\t\t\t\t|'
print '|\tDistributed under GNU GPL 3.0\t\t\t|'
print '|\t\t\t\t\t\t\t|'
print '|\thttps://github.com/michaelglenister/NMA-TASK\t|'
print '|\t\t\t\t\t\t\t|'
print '|=======================================================|'
print '\n'

#==============================================================================
#                            Setting the options
#==============================================================================

parser = argparse.ArgumentParser( usage='%(prog)s -t <MD trajectory> -p <topology file>')
#"Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >"

parser.add_argument("-t", "--trj", dest="trj", help="file name of the MD trajectory", action="store")

parser.add_argument("-p", "--top", dest="topology", help="topology file")      

parser.add_argument("-at", "--ag", dest="atm_grp", help="group of atom for PCA. Default is C alpha atoms. Other options are :"				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")	

parser.add_argument("-r", "--ref", dest="reference", help="reference structure for RMSD") 

parser.add_argument("-m", "--pca_type", dest="pca_type", help="PCA method. Default is svd (Single Value Decomposition) PCA. Options are:\
				  KernelPCA, svd, ipca. If svd is selected, additional arguments can be passed by flag -s. If KernelPCA is selected kernel type can also be defined by flag -k") 
				  
parser.add_argument("-nc", "--comp", type=int, dest="comp", help="Number of components to keep in a PCA object. If not set, by default all the components will be kept.")
				  
parser.add_argument("-kt", "--kernel_type", dest="kernel_type", help="Type of kernel for KernalPCA. default is linear. Options are :"
				  "linear, poly, rbf, sigmoid, cosine, precomputed") 

parser.add_argument("-st", "--svd_solver", dest="svd_solver", help="Type of svd_solver for SVD (Single Value Decomposition) PCA. Default is auto. Options are :"				  "auto, full, arpack, randomized") 
				  

args = parser.parse_args()
                     
atm_name = args.atm_grp



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
#=======================================
# assign the passed arguments and read the trajectory 
#=======================================

traj = args.trj
topology = args.topology
ref = args.reference
ptype=args.pca_type
comp = args.comp
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
	
#print(pca_traj)
top = pca_traj.topology
#print top

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

# print trajectory informations
def trajectory_info():
	'Prints various information of MD trajectory'
	print '\n\nTrajectory info:\n'
	print "Total",pca_traj.n_frames,"frames read from", traj
	print "MD time is from ", pca_traj.time[0],'to',pca_traj.time[-1],'ps'
	print pca_traj.n_atoms, "atoms and ", pca_traj.n_residues, "residues in the trajectory"
	print "Atom group selected for PCA:", atm_name, "\n"
	
	if args.reference == None: 
		print "Reference for RMSD calculation is first frame of trajectory"
	else:
		print "Reference for RMSD calculation is: ", args.reference
		
	print "Total", len(sele_grp), atm_name,'atoms selected for analysis\n'
	
	return;

	
trajectory_info()

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
	pca_sele_traj = PCA(n_components=2)
	pca_sele_traj.fit(sele_traj_reshaped)
	pca_sele_traj_reduced = pca_sele_traj.transform(sele_traj_reshaped)
	#print pca_sele_traj_reduced
	print "Trace of the covariance matrix is: ", np.trace(pca_sele_traj.get_covariance())
	print "Wrote covariance matrix..."
	np.savetxt('cov.dat', pca_sele_traj.get_covariance())
		
	# write the plots 
	write_plots('pca_projection', pca_sele_traj_reduced)
	#write the pcs variance
	#print type(pca_sele_traj.explained_variance_ratio_)
	write_pcs('pca_variance', pca_sele_traj)
	
	# print spinner ..need to fix
	print "processing...\\",
	syms = ['\\', '|', '/', '-']
	bs = '\b'
	for _ in range(10):
		for sym in syms:
			sys.stdout.write("\b%s" % sym)
			sys.stdout.flush()
			
			
			time.sleep(.5)
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
	
	kpca = KernelPCA(kernel = kernel, fit_inverse_transform=True, gamma=10)
	kpca.fit(sele_traj_reshaped)
	#print "Trace of the covariance matrix is: ", np.trace(kpca.get_covariance())
	kpca_reduced = kpca.transform(sele_traj_reshaped)
	
	#write plots
	write_plots('kpca_projection', kpca_reduced)
	
	#write variance
	np.savetxt('kpca_variance', kpca.lambdas_)
	return;


#=============================================================
#
# Incremental PCA
#
#=============================================================

def incremental_pca():
	' normal PCA is not very memory intesive. It can be problemetic for large dataset, \
	since dataset is stored in memory. Incremental principal component analysis (IPCA) is \
	typically used for such cases. '
	
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	
	ipca = IncrementalPCA()
	ipca = ipca.fit(sele_traj_reshaped)
	ipca_reduced=ipca.transform(sele_traj_reshaped)
	
	#write plots
	write_plots('ipca_projection', ipca_reduced)
	
	#write variance
	#np.savetxt('ipca_variance', kpca.lambdas_)

	return;

#===============================================================
#
#  My own method
#
#=================================================================
def my_pca():
	
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
		
	sele_traj_reshaped = sele_traj_reshaped.astype(float) ## to avoid numpy Conversion Error during scaling
	sele_traj_reshaped_scaled = preprocessing.scale(sele_traj_reshaped)
	arr = sele_traj_reshaped_scaled
	
	
	#===============================================
	# covariance matrix of selected coloumns
	cov_mat = np.cov(arr, rowvar=False)
	cov_mat = np.cov(cov_mat)
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
	
	pca = np.empty_like(trj_evec)
	
	# join first two eigenvector into a single matrix and write plot
	eivec_1 = trj_evec.real[:,0].reshape(len(trj_evec[:,0]),1)
	eivec_2 = trj_evec.real[:,1].reshape(len(trj_evec[:,1]),1)
	pca = np.concatenate((eivec_1, eivec_2), axis=1)

	write_plots('my_method_proj', pca)
	

#=============================================
# sort the eigenvales
	e_p = []
	print type(e_p)
	for i in range(len(trj_eval)):
		eig_pairs = [np.abs(trj_eval[i]), trj_evec[:,i]]
		e_p.append(eig_pairs)
	#print (e_p[2])
	e_p.sort(key=lambda x: x[0], reverse=True)
	
	
	# sorted eigenvalues and variation explained
	print ('sorted eigenvalues')
	tot_var = 0
	for i in e_p:
		tot_var +=i[0]
	variation = []
	print tot_var
	cum = []
	j = 0
	eigv = []
	for i in e_p:
		#print (i[0])
		eigv.append(i[0])
		variation.append(i[0]/tot_var)
		#print ("variation explained:",variation[j]*100)
		cum = np.cumsum(variation)
		#print ('cumulative: ', cum[j]*100 )
		j +=1
	return;


if ptype == 'KernelPCA':
	kernel = ''
	kernel = args.kernel_type
	if args.kernel_type:
		print "Performing Kernel PCA with", kernel, 'kernel'
		my_kernelPCA(kernel)
	else:
		print "Performing Kernel PCA with default linear kernel"
		my_kernelPCA('linear')

if ptype == 'svd':
	svd=''
	svd = args.svd_solver
	if svd:
		print "Performing SVD (Single Value Decomposition) PCA with ",svd,"svd_solver"
		svd_pca(svd)
	else:
		print "Performing SVD (Single Value Decomposition) PCA with 'auto' svd_solver"
		svd_pca(svd)
if ptype == 'ipca':
	print "Performing Incremental_pca (IPCA)"
	incremental_pca()

if ptype == 'my':
	my_pca()

