#!/usr/bin/python
#filename: pca.py
from itertools import combinations
import mdtraj as md
import scipy.integrate
import numpy as np
import math
from sklearn import preprocessing
from time import sleep, gmtime, strftime
#==============================================================================#
#											
#			This programe is the part of PCA MD. It read and prints the information about trajectory 
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#==============================================================================#
def man():

	return;

# print trajectory informations

def trajectory_info(pca_traj, traj, atm_name, sele_grp):
	'Prints various information of MD trajectory'
	print '\n\nTrajectory info:\n'
	print "Total",pca_traj.n_frames,"frames read from", traj
	print "MD time is from ", pca_traj.time[0],'to',pca_traj.time[-1],'ps'
	print pca_traj.n_atoms, "atoms and ", pca_traj.n_residues, "residues in the trajectory"
	print "Atom group selected for PCA:", atm_name
	print "Total", len(sele_grp), atm_name,'atoms selected for analysis\n'
	
	return;

#===========================================================
#
# Internal cordinate type
#
#===========================================================
def get_internal_cordinates(top, cordinate_type, pca_traj, atom_indices):
	'get the different types of internal cordinates as per user selections'
	calpha_idx=top.select_atom_indices(atom_indices)
	if cordinate_type == 'distance':
		print "Pair wise atomic distance selected\n "
		atom_pairs = list(combinations(calpha_idx, 2)) # all unique pairs of elements 
		pairwise_distances = md.geometry.compute_distances(pca_traj, atom_pairs)
		int_cord=pairwise_distances
		#print int_cord.shape
	if cordinate_type == 'phi':
		print  "phi torsions  selected\n"
		atom_pairs = list(combinations(calpha_idx, 3)) 
		angle=md.compute_phi(pca_traj)
		
		int_cord=angle[1] ## apparently compute_phi returns tupple of atoms indices and phi angles, index 1 has phi angles 
		#print np.array(angle[1]).shape
		#print int_cord[0]
	
	if cordinate_type == 'psi':
		print "psi torsions  selected\n"
		atom_pairs = list(combinations(calpha_idx, 3)) 
		angle=md.compute_psi(pca_traj)
		
		int_cord=angle[1] ## apparently compute_psi returns tupple of atoms indices and psi angles, index 1 has psi angles 
		#print np.array(angle[1]).shape
		
	if cordinate_type == 'angle':
		print "1-3 angle selected between C,CA and CB"
		cbeta_idx=top.select_atom_indices('minimal')
		#print cbeta_idx
		
	return int_cord;
	
#==========================================================================
#
#		selecting the atoms 
#
# User passes the arguements to select the subset of atoms for PCA
#===========================================================================


def get_trajectory(atm_name, top):
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
##============================================
#
#	cosine content 
#
#=================================================

def get_cosine(pca_sele_traj_reduced, pc_idx):
	i=pc_idx
	t = np.arange(len(pca_sele_traj_reduced))
	T = len(pca_sele_traj_reduced)
	cos = np.cos(np.pi * t * (i+1 ) / T)
	cos=((2.0 / T) * (scipy.integrate.simps(cos*pca_sele_traj_reduced[:, i])) ** 2/scipy.integrate.simps(pca_sele_traj_reduced[:, i] ** 2))
	#print cos
	return cos;

#=====================================
#
# Kaiser-Meyer-Olkein index
#
#==================================

def get_kmo(input):
	'compute the global KMO (kaiser-meyer-olkein) measure for the input dataset. KMO measure is the indication of sampling significace. KMO value range from 1 to 0, 1 indicating a perfect sampling.\
	adopted from book "Comparitive approach to using R and python for statistical data analysis"'
	
	# get the inverse of correlation / covariance matrix
	cov_inv = np.linalg.inv(input)
	cov_inv=np.absolute(cov_inv)
	input_nrow, input_ncol = input.shape

	# partial correlation matrix
	B = np.ones((input_nrow, input_ncol))
	for i in range(0,input_nrow,1):
		for j in range(i, input_ncol, 1):
			
			# above the diagonal
			B[i,j]=-(cov_inv[i,j])/(math.sqrt(cov_inv[i,i]*cov_inv[j,j]))
			
			# below the diagonal
			B[j,i]=B[i,j]
	
	# KMO 
	numerator_kmo=np.sum(np.square(input))-np.sum(np.square(np.diagonal(input)))
	deno_kmo=numerator_kmo+np.sum(np.square(B))-np.sum(np.square(np.diagonal(B)))
	kmo=numerator_kmo/deno_kmo
	print 'KMO for input trajectory is ', kmo
	return;
	
def print_kmo(pca_traj, traj, atm_name, sele_grp):
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	arr = sele_traj_reshaped
	#===============================================
	# covariance matrix 
	#===============================================
	cov_mat = np.cov(arr, rowvar=False)
	
	get_kmo(cov_mat)
	return;
	
def get_rmsf(pca_traj, sele_grp, trj_eval):
	'calculates the RMSD mode'
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	arr = sele_traj_reshaped
	
	iter=len(sele_grp)*3
	iter=iter
	nrow=pca_traj.n_frames
	ncol = len(sele_grp)

	# matrix for storing RMSF
	B = np.ones((nrow, ncol))
	C = np.ones((nrow, ncol))
	D = np.ones((nrow, ncol))
	#print len(trj_eval)
	for i in range(0,pca_traj.n_frames):
		k=0
		for j in range(0,iter,3):   ## iterate over every third coloumn (x,y,z of each atom)
			B[i,k]=math.sqrt(pow(sele_traj_reshaped[i,j],2)+pow(sele_traj_reshaped[i,j+1],2)+pow(sele_traj_reshaped[i,j+2],2)*trj_eval.real[0]) # for first eigenvector, multiply by first eigenvalue
			C[i,k]=math.sqrt(pow(sele_traj_reshaped[i,j],2)+pow(sele_traj_reshaped[i,j+1],2)+pow(sele_traj_reshaped[i,j+2],2)*trj_eval.real[1])
			D[i,k]=math.sqrt(pow(sele_traj_reshaped[i,j],2)+pow(sele_traj_reshaped[i,j+1],2)+pow(sele_traj_reshaped[i,j+2],2)*trj_eval.real[2])
			k=k+1
	
	## write the RMSF mode file for first PC
	
	np.savetxt('pc1_rmsf.agr',np.average(B, axis=0))
	rf = open('pc1_rmsf.agr', 'r')
	rf_cont = rf.read()
	rf.close()
	
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends = '@    title "RMSD modes"\n\
	@    xaxis  label "Residue"\n\
	@    yaxis  label "RMSD"\n\
	@	TYPE xy\n'
	
	pf = open('pc1_rmsf.agr', 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+rf_cont)
	pf.close()

	## write the RMSF mode file for second PC
	np.savetxt('pc2_rmsf.agr',np.average(C, axis=0))
	rf = open('pc2_rmsf.agr', 'r')
	rf_cont = rf.read()
	rf.close()
	
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends = '@    title "RMSD modes"\n\
	@    xaxis  label "Residue"\n\
	@    yaxis  label "RMSD"\n\
	@	TYPE xy\n'
	
	pf = open('pc2_rmsf.agr', 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+rf_cont)
	pf.close()
	
	## write the RMSF mode file for third PC
	np.savetxt('pc3_rmsf.agr',np.average(D, axis=0))
	rf = open('pc3_rmsf.agr', 'r')
	rf_cont = rf.read()
	rf.close()
	
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends = '@    title "RMSD modes"\n\
	@    xaxis  label "Residue"\n\
	@    yaxis  label "RMSD"\n\
	@	TYPE xy\n'
	
	pf = open('pc3_rmsf.agr', 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+rf_cont)
	pf.close()
	
	return;
	
	
	
if __name__=="__main__":
	main()