#!/usr/bin/python
#filename: pca.py
from itertools import combinations
import mdtraj as md
import scipy.integrate
import numpy as np
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

if __name__=="__main__":
	main()