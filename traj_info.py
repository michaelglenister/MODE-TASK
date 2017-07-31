#!/usr/bin/python
#filename: pca.py
from itertools import combinations
import mdtraj as md
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
def get_internal_cordinates(top, cordinate_type, pca_traj):
	'get the different types of internal cordinates as per user selections'
	calpha_idx=top.select_atom_indices('alpha')
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

if __name__=="__main__":
	main()