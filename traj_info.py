#!/usr/bin/python
#filename: pca.py
#==============================================================================#
#											
#			This programe is the part of PCA MD. It prints the information about trajectory 
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

if __name__=="__main__":
	main()