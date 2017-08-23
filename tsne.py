#!/usr/bin/env python
#filename: tsne.py
import os, sys
import argparse
import mdtraj as md
import numpy as np
from sklearn.metrics import euclidean_distances
from sklearn.manifold import TSNE
from write_plot import write_plots, write_pcs, write_fig
from traj_info import trajectory_info, get_internal_cordinates, get_trajectory
from welcome_msg import welcome_msg

def main():
	
	return;

#==============================================================================#
#											TSNE MD
#
#			This programe performs the TSNE on a MD trajectory
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#==============================================================================#

##===============================================================================
##								 Welcome message
##===============================================================================
title='TSNE MD'
welcome_msg(title)
#==============================================================================
#                            Setting the options
#==============================================================================

def get_options():
	
	parser = argparse.ArgumentParser(usage='%(prog)s -t <MD trajectory> -p <topology file>')
	
	parser.add_argument("-t", "--trj", dest="trj",				help="file name of the MD trajectory")
	parser.add_argument("-p", "--top", dest="topology",				help="topology file")
	parser.add_argument("-out", "--out", dest="out_dir", help="Name of the output directory. Default is out")
	parser.add_argument("-ag", "--ag", dest="atm_grp", help="group of atom for MDS. Default is C alpha atoms. Other options are :"				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")	
	parser.add_argument("-ct", "--cordinate_type",  dest="cordinate_type",				help="Type of cordinates to use for distance calculation")
	parser.add_argument("-dt", "--dissimilarity_type",  dest="dissimilarity_type",				help="Type of dissimilarity matrix to use. Euclidean distance between internal cordinates or pairwise RMSD")
	args = parser.parse_args()
	
	if args.out_dir == None:
		out=args.trj
		args.out_dir=out
		
	if args.atm_grp == None:
		print 'No atom selected. MDS will be performed on C alpha atoms '
		args.atm_grp = 'CA'  # set to default C-alpha atoms
		
	if args.atm_grp not in  ('all', 'CA', 'backbone', 'protein'):
		print 'ERROR: no such option as', args.atm_grp, 'for flag -at \nPlease see the usage\n\n '
		sys.exit(1)
	return args;
	
args=get_options()

	


#====================================================================
# if no arguments are passed
#====================================================================
if args.trj is None:
	print 'Missing trajectory arguments :(\nPlease see the help by running \n\nsystem_setup.py -h\n\n '
	parser.print_help()
	sys.exit(1)

if args.topology is None:
	print 'Missing topology !!\nPlease see the help by running \n\nsystem_setup.py -h\n\n '
	parser.print_help()
	sys.exit(1)
if args.cordinate_type == None:
	args.cordinate_type = "distance"
	print "distance is the cordinate type by default"

#=======================================
# assign the passed arguments and read the trajectory
#=======================================

traj = args.trj
topology = args.topology

#pca_traj = md.load(traj, top=topology)

print 'Reading trajectory ', args.trj, '...' 

try:
	pca_traj = md.load(traj, top=topology)
except:
	raise IOError('Could not open trajectory {0} for reading. \n' .format(trj))
top = pca_traj.topology
#sele_grp=top.select("name CA")
#atm_name='CA'
atm_name=args.atm_grp
sele_grp = get_trajectory(atm_name, top)

# take the input trj name for output directory
out_dir=args.out_dir
out_dir=out_dir.split('/')
out_dir=out_dir[-1]
out_dir='out_'+out_dir

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
#else:
#	print out_dir, 'already exist. Can not overwrite the output directory!\n'
#	sys.exit(1)
print 'Results will be written in ', out_dir

## =============================
# trajectory info
#===================================
trajectory_info(pca_traj, traj, atm_name, sele_grp)
# ==================================================
#	pair wise RMSD 
#
#===================================================

def get_pair_rmsd(pca_traj, sele_grp):
	'pair wise RMSD over all the frames, return a square matrix of pairwise rmsd'
	pair_rmsd=np.empty((pca_traj.n_frames, pca_traj.n_frames))
	for i in range(pca_traj.n_frames):
		pair_rmsd[i]=md.rmsd(pca_traj, pca_traj, i, atom_indices=sele_grp)
	pair_rmsd=(pair_rmsd+pair_rmsd.transpose())/2  ## due to precision level matrix might not evaluate as symmetric, hence to make it symmetric
	return pair_rmsd;


#============================================
#
#  TSNE
#
#=============================================

def tsne(input):
	't-distributed Stochastic Neighbor Embedding'
	seed = np.random.RandomState(seed=1)
	my_tsne = TSNE(n_components=10, n_iter=3000, random_state=seed, init='pca')
	print "Performing TSNE..."
	mpos = my_tsne.fit_transform(input)
	#print mpos.shape
	write_plots('tsne_projection', mpos, out_dir)
	write_fig('tsne_projection', mpos, out_dir)

	return;


if args.dissimilarity_type == 'rmsd' or args.dissimilarity_type == None:
	pair_rmsd=get_pair_rmsd(pca_traj, sele_grp)
	tsne(pair_rmsd)
	print 'FINISHED!'
if args.dissimilarity_type == 'distance':
	int_cord=get_internal_cordinates(top, args.cordinate_type, pca_traj)
	tsne(int_cord)
	print 'FINISHED!'
	


if __name__=="__main__":
	main()