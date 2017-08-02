#!/usr/bin/python
#filename: mds.py
import os, sys
import shlex, subprocess, time, re, traceback, math, matplotlib
import argparse
from time import sleep, gmtime, strftime
from datetime import datetime
import mdtraj as md
import numpy as np
from matplotlib import cm
from sklearn.metrics import euclidean_distances
from sklearn.manifold import TSNE
from sklearn import preprocessing
from write_plot import write_plots, write_pcs
from traj_info import trajectory_info, get_internal_cordinates
from itertools import combinations
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
print '\n\n'
print '|=======================================================|'
print '|\t\t\t\t\t\t\t|'
print '|\t :-) >>-----------> TSNE MD <-----------<< (-:	|'
print '|\t\t\t\t\t\t\t|'
print '|\t\t\t\t\t\t\t|'
print '|   This programe performs the TSNE (t-distributed  \t| \n|             Stochastic Neighbor Embedding.) on a MD trajectory\t\t|'
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


parser = argparse.ArgumentParser(usage='%(prog)s -t <MD trajectory> -p <topology file>')

parser.add_argument("-t", "--trj", dest="trj",				help="file name of the MD trajectory")
parser.add_argument("-p", "--top", dest="topology",				help="topology file")
parser.add_argument("-ct", "--cordinate_type",  dest="cordinate_type",				help="Type of cordinates to use for distance calculation")
parser.add_argument("-dt", "--dissimilarity_type",  dest="dissimilarity_type",				help="Type of dissimilarity matrix to use. Euclidean distance between internal cordinates or pairwise RMSD")



args = parser.parse_args()	


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

pca_traj = md.load(traj, top=topology)
top = pca_traj.topology
sele_grp=top.select("name CA")
#trajectory_info(pca_traj, traj, atm_name, sele_grp)

# pair wise RMSD 
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
calpha_idx=top.select_atom_indices('alpha')
atom_pairs = list(combinations(calpha_idx, 2)) # all unique pairs of elements 
#pairwise_distances = md.geometry.compute_distances(pca_traj, atom_pairs)
#int_cord=pairwise_distances

def tsne(input):
	't-distributed Stochastic Neighbor Embedding'
	seed = np.random.RandomState(seed=1)
	#similarities = euclidean_distances(input)
	my_tsne = TSNE(n_iter=3000, random_state=seed, init='pca')
	print "Performing TSNE.."
	mpos = my_tsne.fit_transform(input)
	write_plots('tsne_projection', mpos)
	return;

#int_cord=get_internal_cordinates()


if args.dissimilarity_type == 'rmsd' or args.dissimilarity_type == None:
	pair_rmsd=get_pair_rmsd(pca_traj, sele_grp)
	tsne(pair_rmsd)
if args.dissimilarity_type == 'distance':
	int_cord=get_internal_cordinates(top, args.cordinate_type, pca_traj)
	tsne(int_cord)
	


if __name__=="__main__":
	main()