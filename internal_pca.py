#!/usr/bin/python
#filename: pca.py

import os, sys
import shlex, subprocess, time, re, argparse, traceback, math, matplotlib
from optparse import OptionParser
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
#											Internal PCA MD 
#
#			This programe performs the PCA on internal cordinates of a MD trajectory
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
print '|\t :-) >>-------> Internal PCA MD <-------<< (-:	|'
print '|\t\t\t\t\t\t\t|'
print '|\t\t\t\t\t\t\t|'
print '|   This programe performs the PCA \t| \n|  on internal cordinates of a MD trajectory\t\t|'
print '|\t\t\t\t\t\t\t|\n', '|\tAuthors:  Bilal Nizami\t\t\t\t|\n','|\tResearch Unit in Bioinformatics (RUBi)\t\t|\n', '|\tRhodes University, 2017\t\t\t\t|'
print '|\tDistributed under GNU GPL 3.0\t\t\t|'
print '|\t\t\t\t\t\t\t|'
print '|\thttps://github.com/michaelglenister/NMA-TASK\t|'
print '|\t\t\t\t\t\t\t|'
print '|=======================================================|'
print '\n'

def set_option():
	parser = argparse.ArgumentParser( usage='%(prog)s -t <MD trajectory> -p <topology file>')
	#"Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >"
	
	parser.add_argument("-t", "--trj", dest="trj", help="file name of the MD trajectory", action="store")
	parser.add_argument("-p", "--top", dest="topology", help="topology file")      
	parser.add_argument("-at", "--ag", dest="atm_grp", help="group of atom for PCA. Default is C alpha atoms. Other options are :"				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")	
	parser.add_argument("-ct", "--ref", dest="cordinate_type", help="nternal cordinate type. Options are: distance, angles, dihedral") 
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
		
	return args
	
args = set_option()
#==============================================================================
#                            Setting the options
#==============================================================================

#parser = OptionParser("Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >")
#
#parser.add_option("-t", "--trj", type='string', dest="trj",
#                  help="file name of the MD trajectory")
#
#parser.add_option("-p", "--top", type='string', dest="topology",
#                  help="topology file")      
#
#parser.add_option("-a", "--ag", type='string', dest="atm_grp",
#                  help="group of atom for PCA. Default is C alpha atoms. Other options are :"
#				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")
#
#parser.add_option("-c", "--ct", type='string', dest="cordinate_type",
#                  help="Internal cordinate type. Options are: distance, angles, dihedral")
#
#(options, args) = parser.parse_args()
#                     
#
#atm_name = options.atm_grp

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

#=======================================
# assign the passed arguments and read the trajectory 
#=======================================

traj = args.trj
topology = args.topology
pca_traj = md.load(traj, top=topology)
top = pca_traj.topology

#==============================================
#
# Setting the default options
#
#===============================================

if args.atm_grp == None:
	print 'No atom has been selected. PCA will be performed on C alpha atoms '
	atm_name = 'CA'  # set to default C-alpha atoms


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
	
		
	print "Total", len(sele_grp), atm_name,'atoms selected for analysis\n'
	
	return;

	
trajectory_info()


#===========================================================
#
# Internal cordinate type
#
#===========================================================
def get_internal_cordinates():
	'get the different types of internal cordinates as per user selections'
	calpha_idx=top.select_atom_indices('alpha')
	if args.cordinate_type == 'distance':
		print "Pair wise atomic distance selected\n "
		atom_pairs = list(combinations(calpha_idx, 2)) # all unique pairs of elements 
		pairwise_distances = md.geometry.compute_distances(pca_traj, atom_pairs)
		int_cord=pairwise_distances
		#print int_cord.shape
	if args.cordinate_type == 'phi':
		print  "phi torsions  selected\n"
		atom_pairs = list(combinations(calpha_idx, 3)) 
		angle=md.compute_phi(pca_traj)
		
		int_cord=angle[1] ## apparently compute_phi returns tupple of atoms indices and phi angles, index 1 has phi angles 
		print np.array(angle[1]).shape
		#print int_cord[0]
	
	if args.cordinate_type == 'psi':
		print "psi torsions  selected\n"
		atom_pairs = list(combinations(calpha_idx, 3)) 
		angle=md.compute_psi(pca_traj)
		
		int_cord=angle[1] ## apparently compute_psi returns tupple of atoms indices and psi angles, index 1 has psi angles 
		print np.array(angle[1]).shape
		
	if args.cordinate_type == 'angle':
		print "1-3 angle selected between C,CA and CB"
		cbeta_idx=top.select_atom_indices('minimal')
		print cbeta_idx
		
	return int_cord;

#===========================================================
#
#  Internal Distance Coordinate Based PCA
#
#===========================================================
def distance_pca(int_cord1):
	'Internal Coordinate Based PCA'
	
	pca = PCA()
	
	dpca = pca.fit(int_cord1)
	dpca_reduced=dpca.transform(int_cord1)
	write_plots('dpca_projection', dpca_reduced)
	write_pcs('dpca_pcs', dpca)
	return;

int_cord=get_internal_cordinates()
distance_pca(int_cord)
