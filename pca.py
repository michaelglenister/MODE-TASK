<<<<<<< HEAD
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
from sklearn.decomposition import PCA
from sklearn import preprocessing


#==============================================================================#
#											PCA MD 
#
#			This programe performs the PCA on a MD trajectory
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
=======
#!/usr/bin/env python
#filename: pca.py
import os, sys
import shlex
import subprocess
import time
import re
from optparse import OptionParser
from time import sleep
from datetime import datetime
import mdtraj as md
import numpy as np
#from lib.utils import *
import argparse, traceback, math, matplotlib
from matplotlib import cm
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from time import gmtime, strftime
#from calc_correlation import parse_traj

#==============================================================================#
#							PCA MD 
#
#		This programe performs the PCA on a MD trajectory
#
# 							Author : Bilal Nizami
# 						   Rhodes University, 2017
>>>>>>> 54d3161fecd5b6c9572fc5d53c16923574550627
#==============================================================================#

##===============================================================================
##								 Welcome message
##===============================================================================
print '\n\n'
print '|=======================================================|'
print '|\t\t\t\t\t\t\t|'
print '|\t\t :-) PCA MD (-:	\t\t\t|\n','|\t\t\t\t\t\t\t|\n', '|\t\twritten by Bilal Nizami\t\t\t|\n', '|\t\tRhodes University, 2017\t\t\t|'
print '|\t    distributed under GNU GPL 3.0\t\t|'
print '|\t\t\t\t\t\t\t|'
print '|=======================================================|'
print '\n'

#==============================================================================
#							   Setting the options
#==============================================================================

<<<<<<< HEAD
parser = OptionParser("Usage: pca.py -t <MD trajecotry> -p <topology file>  -a <atom group >")

parser.add_option("-t", "--trj", type='string', dest="trj",
                  help="file name of the MD trajectory")

parser.add_option("-p", "--top", type='string', dest="topology",
                  help="topology file")      

parser.add_option("-a", "--ag", type='string', dest="atm_grp",
                  help="group of atom for PCA. Options are:\
				  all= all atoms, backbone = backbone atoms, CA= C alpha atoms" )				  
=======
parser = OptionParser("Usage: pca.py -t <MD trajecotry> -p <topology file> ")

parser.add_option("-t", "--trj", type='string', dest="trj",
                  help="Name of the MD trajectory")

parser.add_option("-p", "--top", type='string', dest="topology",
                  help="topology file")                                            
>>>>>>> 54d3161fecd5b6c9572fc5d53c16923574550627

(options, args) = parser.parse_args()
                     

#====================================================================
# if no arguments are passed
#====================================================================
if options.trj is None: 
	print 'Missing trajectory arguments :(\nPlease see the help by running \n\nsystem_setup.py -h\n\n '
	parser.print_help()
	sys.exit(1)

if options.topology is None:
	print 'Missing topology !!\nPlease see the help by running \n\nsystem_setup.py -h\n\n '
	parser.print_help()
	sys.exit(1)

# assign the passed arguments
traj = options.trj
topology = options.topology


<<<<<<< HEAD
# read the trajecotry 
pca_traj = md.load(traj, top=topology)
#print(pca_traj)
top = pca_traj.topology
#print top
#==========================================================================
#
#		selecting the atoms 
#
# User passes the arguements to select the subset of atoms for PCA
#===========================================================================
if options.atm_grp is None:
	print 'No atom group selected. Please select a atom group for PCA analysis\n'
	parser.print_help()
	sys.exit(1)
	
atm_name = options.atm_grp
print "Group selected for PCA:", atm_name, "\n"

def select_atoms():

	if options.atm_grp is 'CA':
		print 'Selected CA for PCA\n'
		
		
	if options.atm_grp is 'BB':
		print 'Selected Backbone atoms for PCA\n'
		
	return;


sele_grp = top.select("name CA")

print len(sele_grp)
#sele_trj = pca_traj.xyz[:,sele_grp,:]
#print "Number of frames in the trajectory selected for PCA:", len(sele_trj), "\n"


#===============================================================
#
#  RMSD in reference with first frame
#
#===============================================================
rmsd = md.rmsd(pca_traj, pca_traj, 0, atom_indices=sele_grp)
print "RMSD written to rmsd.agr \n "
## write the RMSD file
np.savetxt('rmsd.agr', rmsd)
rf = open('rmsd.agr', 'r')
rf_cont = rf.read()
rf.close()

time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
title = '\tcreated by pca.py\t'
legends = '@    title "RMSD"\n\
@    xaxis  label "Time"\n\
@    yaxis  label "RMSD"\n\
@	TYPE xy\n'

pf = open('rmsd.agr', 'w')
pf.write('#'+title+'\ton\t'+time+'\n'+legends+'\n'+rf_cont)
pf.close()

#===============================================================
#
# PCA using sci-kit learn library
#===============================================================


def my_pca1():
	#pca_align = PCA() 		## initialize PCA object 
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	#pca_traj_reshaped = pca_traj.xyz.reshape(pca_traj.n_frames, pca_traj.n_atoms * 3) # reshape the cordinates of trajectory
	#pca_align.fit(pca_traj_reshaped)							# Fit the model
	#pca_reduced = pca_align.transform(pca_traj_reshaped)		 #apply the dimensionality reduction.
	
	# with selected trajecotry
	#sele_trj.superpose(pca_traj, 0, atom_indices=sele_grp)   # selected trajecotry
	
	
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	pca_sele_traj = PCA()
	pca_sele_traj.fit(sele_traj_reshaped)
	pca_sele_traj_reduced = pca_sele_traj.transform(sele_traj_reshaped)
	print "Trace of the covariance matrix is: ", np.trace(pca_sele_traj.get_covariance())
	print "Wrote covariance matrix..."
	np.savetxt('cov.dat', pca_sele_traj.get_covariance())
	np.savetxt('pca_projection.agr', pca_sele_traj_reduced)
	
		
	##print pca_align.get_covariance()
	#print(pca_reduced.shape)
	#print (type(pca_reduced))
	#np.savetxt('pc_test.agr', np.array(pca_reduced ))
	## write the eigenvectors to a file
	#np.savetxt('pc.agr', pca_reduced )
	
	pf = open('pca_projection.agr', 'r')
	pf_cont = pf.read()
	pf.close()
	
	time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends = '@    title "Projection of PC"\n\
	@    xaxis  label "PC1"\n\
	@    yaxis  label "PC2"\n\
	@	TYPE xy\n'
	
	pf = open('pca_projection.agr', 'w')
	pf.write('#'+title+'\ton\t'+time+'\n'+legends+'\n'+pf_cont)
	pf.close()
	
	## write PCs and explained_variance_ratio_
	expl_ratio = pca_sele_traj.explained_variance_ratio_
	expl_ratio = expl_ratio*100   # to make it percent
	
	print expl_ratio.reshape((1,101)).shape
	np.savetxt("explained_variance.agr", expl_ratio)
	
	ef = open('explained_variance.agr', 'r')
	ef_cont = ef.read()
	ef.close()
	
	legends = '@    title "explained_variance of PCs"\n\
	@    xaxis  label "PCs"\n\
	@    yaxis  label "% Variance"\n\
	@	TYPE xy\n'
	
	ef = open('explained_variance.agr', 'w')
	ef.write('#'+title+'\ton\t'+time+'\n'+legends+'\n'+ef_cont)
	ef.close()
	
	return;

my_pca1()

# ==============================================================
#
#
#  My own method
#
#=================================================================
def my_pca():
	print "shape of the trajectory object: ", sele_trj.shape, "\n"
	print pca_traj.xyz[:,sele_grp,:].reshape(pca_traj.n_frames, 860 * 3).T.shape

	pca_traj.xyz[0,sele_grp,:].reshape(1, 860 * 3).tofile('ca_cord_tr.xyz', sep = ' ',format = '%s')
	pca_traj.xyz[0, sele_grp].transpose().tofile('ca_cord.xyz', sep = '\t', format = '%s')
	x= pca_traj.xyz[:,sele_grp,:].reshape(pca_traj.n_frames, 860 * 3).T
	x = x.astype(float) ## to avoid numpy Conversion Error during scaling
	x_scaled = preprocessing.scale(x)
	#print x_scaled
	cov_mat = np.cov(x_scaled)
	print '\nwrote covariance matrix\n '
	np.savetxt('cov.txt', cov_mat )
	print 'covariance matrix shape is: ', cov_mat.shape, "\n"
	print "Trace of the covariance matrix is" , np.trace(cov_mat), "\n"
	return;

#my_pca()

### another methods
def my_pca2():
	sele_trj = pca_traj.xyz[:,sele_grp,:]
	pca2 = PCA()
	
	from itertools import combinations
	# this python function gives you all unique pairs of elements from a list
	
	atom_pairs = list(combinations(range(len(sele_trj)), 2))
	print list(combinations(range(4),2))
	
	pairwise_distances = md.geometry.compute_distances(pca_traj, atom_pairs)
	print(pairwise_distances.shape)
	reduced_distances = pca2.fit_transform(pairwise_distances)
	
	#np.savetxt('pc1.agr',reduced_distances)
	return;
=======
# read the trajecotry
pca_traj = md.load(traj, top=topology)
print(pca_traj)

## PCA using sci-kit learn library

#pca_align = PCA(n_components=2)

# Superpose each conformation in the trajectory upon first frame

#pca_traj.superpose(pca_traj, 0)



#reduced_cartesian = pca_align.fit_transform(pca_traj.xyz.reshape(pca_traj.n_frames, pca_traj.n_atoms * 3)) ## Fit the model with and apply the dimensionality reduction.

#print(reduced_cartesian.shape)
#print (np.array(reduced_cartesian))


## write the eigenvectors to a file
#np.savetxt('pc.agr', reduced_cartesian )
#pf = open('pc.agr', 'r')
#pf_cont = pf.read()
#pf.close()

time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
title = '\tcreated by pca.py\t'
legends = '@    title "Projection of PC"\n\
@    xaxis  label "PC1"\n\
@    yaxis  label "PC2"\n\
@	TYPE xy\n'

#pf = open('pc.agr', 'w')
#pf.write('#'+title+'\ton\t'+time+'\n'+legends+'\n'+pf_cont)


#print pf_cont
#pf.close()

### another methods

pca2 = PCA(n_components=2)

from itertools import combinations
# this python function gives you all unique pairs of elements from a list

atom_pairs = list(combinations(range(pca_traj.n_atoms), 2))
print list(combinations(range(4),2))

#pairwise_distances = md.geometry.compute_distances(pca_traj, atom_pairs)
#print(pairwise_distances.shape)
#reduced_distances = pca2.fit_transform(pairwise_distances)

#np.savetxt(pc1.agr,reduced_distances)
>>>>>>> 54d3161fecd5b6c9572fc5d53c16923574550627



