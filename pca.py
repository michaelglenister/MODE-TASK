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
from sklearn import preprocessing


#==============================================================================#
#											PCA MD 
#
#			This programe performs the PCA on a MD trajectory
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#=============================================================================================+=======
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

<<<<<<< HEAD
parser = OptionParser("Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >")
=======

parser = OptionParser("Usage: pca.py -t <MD trajecotry> -p <topology file>  -a <atom group >")
>>>>>>> 64b6917d827745068dd25b6ff1ef5996f3637ad0

parser.add_option("-t", "--trj", type='string', dest="trj",
                  help="file name of the MD trajectory")

parser.add_option("-p", "--top", type='string', dest="topology",
                  help="topology file")      

parser.add_option("-a", "--ag", type='string', dest="atm_grp",
<<<<<<< HEAD
                  help="group of atom for PCA. Default is C alpha atoms. Other options are :"
				  "all= all atoms, backbone = backbone atoms, CA= C alpha atoms, protein= protein's atoms")	

parser.add_option("-r", "--ref", type='string', dest="reference",
                  help="reference structure for RMSD") 

parser.add_option("-m", "--pca_type", type='string', dest="pca_type",
                  help="PCA method. Default is normal PCA. Options are:\
				  KernelPCA, normal, ipca. If normal additional arguments can be passed by flag -svd. If KernelPCA is selected kernel type can also be defined by flag -k") 
				  
parser.add_option("-k", "--kernel_type", type='string', dest="kernel_type",
                  help="Type of kernel for KernalPCA. default is linear. Options are :"
				  "linear, poly, rbf, sigmoid, cosine, precomputed") 

parser.add_option("-s", "--svd_solver", type='string', dest="svd_solver",
                  help="Type of svd_solver  for PCA. Default is auto. Options are :"
				  "auto, full, arpack, randomized") 
=======
                  help="group of atom for PCA. Options are:\

				  all= all atoms, backbone = backbone atoms, CA= C alpha atoms\
				  protein= protein's atoms")				  
#=======
				  all= all atoms, backbone = backbone atoms, CA= C alpha atoms" )				  
#=======
parser = OptionParser("Usage: pca.py -t <MD trajecotry> -p <topology file> ")

parser.add_option("-t", "--trj", type='string', dest="trj",
                  help="Name of the MD trajectory")

parser.add_option("-p", "--top", type='string', dest="topology",
                  help="topology file")                                            

>>>>>>> 64b6917d827745068dd25b6ff1ef5996f3637ad0

(options, args) = parser.parse_args()
                     
atm_name = options.atm_grp
print "Atom group selected for PCA:", atm_name, "\n"


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

#=======================================
# assign the passed arguments and read the trajectory 
#=======================================

traj = options.trj
topology = options.topology
ref = options.reference
ptype=options.pca_type
if ref:
	ref = md.load(options.reference)

<<<<<<< HEAD
=======


# read the trajecotry 
>>>>>>> 64b6917d827745068dd25b6ff1ef5996f3637ad0
pca_traj = md.load(traj, top=topology)
#print(pca_traj)
top = pca_traj.topology
#print top

#==============================================
#
# Setting the default options
#
#===============================================

if options.atm_grp == None:
	print 'No atom has been selected. PCA will be performed on C alpha atoms '
	atm_name = 'CA'  # set to default C-alpha atoms
if options.reference == None:
	print "No reference structure given, RMSD will be computed to the first frame in the trajectory"
	ref = pca_traj # set reference to current trajectory
if options.pca_type == None:
	ptype = 'normal'

if options.svd_solver == None:
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
	
	if options.reference == None: 
		print "Reference for RMSD calculation is first frame of trajectory"
	else:
		print "Reference for RMSD calculation is: ", options.reference
		
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

## write plots
def write_plots(file_name, pca):
	'function to write pca plots. takes name of the file to write and pca object name'
	fname = ''
	fname = file_name+'.agr'
	np.savetxt(fname, pca)
	pf = open(fname, 'r')
	pf_cont = pf.read()
	pf.close()
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends = '@    title "Projection of PC"\n\
	@    xaxis  label "PC1"\n\
	@    yaxis  label "PC2"\n\
	@	TYPE xy\n'
	
	pf = open(fname, 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+pf_cont)
	pf.close()
	
	return;

def write_pcs(file_name, pca):
	'write PCs and explained_variance_ratio_. takes name of the file to write and pca object name'
	fname = ''
	fname = file_name+'.agr'
	#print type(pca)
	e_ratio = pca.explained_variance_ratio_
	e_ratio = e_ratio*100   # to make it percent
	
	print e_ratio.reshape((1,101)).shape
	np.savetxt(fname, e_ratio)
	
	ef = open(fname, 'r')
	ef_cont = ef.read()
	ef.close()
	
	title = '\tcreated by pca.py\t'
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	legends = '@    title "explained_variance of PCs"\n\
	@    xaxis  label "PCs"\n\
	@    yaxis  label "% Variance"\n\
	@	TYPE xy\n'
	
	ef = open(fname, 'w')
	ef.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+ef_cont)
	ef.close()
	return;

#===============================================================
#
# PCA using sci-kit learn library
#===============================================================

def my_pca(svd):
	 
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	pca_sele_traj = PCA()
	pca_sele_traj.fit(sele_traj_reshaped)
	pca_sele_traj_reduced = pca_sele_traj.transform(sele_traj_reshaped)
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


if ptype == 'KernelPCA':
	kernel = ''
	kernel = options.kernel_type
	if options.kernel_type:
		print "Performing Kernel PCA with", kernel, 'kernel'
		my_kernelPCA(kernel)
	else:
		print "Performing Kernel PCA with default linear kernel"
		my_kernelPCA('linear')

if ptype == 'normal':
	svd=''
	svd = options.svd_solver
	if svd:
		print "Performing normal PCA with ",svd,"svd_solver"
		my_pca(svd)
	else:
		print "Performing normal PCA with 'auto' svd_solver"
		my_pca(svd)
if ptype == 'ipca':
	print "Performing Incremental_pca (IPCA)"
	incremental_pca()


#===============================================================
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




