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

parser = OptionParser("Usage: pca.py -t <MD trajectory> -p <topology file>  -a <atom group >")

parser.add_option("-t", "--trj", type='string', dest="trj",
                  help="file name of the MD trajectory")

parser.add_option("-p", "--top", type='string', dest="topology",
                  help="topology file")      

parser.add_option("-a", "--ag", type='string', dest="atm_grp",
                  help="group of atom for PCA. Default is C alpha atoms. Other options are:\
				  all= all atoms\
				  backbone = backbone atoms\
				  CA= C alpha atoms\
				  protein= protein's atoms")	

parser.add_option("-r", "--ref", type='string', dest="reference",
                  help="reference structure for RMSD") 

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
if ref:
	ref = md.load(options.reference)

pca_traj = md.load(traj, top=topology)
#print(pca_traj)
top = pca_traj.topology
#print top

#==============================================
#
# Setting the default options
#
#===============================================

if options.atm_grp is None:
	print 'No atom has been selected. PCA will be performed on C alpha atoms '
	atm_name = 'CA'  # set to default C-alpha atoms
if options.reference == None:
	print "No reference structure given, RMSD will be computed to the first frame in the trajectory"
	ref = pca_traj # set reference to current trajectory


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
		
	if atm_name == None:
		print "Total", len(sele_grp), 'C-alpha atoms selected for analysis\n'
	else:
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


def my_pca1():
	 
	pca_traj.superpose(pca_traj, 0, atom_indices=sele_grp) 			# Superpose each conformation in the trajectory upon first frame
	sele_trj = pca_traj.xyz[:,sele_grp,:]												# select cordinates of selected atom groups
	sele_traj_reshaped = sele_trj.reshape(pca_traj.n_frames, len(sele_grp) * 3)
	pca_sele_traj = PCA()
	pca_sele_traj.fit(sele_traj_reshaped)
	pca_sele_traj_reduced = pca_sele_traj.transform(sele_traj_reshaped)
	print "Trace of the covariance matrix is: ", np.trace(pca_sele_traj.get_covariance())
	print "Wrote covariance matrix..."
	
	print "processing...\\",
	syms = ['\\', '|', '/', '-']
	bs = '\b'

	# print spinner ..need to fix
	for _ in range(10):
		for sym in syms:
			sys.stdout.write("\b%s" % sym)
			sys.stdout.flush()
			time.sleep(.5)
	np.savetxt('cov.dat', pca_sele_traj.get_covariance())
	np.savetxt('pca_projection.agr', pca_sele_traj_reduced)
	
	
	pf = open('pca_projection.agr', 'r')
	pf_cont = pf.read()
	pf.close()
	
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	title = '\tcreated by pca.py\t'
	legends = '@    title "Projection of PC"\n\
	@    xaxis  label "PC1"\n\
	@    yaxis  label "PC2"\n\
	@	TYPE xy\n'
	
	pf = open('pca_projection.agr', 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+pf_cont)
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
	ef.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+ef_cont)
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



