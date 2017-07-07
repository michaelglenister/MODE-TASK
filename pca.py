#!/usr/bin/python
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
print '|\t\t :-) PCA MD (-:	\t\t\t|\n','|\t\t\t\t\t\t\t|\n', '|\t\twritten by Bilal Nizami\t\t\t|\n', '|\t\tRhodes University, 2017\t\t\t|'
print '|\t    distributed under GNU GPL 3.0\t\t|'
print '|\t\t\t\t\t\t\t|'
print '|=======================================================|'
print '\n'

#==============================================================================
#							   Setting the options
#==============================================================================

parser = OptionParser("Usage: pca.py -t <MD trajecotry> -p <topology file> ")

parser.add_option("-t", "--trj", type='string', dest="trj",
                  help="Name of the MD trajectory")

parser.add_option("-p", "--top", type='string', dest="topology",
                  help="topology file")                                            

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


# read the trajecotry
pca_traj = md.load(traj, top=topology)
print(pca_traj)

## PCA using sci-kit learn library

pca_align = PCA(n_components=2)

# Superpose each conformation in the trajectory upon first frame

pca_traj.superpose(pca_traj, 0)

reduced_cartesian = pca_align.fit_transform(pca_traj.xyz.reshape(pca_traj.n_frames, pca_traj.n_atoms * 3))
print(reduced_cartesian.shape)
#print (reduced_cartesian)

## Project

plt.figure(figsize=(8,8))
plt.plot(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='x', c=pca_traj.time)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Cartesian coordinate PCA: alanine dipeptide')
cbar = plt.colorbar()
cbar.set_label('Time [ps]')
plt.show()


