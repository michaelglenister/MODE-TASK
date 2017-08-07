#!/usr/bin/python
#filename: pca.py
import numpy as np
from time import sleep, gmtime, strftime

#==============================================================================#
#											
#			This programe is the part of PCA MD. It writes the PCA plots in xmgrace formatted .agr file 
#
# 								Author : Bilal Nizami
# 						  	 Rhodes University, 2017
#==============================================================================#

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
	@	TYPE xy\n\
	@    s0 line type 0\n\
	@    s0 line linestyle 1\n\
	@    s0 line linewidth 1.0\n\
	@    s0 line color 1\n\
	@    s0 line pattern 1\n\
	@    s0 baseline type 0\n\
	@    s0 baseline off\n\
	@    s0 dropline off\n\
	@    s0 symbol 1\n\
	@    s0 symbol size 0.250000\n\
	@    s0 symbol color 1\n\
	@    s0 symbol pattern 1\n\
	@    s0 symbol fill color 1\n\
	@    s0 symbol fill pattern 1\n\
	@    s0 symbol linewidth 1.0\n\
	@    s0 symbol linestyle 1\n\
	@    s0 symbol char 25\n\
	@    s0 symbol char font 0\n\
	@    s0 symbol skip 0\n'
	
	pf = open(fname, 'w')
	pf.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+pf_cont)
	pf.close()
	
	return;
	
## write PCs 
def write_pcs(file_name, pca):
	'write PCs and explained_variance_ratio_. takes name of the file to write and pca object name'
	fname = ''
	fname = file_name+'.agr'
	#print type(pca)
	e_ratio = pca.explained_variance_ratio_
	e_ratio = e_ratio*100   # to make it percent
	
	np.savetxt(fname, e_ratio)
	
	ef = open(fname, 'r')
	ef_cont = ef.read()
	ef.close()
	title = '\tcreated by pca.py\t'
	my_time = strftime("%Y-%m-%d  %a  %H:%M:%S", gmtime())
	legends = '@    title "explained_variance of PCs"\n\
	@    xaxis  label "PCs"\n\
	@    yaxis  label "% Variance"\n\
	@	TYPE xy\n\
	@    s0 symbol 1\n\
	@    s0 symbol size 0.250000\n\
	@    s0 symbol color 1\n\
	@    s0 symbol pattern 1\n\
	@    s0 symbol fill color 1\n\
	@    s0 symbol fill pattern 1\n\
	@    s0 symbol linewidth 1.0\n\
	@    s0 symbol linestyle 1\n\
	@    s0 symbol char 25\n\
	@	s0 symbol fill color 2\n\
	@	s0 symbol color 2\n\
	@    s0 symbol char font 0\n\
	@    s0 symbol skip 0\n'
	
	ef = open(fname, 'w')
	ef.write('#'+title+'\ton\t'+my_time+'\n'+legends+'\n'+ef_cont)
	ef.close()
	return;
	
