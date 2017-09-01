#!/usr/bin/env python

# coformationalMode.py
# Identifies normal modes that act in the direction of a conformational change
# Author: Caroline Ross: caroross299@gmail.com
# August 2017

import os
import sys
import argparse
from datetime import datetime
from lib.utils import *
import numpy as np
from math import sqrt
import lib.sdrms
import numpy as np

def main(args):

    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
	print '\n**************************************\nUnrecognised atom type\nInput Options:\nCa: to select alpha carbon atoms\CB: to select beta carbon atoms\n**************************************'
	sys.exit()



    if args.pdbConf==args.pdbANM:
	print '\n**************************************\nWARNING!!!\nConformational change PDB files are the same:\n--pdbANM: '+args.pdbANM+'\n--pdbConf: '+args.pdbConf+'\n**************************************\n'
	

    try:
        f = open(args.pdbConf, 'r')
        lines_empty = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.pdbConf+' NOT FOUND:\n**************************************\n'
	sys.exit()



    try:
        f = open(args.pdbANM, 'r')
        nma = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.pdbANM+' NOT FOUND:\n**************************************\n'
	sys.exit()

    empty_residues = {}
    full_residues = {}

    common_residues = {}

    #determine the number of assymetric units
    number_of_protomers = 0
    currentResidue = 0

    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
		if currentResidue == 0:
			currentResidue = res
			number_of_protomers+=1	
		elif res==currentResidue:
			number_of_protomers+=1
		else:
			break
   

    #determine co-ords of residues common in both PDB files		

    for line in lines_empty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
		if chain in empty_residues:
                    if res not in empty_residues[chain]:
                        empty_residues[chain].append(res)
		else:
			empty_residues[chain] = [res]

    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
	        if chain in full_residues:
                    if res not in full_residues[chain]:
                        full_residues[chain].append(res)
		else:
			 full_residues[chain]=[res]


    for ch in full_residues:
        for r in full_residues[ch]:
	    if ch in empty_residues:
                if r in empty_residues[ch]:
		    if ch in common_residues:
                    	common_residues[ch].append(r)
		    else:
			common_residues[ch] = [r]
	    else:
		break

    
    
   
    
    selected_full = {}
    for ch in common_residues:
	selected_full[ch] = []

    count_common = 0
    empty_cords = []
    full_cords = []
    for line in lines_empty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
	    if chain in common_residues:
		
                if res in common_residues[chain]:
                    if atype == atomT or (atype == "CA" and res_type == "GLY"):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        cod = [x, y, z]
                        empty_cords.append(cod)
			count_common+=1

    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
	    if chain in common_residues:		
                if res in common_residues[chain]:
                    if res not in selected_full[chain]:
                        selected_full[chain].append(res)
                    if atype == atomT or (atype == "CA" and res_type == "GLY"):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        cod = [x, y, z]
                        full_cords.append(cod)


    target = np.zeros((count_common,3))
    for i,res in enumerate(empty_cords):
	for j,c in enumerate(res):
		target[i,j] = c
	

    struc = np.zeros((count_common,3))
    for i,res in enumerate(full_cords):
	for j,c in enumerate(res):
		struc[i,j] = c
    confAligned = sdrms.superpose3D(struc, target)
    alignedFull = confAligned[0]
    rmsd = confAligned[1]
    
    print '\nRMSD between '+args.pdbANM+' and '+args.pdbConf+' = '+str(rmsd)+'\n'
    # Calculate deltaR
    delta_r = []
    for j in range(len(alignedFull)):
      full = alignedFull[j]
      empty = target[j]
      rx = empty[0] - full[0]
      ry = empty[1] - full[1]
      rz = empty[2] - full[2]
      delta_r.append(rx)
      delta_r.append(ry)
      delta_r.append(rz)

    # Calculate the magnitude
    correlationDR = []
    csum=0
    mag_d_r = 0
    for i,r in enumerate(delta_r):
        mag_d_r += r * r
	csum = csum+r*r
	if i%3==2:
            csum = sqrt(csum)
	    correlationDR.append(csum)
	    csum=0	
    mag_d_r = sqrt(mag_d_r)

    # Get Eigenvectors for a mode
    # Calculate Indexes of vectors to be selected
    interface_index = []
    count = 0
    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
                if chain in common_residues:
                    if res in common_residues[chain]:
                        interface_index.append(count)
                count += 1


    try:   
        f = open(args.vtMatrix, 'r')
        vectors = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.vtMatrix+' NOT FOUND:\n**************************************\n'
	sys.exit()
    mode_range = range(len(vectors)-6) #excludes the trivial modes
	
    
    output = {}
    try:
        for mode in mode_range:
	    correlationMode = []
            overlap = 0
            common_vector = []
            vector = vectors[mode].split()
            for res in interface_index:
                csum = 0
                for i in range(3):
                    ele = float(vector[res * 3 + i])
		    csum = csum+ele*ele
                    common_vector.append(ele)
	        csum = sqrt(csum)
	        correlationMode.append(csum)


            # Calculate the magnitude
            mag_mode = 0
            for r in common_vector:
                mag_mode += r * r
	
            mag_mode = sqrt(mag_mode)
            # Calculate Dot Product
            if len(common_vector) == len(delta_r):
                # print "Vectors Match"
                C = abs(np.corrcoef(correlationMode,correlationDR)[0,1])
                for i in range(len(common_vector)):
                    overlap += common_vector[i] * delta_r[i]

                overlap = overlap / (mag_d_r * mag_mode)
	    
                spaces = len("mode: "+str(mode+1))
	    
	        spaces = 15-spaces
	    
	        if abs(overlap) in output:
                    output[abs(overlap)].append("Mode: " + str(mode+1) + ' '*spaces+ str(overlap) +'      '+str(C)+'\n')
	        else:
		    output[abs(overlap)]=["Mode: " + str(mode+1) + ' '*spaces + str(overlap)+'      '+ str(C)+'\n']

        overlap_list = output.keys()
        overlap_list.sort()
        overlap_list.reverse()

        w = open(args.outdir + "/" + args.output, 'w')
        w.write('MODE           Overlap              Correlation\n\n')
        for out in overlap_list:
	    for o in output[out]:
                w.write(o)

        w.close()
    except IndexError:
        print '\n**************************************\nFILE '+args.vtMatrix+' IS NOT A VALID EIGENVECTOR FILE:\n**************************************\n'
        sys.exit()
silent = False
stream = sys.stdout


def log(message):
    global silent
    global stream

    if not silent:
        print >> stream, message


if __name__ == "__main__":
    # parse cmd arguments
    parser = argparse.ArgumentParser()

    # standard arguments for logging
    parser.add_argument("--silent", help="Turn off logging", action='store_true', default=False)
    parser.add_argument("--welcome", help="Display welcome message (true/false)", default="true")
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    parser.add_argument("--outdir", help="Output directory", default="output")

    # custom arguments
    parser.add_argument("--pdbConf", help="")
    parser.add_argument("--pdbANM", help="")
    parser.add_argument("--vtMatrix", help="")  # note: change this from vtProtomer
    parser.add_argument("--output", help="Output file", default="ModesOfConformtionalChange.txt")
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='CA')

    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Conformation mode", "Caroline Ross (caroross299@gmail.com)")

    print "!=====================================================================================!"
    print "! Please check the following:                                                          !"
    print "! --pdbANM must be the PDB file that NMA was performed on                             !"
    print "! --pdbProtAligned must be a PDB of your complex aligned to the conformational change !"
    print "!=====================================================================================!"

    # Check if required directories exist
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Check if args supplied by user
    if len(sys.argv) > 1:
        # set up logging
        silent = args.silent

        if args.log_file:
            stream = open(args.log_file, 'w')

        start = datetime.now()
        log("Started at: %s" % str(start))

        # run script
        main(args)

        end = datetime.now()
        time_taken = format_seconds((end - start).seconds)

        log("Completed at: %s" % str(end))
        log("- Total time: %s" % str(time_taken))

        # close logging stream
        stream.close()
    else:
        print "No arguments provided. Use -h to view help"
