#!/usr/bin/env python
import os
import sys
import argparse
from datetime import datetime
from utils import *
import numpy as np
from math import sqrt

def main(args):
    # Read In PDB files
    f = open(args.pdbConfAligned, 'r')
    lines_empty = f.readlines()
    f.close()

    f = open(args.pdbProtAligned, 'r')
    lines_full = f.readlines()
    f.close()

    empty_residues = {}
    full_residues = {}

    common_residues = {}

    #determine the number of assymetric units
    number_of_protomers = 0
    currentResidue = 0

    for line in lines_full:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            res = int(info[5].strip())
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
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
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
		if chain in empty_residues:
                    if res not in empty_residues[chain]:
                        empty_residues[chain].append(res)
		else:
			empty_residues[chain] = [res]

    for line in lines_full:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):	
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
                    if atype == "CB" or (atype == "CA" and res_type == "GLY"):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        cod = [x, y, z]
                        empty_cords.append(cod)

    for line in lines_full:
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
                    if atype == "CB" or (atype == "CA" and res_type == "GLY"):
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        cod = [x, y, z]
                        full_cords.append(cod)

    # Calculate deltaR

    delta_r = []
    

    for i in range(number_of_protomers):
        for j in range(len(empty_cords) / number_of_protomers):
            full = full_cords[i + (j * number_of_protomers)]
            empty = empty_cords[i + (j * number_of_protomers)]
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
    f = open(args.pdbANM, 'r')
    nma = f.readlines()
    f.close()
    count = 0
    for line in nma:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
                if chain in common_residues:
                    if res in common_residues[chain]:
                        interface_index.append(count)
                count += 1

 
    f = open(args.vtMatrix, 'r')
    vectors = f.readlines()
    f.close()
    mode_range = range(len(vectors)-6) #excludes the trivial modes
	
    
    output = {}
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
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    parser.add_argument("--outdir", help="Output directory", default="output")

    # custom arguments
    parser.add_argument("--pdbConfAligned", help="")

    parser.add_argument("--pdbProtAligned", help="")
    parser.add_argument("--pdbANM", help="")
    parser.add_argument("--vtMatrix", help="")  # note: change this from vtProtomer
    parser.add_argument("--output", help="Output file", default="ModesOfConformtionalChange.txt")

    args = parser.parse_args()

    print "!=====================================================================================!"
    print "! Pleas check the following:                                                          !"
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
