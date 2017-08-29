#!/usr/bin/env python


# commonResidues.py
# Determines the common residues between two sets of PDB files 
# Eg. a protein with two levels of coarse graining 
# Eg. a protein with two conformations
# Author: Caroline Ross: caroross299@gmail.com
# August 2017

import os
import sys
import argparse
from datetime import datetime

from utils import *


def main(args):

    atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
	print '\n**************************************\nUnrecognised atom type\nInput Options:\nCa: to select alpha carbon atoms\CB: to select beta carbon atoms\n**************************************'
	sys.exit()


    #getCommonResidues:
    # Takes two pdb models and determines the common residues
    #####################################
    try:	
        f = open(args.conf1, 'r')
        lines_c1 = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.conf1+' NOT FOUND:\n**************************************\n'
	sys.exit()
    try:
        f = open(args.conf2, 'r')
        lines_c2 = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.conf2+' NOT FOUND:\n**************************************\n'
	sys.exit()
    #####################################

    c1_residues = {}
    c2_residues = {}

    common_residues = {}

    for line in lines_c1:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
                
                if chain in c1_residues:
                    if res not in c1_residues[chain]:
                        c1_residues[chain].append(res)
                else:
	            c1_residues[chain] = [res]


    for line in lines_c2:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res_type == "GLY"):
                if chain in c2_residues:
                    if res not in c2_residues[chain]:
                        c2_residues[chain].append(res)
                else:
	            c2_residues[chain] = [res]
	    else:
                print "NO " +str(res)
                print atype
                print atomT
                print res_type


    for ch in c2_residues:
        for r in c2_residues[ch]:
	    if ch in c1_residues:
                if r in c1_residues[ch]:
		    if ch in common_residues:
                    	common_residues[ch].append(r)
		    else:
			common_residues[ch] = [r]
	    else:
		break
    '''atomT = args.atomType.upper()
    if atomT!='CA' and atomT!='CB':
	print '\n**************************************\nUnrecognised atom type\nInput Options:\nCA: to select alpha carbon atoms\nCB: to select beta carbon atoms\n**************************************\n'
	sys.exit()
 
    #####################################
    try:
        f = open(args.conf1, 'r')
        lines_c1 = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.conf1+' NOT FOUND:\n**************************************\n'
	sys.exit()
 
    try:
        f = open(args.conf2, 'r')
        lines_c2 = f.readlines()
        f.close()
    except IOError:
        print '\n**************************************\nFILE '+args.conf2+' NOT FOUND:\n**************************************\n'
	sys.exit()
 
    #####################################

    c1_residues = {}
    c2_residues = {}

    common_residues = {}

    for line in lines_c1:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res == "GLY"):
                if chain not in c1_residues:
                    c1_residues[chain] = [res]
                else:
                    if res not in c1_residues[chain]:
                        c1_residues[chain].append(res)
    for line in lines_c2:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == atomT or (atype == "CA" and res == "GLY"):
                if chain not in c2_residues:
                    c2_residues[chain] = [res]
                else:
                    if res not in c2_residues[chain]:
                        c2_residues[chain].append(res)

    for ch in c2_residues:
        for r in c2_residues[ch]:
	    if ch in c1_residues:
                if r in c1_residues[ch]:
                    if ch in common_residues:
                    	common_residues[ch].append(r)
		    else:
			common_residues[ch] = [r]
	    else:
		break'''

    w = open(args.outdir + "/common_residues", 'w')
    w.write(str(common_residues))
    w.close()

    # print CommonResidues


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
    parser.add_argument("--atomType", help="Enter CA to select alpha carbons or CB to select beta carbons", default='CA')

    # custom arguments
    parser.add_argument("--conf1", help="commonResidues.py extracts the common residues between two PDB files eg, two protein conformations, two PDB of different Coarse grained resolution\n--conf1 takes the path to the first of two pdb files")  
    parser.add_argument("--conf2", help="commonResidues.py extracts the common residues between two PDB files eg, two protein conformations, two PDB of different Coarse grained resolution\n--conf2 takes the path to the second of two pdb files")  

    args = parser.parse_args()

    if args.welcome == "true":
        welcome_msg("Common residues", "Caroline Ross (caroross299@gmail.com)")

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
