#!/usr/bin/env python
import os
import sys
import argparse
from datetime import datetime

from utils import *


def main(args):

    # Takes two pdb models and determines the common residues
    #####################################
    f = open(args.conf1, 'r')
    lines_c1 = f.readlines()
    f.close()

    f = open(args.conf2, 'r')
    lines_c2 = f.readlines()
    f.close()
    #####################################

    c1_residues = {}
    c2_residues = {}

    common_residues = {'A': [], 'B': [], 'C': [], 'D': []}

    for line in lines_c1:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
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
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
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
		break

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
    parser.add_argument("--log-file", help="Output log file (default: standard output)", default=None)
    parser.add_argument("--outdir", help="Output directory", default="output")

    # custom arguments
    parser.add_argument("--conf1", help="commonResidues.py extracts the common residues between two PDB files eg, two protein conformations, two PDB of different Coarse grained resolution\n--conf1 takes the path to the first of two pdb files")  
    parser.add_argument("--conf2", help="commonResidues.py extracts the common residues between two PDB files eg, two protein conformations, two PDB of different Coarse grained resolution\n--conf2 takes the path to the second of two pdb files")  

    args = parser.parse_args()

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
