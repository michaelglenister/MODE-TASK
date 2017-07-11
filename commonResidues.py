#!/usr/bin/env python
import os
import sys
import argparse
from datetime import datetime

from utils import *


def main(args):
    # Takes two pdb models and determines the common residues
    #####################################
    f = open(args.fullCapsid, 'r')
    linesEmpty = f.readlines()
    f.close()

    f = open(args.protomer, 'r')
    linesFull = f.readlines()
    f.close()
    #####################################

    EmptyResidues = {'A': [], 'B': [], 'C': [], 'D': []}
    FullResidues = {'A': [], 'B': [], 'C': [], 'D': []}

    CommonResidues = {'A': [], 'B': [], 'C': [], 'D': []}

    for line in linesEmpty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            resType = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if (atype == "CB" or (atype == "CA" and resType == "GLY")):
                if res not in EmptyResidues[chain]:
                    EmptyResidues[chain].append(res)
    for line in linesFull:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            resType = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if (atype == "CB" or (atype == "CA" and resType == "GLY")):
                if res not in FullResidues[chain]:
                    FullResidues[chain].append(res)

    for ch in FullResidues:
        for r in FullResidues[ch]:
            if r in EmptyResidues[ch]:
                CommonResidues[ch].append(r)

    print CommonResidues


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
    parser.add_argument("--silent", help="Turn off logging",
                        action='store_true', default=False)
    parser.add_argument(
        "--log-file", help="Output log file (default: standard output)", default=None)

    # custom arguments
    parser.add_argument("--fullCapsid", help="")  # '3VBSFull4_SCA.pdb'
    parser.add_argument("--protomer", help="")  # '3VBSProtomer3_SCA.pdb'

    args = parser.parse_args()

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
		print "No argumeants provided. Use -h to view help"
