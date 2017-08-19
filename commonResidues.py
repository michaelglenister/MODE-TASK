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
    lines_empty = f.readlines()
    f.close()

    f = open(args.protomer, 'r')
    lines_full = f.readlines()
    f.close()
    #####################################

    empty_residues = {}
    full_residues = {}

    common_residues = {'A': [], 'B': [], 'C': [], 'D': []}

    for line in lines_empty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
                if chain not in empty_residues:
                    empty_residues[chain] = [res]
                else:
                    if res not in empty_residues[chain]:
                        empty_residues[chain].append(res)
    for line in lines_full:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
                if chain not in full_residues:
                    full_residues[chain] = [res]
                else:
                    if res not in full_residues[chain]:
                        full_residues[chain].append(res)

    for ch in full_residues:
        for r in full_residues[ch]:
            if r in empty_residues[ch]:
                common_residues[ch].append(r)

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
    parser.add_argument("--fullCapsid", help="")  # '3VBSFull4_SCA.pdb'
    parser.add_argument("--protomer", help="")  # '3VBSProtomer3_SCA.pdb'

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
