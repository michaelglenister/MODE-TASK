#!/usr/bin/env python
# Calculates and Returns Diagonals of Correlated Matrix for a given set of modes
import os
import sys
import argparse
from datetime import datetime

from utils import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import numpy as np


# Lets say that the user has performed NMA on two coarse grained models of the same protein and now wants to compare
# to see if the additional coarse graining decreased the accuracy. If we obtain the same mean square fluctuations for
# each residue then in each model then we can say that the results are comparable regardless of the coarse graining
# level. But obviously must compare only the residues that are common in each model. hence we specify commonResidues
# (There is a script that works these out as well so we must just link the output of that) here i have specified by
# chain

def main(args):
    # residues common to multiple pdb files
    with open(args.commonResidues, 'r') as inf:
        common_residues = eval(inf.read())
    # print common_residues

    interface_index = []

    # The pdb file on which NMA analysis was performed, this script handles on model at a time but we can change this
    f = open(args.pdbProtomer, 'r')
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
            last_res = int(info[5].strip())
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
                if res in common_residues[chain]:
                    # gets the index of the common atoms, as they would appear in the output W, U and VT matrix. As
                    # these matrices contain info on all atoms
                    interface_index.append(count)
                count += 1

    protein_name = args.pdbProtomer
    protein_name = protein_name[protein_name.rfind("/")+1:protein_name.rfind("/")+5] + "_Protomer"

    # Specify modes
    total_modes = 0 #args.totalModes  # 2526 #number of residues in protein *3
    f = open(args.wMatrix, 'r')
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
        total_modes += 1

    first_mode = args.firstMode  # 2519
    last_mode = args.lastMode  # 2519

    # If user fails to provide first and last modes get default values
    if args.firstMode == 0 and args.lastMode == 0:
        first_mode = total_modes - 6
        last_mode = total_modes - 1

    # Llama

    first_res = 0
    # Specify Residue Indexes
    # Get first residue number
    for line in nma:  # Future work to include this in existing pdbProtomer file read above
        if line.startswith("ATOM"):
            if first_res == 0:
                info = line.split()
                first_res = int(info[1].strip())
    res_range = range(first_res - 1, last_res)
    print "Residue range: " + str(res_range)

    print "Residue count: " + str(last_res - first_res + 1)

    mode_range = range(first_mode, last_mode + 1)
    print "Mode range: " + str(mode_range)

    # '3VBSProtomerW.txt' # W matrix input file that was output from C++ Scripts
    fw = open(args.wMatrix, 'r')
    eigen_values = fw.readlines()
    fw.close()

    # Create A Full W Inverse Matrix (This is if we want correlation averaged over all modes)
    w_inv = np.zeros((total_modes, total_modes))
    for i in range(total_modes):
        if i < total_modes - 6:
            w_inv[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))
    # print "W in "

    # Create Filtered W Inverse Matrix (This is if we want correlation for a specific mode)
    w_f = np.zeros((total_modes, total_modes))
    for i in mode_range:
        w_f[i, i] = 1 / (float(eigen_values[i].split()[1].strip()))
    # print "WF in "

    # Read In U and VT full Matrix as U is the transpose of VT I only read in VT and create U from the VT matrix
    # info. So we can exclude U output from C++ script for faster analysis
    fvt = open(args.vtMatrix, 'r')  # '3VBSProtomerVT.txt'
    eigen_vectors = fvt.readlines()
    fvt.close()
    # print "U and VT file in "

    v_t = np.zeros((total_modes, total_modes))
    u = np.zeros((total_modes, total_modes))

    for i in range(total_modes):
        vectors = eigen_vectors[i].split()
        for j in range(total_modes):
            vector = float(vectors[j].strip())
            v_t[i, j] = vector
            u[j, i] = vector

    # print "U and VT read"
    # Calculate Correlation Matrices
    # Full C matrix
    w_v_t = np.dot(w_inv, v_t)
    # print "Correlations Calculated"
    c = np.dot(u, w_v_t)
    # print "Correlations Calculated"

    # Mode Specific C Matrix
    w_v_tm = np.dot(w_f, v_t)
    # print "Correlations Calculated"
    CM = np.dot(u, w_v_tm)

    # print "Correlations Calculated"

    # Calculate Trace of the Correlation Matrices
    trace_c = np.zeros((total_modes / 3, total_modes / 3))
    trace_c_m = np.zeros((total_modes / 3, total_modes / 3))

    for i in range(0, total_modes, 3):
        for j in range(0, total_modes, 3):
            trace = 0
            trace_m = 0
            for k in range(3):
                trace = trace + c[i + k, j + k]
                trace_m = trace_m + CM[i + k, j + k]
            trace_c[i / 3, j / 3] = trace
            trace_c_m[i / 3, j / 3] = trace_m

    # Print the diagonal values per residue
    w = open(args.outdir + "/" + protein_name + str(first_mode) + "BetaValues.txt", 'w')
    w.write("Full Correlation\n")
    w.write("Res\tBetaValue\n")
    for i in res_range:
        w.write(str(i + 1) + "\t" + str(trace_c[i, i]) + "\n")
    w.write("Common Residues\n")
    for i in interface_index:
        w.write(str(i + 1) + "\t" + str(trace_c[i, i]) + "\n")

    w.write("\nFiltered Correlation\n")
    w.write("Res\tBetaValue\n")
    for i in res_range:
        w.write(str(i + 1) + "\t" + str(trace_c_m[i, i]) + "\n")
    w.write("Common Residues\n")
    for i in interface_index:
        w.write(str(i + 1) + "\t" + str(trace_c_m[i, i]) + "\n")
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
    parser.add_argument("--silent", help="Turn off logging",
                        action='store_true', default=False)
    parser.add_argument(
        "--log-file", help="Output log file (default: standard output)", default=None)
    parser.add_argument(
        "--outdir", help="Output directory", default="output")

    # custom arguments
    parser.add_argument("--commonResidues", help="Files containing a dictionary like data set of common residues")  # '3VBSProtomer.pdb'
    parser.add_argument("--pdbProtomer", help="Input")  # '3VBSProtomer.pdb'
    #parser.add_argument("--totalModes", help="[int]", default=810, type=int)  # Generalise
    parser.add_argument("--firstMode", help="[int]", default=0, type=int)  # last mode - 7
    parser.add_argument("--lastMode", help="[int]", default=0, type=int)  # use last mode if no user input provided
    #parser.add_argument("--firstResidue", help="[int]", default=1, type=int)
    #parser.add_argument("--lastResidue", help="[int]", default=270, type=int)
    parser.add_argument(
        "--wMatrix", help="W matrix input file that was output from C++ Scripts")
    parser.add_argument("--vtMatrix", help="U and VT full Matrix")

    args = parser.parse_args()

    # Check if args supplied by user
    if len(sys.argv) > 1:
        # Check modes
        if args.firstMode > args.lastMode:
            print "First mode cannot be greater than last mode"
            sys.exit()


        # Check if required directories exist
        if not os.path.isdir(args.outdir):
            os.makedirs(args.outdir)

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
