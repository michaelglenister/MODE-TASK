#!/usr/bin/env python
import os
import sys
import argparse
from datetime import datetime

from utils import *

from math import sqrt
# Identifies Modes responsible for conformational change for a molecule wth 15 copies of each atom


def main(args):
    # Read In PDB files
    f = open(args.pdbAligned, 'r')
    lines_empty = f.readlines()
    f.close()

    f = open(args.pdbProtomerAligned, 'r')
    lines_full = f.readlines()
    f.close()

    empty_residues = {'A': [], 'B': [], 'C': [], 'D': []}
    full_residues = {'A': [], 'B': [], 'C': [], 'D': []}

    common_residues = {'A': [], 'B': [], 'C': [], 'D': []}

    for line in lines_empty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if atype == "CB" or (atype == "CA" and res_type == "GLY"):
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
                if res not in full_residues[chain]:
                    full_residues[chain].append(res)
    for ch in full_residues:
        for r in full_residues[ch]:
            if r in empty_residues[ch]:
                common_residues[ch].append(r)
    # print full_residues
    # print empty_residues
    # print common_residues

    selected_full = {'A': [], 'B': [], 'C': [], 'D': []}
    empty_cords = []
    full_cords = []
    for line in lines_empty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            res_type = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
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
    number_of_protomers = 1

    for i in range(number_of_protomers):
        for j in range(len(empty_cords) / number_of_protomers):
            full = full_cords[i + (j * number_of_protomers)]
            empty = empty_cords[i + (j * number_of_protomers)]
            # print full
            # print empty
            rx = empty[0] - full[0]
            ry = empty[1] - full[1]
            rz = empty[2] - full[2]
            delta_r.append(rx)
            delta_r.append(ry)
            delta_r.append(rz)

    # print len(delta_r)

    # Calculate the magnitude
    mag_d_r = 0
    for r in delta_r:
        mag_d_r += r * r
    mag_d_r = sqrt(mag_d_r)

    # Get Eigenvectors for a mode
    # Calculate Indexes of vectors to be selected
    interface_index = []
    f = open(args.pdbSca, 'r')
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
                if res in common_residues[chain]:
                    interface_index.append(count)
                count += 1

    # Mode range (Starting with 10 slowest modes)
    # CG = 621
    # Protomer = 2520
    mode_range = range(0, 621)
    f = open(args.vtProtomer, 'r')
    vectors = f.readlines()
    f.close()

    output = []
    overlay_list = []
    for mode in mode_range:
        overlap = 0
        common_vector = []
        vector = vectors[mode].split()
        for res in interface_index:
            for i in range(3):
                ele = float(vector[res * 3 + i])
                common_vector.append(ele)

        # Calculate the magnitude
        mag_mode = 0
        for r in common_vector:
            mag_mode += r * r
        mag_mode = sqrt(mag_mode)
        # Calculate Dot Product
        if len(common_vector) == len(delta_r):
            # print "Vectors Match"
            for i in range(len(common_vector)):
                overlap += common_vector[i] * delta_r[i]

            overlap = overlap / (mag_d_r * mag_mode)

            output.append("Mode: " + str(mode) +
                          " Overlap: " + str(overlap) + '\n')
            overlay_list.append(overlap)

    overlay_list.sort()
    overlay_list.reverse()

    w = open(args.outdir + "/" + args.output, 'w')
    for out in output:
        w.write(out)
    w.write("Sorted Values:\n")
    for o in overlay_list:
        w.write(str(o) + "\n")
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
    parser.add_argument("--pdbAligned", help="")  # '4n43_aligned.pdb'
    # missing pdbSca args?
    parser.add_argument("--pdbProtomerAligned", help="") # '3VBSProtomer_aligned3_SCA.pdb'
    parser.add_argument("--pdbSca", help="")  # '3VBSProtomer3_SCA.pdb'
    parser.add_argument("--vtProtomer", help="")  # 'Protomer3CG_VT.txt'
    parser.add_argument("--output", help="Output file",
                        default="ProtomerCGrained.txt")  # 'Protomer3CG_VT.txt'

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
