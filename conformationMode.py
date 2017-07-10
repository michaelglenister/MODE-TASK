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
    linesEmpty = f.readlines()
    f.close()

    f = open(args.pdbProtomerAligned, 'r')
    linesFull = f.readlines()
    f.close()

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
    print FullResidues
    print EmptyResidues
    print CommonResidues

    SelectedFull = {'A': [], 'B': [], 'C': [], 'D': []}
    EmptyCords = []
    FullCords = []
    for line in linesEmpty:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            resType = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if res in CommonResidues[chain]:
                if (atype == "CB" or (atype == "CA" and resType == "GLY")):
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    cod = [x, y, z]
                    EmptyCords.append(cod)

    for line in linesFull:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            resType = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if res in CommonResidues[chain]:
                if res not in SelectedFull[chain]:
                    SelectedFull[chain].append(res)
                if (atype == "CB" or (atype == "CA" and resType == "GLY")):
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    cod = [x, y, z]
                    FullCords.append(cod)

    # Calculate deltaR

    deltaR = []
    NumberOfProtomers = 1

    for i in range(NumberOfProtomers):
        for j in range(len(EmptyCords) / NumberOfProtomers):
            Full = FullCords[i + (j * NumberOfProtomers)]
            Empt = EmptyCords[i + (j * NumberOfProtomers)]
            print Full
            print Empt
            rx = Empt[0] - Full[0]
            ry = Empt[1] - Full[1]
            rz = Empt[2] - Full[2]
            deltaR.append(rx)
            deltaR.append(ry)
            deltaR.append(rz)

    print len(deltaR)

    # Calculate the magnitude
    magDR = 0
    for r in deltaR:
        magDR += r * r
    magDR = sqrt(magDR)

    # Get Eigenvetors for a mode
    # Calculate Indexes of vectors to be selected
    interfaceIndex = []
    f = open(args.pdb_sca, 'r')
    NMA = f.readlines()
    f.close()
    count = 0
    for line in NMA:
        if line.startswith("ATOM"):
            info = line.split()
            atype = info[2].strip()
            resType = info[3].strip()
            chain = info[4].strip()
            res = int(info[5].strip())
            if (atype == "CB" or (atype == "CA" and resType == "GLY")):
                if res in CommonResidues[chain]:
                    interfaceIndex.append(count)
                count += 1

    # Mode range (Starting with 10 slowest modes)
    #CG = 621
    #Protomer = 2520
    ModeRange = range(0, 621)
    f = open(args.vtProtomer, 'r')
    vectors = f.readlines()
    f.close()

    output = []
    overlayList = []
    for mode in ModeRange:
        Overlap = 0
        CommonVector = []
        vector = vectors[mode].split()
        for res in interfaceIndex:
            for i in range(3):
                ele = float(vector[res * 3 + i])
                CommonVector.append(ele)

        # Calculate the magnitude
        magMode = 0
        for r in CommonVector:
            magMode += r * r
        magMode = sqrt(magMode)
        # Calculate Dot Product
        if len(CommonVector) == len(deltaR):
            print "Vectors Match"
            for i in range(len(CommonVector)):
                Overlap += CommonVector[i] * deltaR[i]

            Overlap = Overlap / (magDR * magMode)

            output.append("Mode: " + str(mode) +
                          " Overlap: " + str(Overlap) + '\n')
            overlayList.append(Overlap)

    overlayList.sort()
    overlayList.reverse()

    w = open(args.output, 'w')
    for out in output:
        w.write(out)
    w.write("Sorted Values:\n")
    for o in overlayList:
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

    # custom arguments
    parser.add_argument("--pdbAligned", help="")  # '4n43_aligned.pdb'
	# missing pdb_sca args?
    parser.add_argument("--pdbProtomerAligned", help="") # '3VBSProtomer_aligned3_SCA.pdb'
    parser.add_argument("--pdbProtomer", help="")  # '3VBSProtomer3_SCA.pdb' used where?
    parser.add_argument("--vtProtomer", help="")  # 'Protomer3CG_VT.txt'
    parser.add_argument("--output", help="Output file",
                        default="ProtomerCGrained.txt")  # 'Protomer3CG_VT.txt'

    args = parser.parse_args()

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
