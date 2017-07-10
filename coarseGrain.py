#!/usr/bin/env python
# Takes a protomer structure and coarse grains to select a set amount of C-Beta
import os
import sys
import argparse
from datetime import datetime

from utils import *

import math


def main(args):
    output = args.pdbFile[:args.pdbFile.index('.')]  # paramter
    f = open(args.pdbFile, 'r')  # paramter
    lines = f.readlines()
    f.close()

    # Get Info to rewrite PDB
    header = []
    atomLines = []
    end = []

    # get index of first atom
    for i in range(len(lines)):
        if lines[i].startswith("ATOM"):
            indexAtom = i
            break
    header = lines[0:indexAtom]
    AllAtoms = lines[indexAtom:]

    CBetaAtoms = []
    for atom in AllAtoms:
        CBeta = []
        if atom.startswith("ATOM"):
            info = atom.split()
            first = info[0].strip()
            res = info[3]
            atomType = info[2]
            if first == "ATOM" and (atomType == "CB" or (atomType == "CA" and res == "GLY")):
                CBetaAtoms.append(atom)

    # Set level of coarsegrain
    CG = args.cg  # paramater

    # Read in all cbeta atoms for
    cbetas = []
    for line in lines:
        if line.startswith("ATOM"):
            coords = []
            res = line.split()[3].strip()
            atomType = line.split()[2].strip()
            if atomType == "CB" or (atomType == "CA" and res == "GLY"):

                coords.append(float(line[30:38].strip().strip()))
                coords.append(float(line[38:46].strip().strip()))
                coords.append(float(line[46:54].strip().strip()))
                cbetas.append(coords)

    IndexOfSelectedAtoms = []

    startingAtom = args.startingAtom  # residue number of starting atoms
    startingAtomI = startingAtom - 1  # Index for starting atom
    IndexOfSelectedAtoms.append(startingAtomI)

    NumberProtomerAtoms = args.protomerAtoms  # paratmer
    ProtomerCBetas = cbetas[0:NumberProtomerAtoms]
    coordsStart = ProtomerCBetas[startingAtomI]
    distancesFromStart = []
    distanceIndex = {}  # holds the index of the atoms in order of distance from

    xstart = coordsStart[0]
    ystart = coordsStart[1]
    zstart = coordsStart[2]

    for i in range(len(ProtomerCBetas)):
        atom = ProtomerCBetas[i]
        if i == startingAtomI:
            continue
        x = atom[0]
        y = atom[1]
        z = atom[2]
        distance = ((xstart - x) * (xstart - x)) + ((ystart - y)
                                                    * (ystart - y)) + ((zstart - z) * (zstart - z))
        distance = math.sqrt(distance)
        distancesFromStart.append(distance)
        distanceIndex[distance] = i
    distancesFromStart.sort()

    # selects atoms which are not within this distance to already selected atoms
    cutoff = distancesFromStart[(CG * (CG - 1)) - CG]
    atomIndex = distanceIndex[cutoff]
    IndexOfSelectedAtoms.append(atomIndex)

    distribution = {}

    # loops through ordered list and selects all suitable atoms
    for dist in distancesFromStart[((CG * (CG - 1)) - CG) + 1:]:
        atomIndex = distanceIndex[dist]
        x = ProtomerCBetas[atomIndex][0]
        y = ProtomerCBetas[atomIndex][1]
        z = ProtomerCBetas[atomIndex][2]
        tooClose = False
        localDistribution = []
        for atom in IndexOfSelectedAtoms:
            x1 = ProtomerCBetas[atom][0]
            y1 = ProtomerCBetas[atom][1]
            z1 = ProtomerCBetas[atom][2]
            surroundingDistance = ((x1 - x) * (x1 - x)) + \
                ((y1 - y) * (y1 - y)) + ((z1 - z) * (z1 - z))
            surroundingDistance = math.sqrt(surroundingDistance)
            localDistribution.append(surroundingDistance)
            if surroundingDistance < cutoff:
                tooClose = True
                break
        if not tooClose:
            IndexOfSelectedAtoms.append(atomIndex)
            localDistribution.sort()
            distribution[atomIndex] = localDistribution

    # print IndexOfSelectedAtoms
    print "No. atoms selected per protomer: " + str(len(IndexOfSelectedAtoms))
    print "No. atoms selected per Capsid: " + str(len(IndexOfSelectedAtoms) * 60)

    # for a in distribution:
    #	print a
    #	print distribution[a]

    IndexOfSelectedAtoms.sort()
    # write a pdb file of the Coarse-Grained Capsid with renumbered atoms
    SelectedCBetaLines = []
    # Includes all cbetas of the first penatmer and then coarse grains the rest of the surrounding capsid
    count = 0
    '''for i in range(5):
		for j in range(NumberProtomerAtoms):
		index = j+i*NumberProtomerAtoms
		a = CBetaAtoms[index]
		spaces = " "*(len(a[6:11])-len(str(count+1)))
		a = a[0:6]+spaces+str(count+1)+a[11:]
		SelectedCBetaLines.append(a)
		count+=1
		ter = "TER"+" "*3+" "*(5-len(str(count)))+str(count)+" "*6+a.split()[3]+" "+a.split()[4]+" "+" "*(3-len(a.split()[5]))+a.split()[5]+" \n"
		SelectedCBetaLines.append(ter)'''

    for i in range(0, 5):  # paramarter
        for j in IndexOfSelectedAtoms:
            index = j + i * NumberProtomerAtoms
            a = CBetaAtoms[index]
            spaces = " " * (len(a[6:11]) - len(str(count + 1)))
            a = a[0:6] + spaces + str(count + 1) + a[11:]
            SelectedCBetaLines.append(a)
            count += 1
        ter = "TER" + " " * 3 + " " * (5 - len(str(count))) + str(count) + " " * 6 + a.split(
        )[3] + " " + a.split()[4] + " " + " " * (3 - len(a.split()[5])) + a.split()[5] + " \n"
        SelectedCBetaLines.append(ter)
    w = open(output + str(CG) + "_SCA.pdb", 'w')
    w.writelines(header)
    w.writelines(SelectedCBetaLines)
    w.write("END")
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
    # parser.add_argument("--output", help="")#output = "3VBSPent"
    parser.add_argument("--pdbFile", help="PDB input file")  # 3VBSPent.pdb
    parser.add_argument("--cg", help="Course grain level [int]", default=4, type=int)
    parser.add_argument(
        "--startingAtom", help="Residue number of starting atoms [int]", default=15, type=int)
    parser.add_argument("--protomerAtoms",
                        help="", default=842, type=int)

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
