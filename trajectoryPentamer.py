#!/usr/bin/env python
# Makes a trajectory of 100 PDB files
import os
import sys
import argparse
from datetime import datetime

from decimal import Decimal


def main(args):
    pdbf = open(args.pdb, 'r')  # CHANGE HERE
    pdblines = pdbf.readlines()
    pdbf.close()

    # CHANGE THE FOLLOWING
    modeF = args.modeF
    modeL = args.modeL
    modes = []
    for i in range(modeF, modeL + 1):
        modes.append(str(i))

    for mode in modes:
        structure = "Protomer"
        # Get the header and co-ords

        header = []
        atomLines = []
        coOrds = []
        end = []

        # get index of first atom
        for i in range(len(pdblines)):
            if pdblines[i].startswith("ATOM"):
                indexAtom = i
            break
        header = pdblines[0:indexAtom]
        AllAtoms = pdblines[indexAtom:]

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
            else:
                if "TER" in atom or "END" in atom:
                    CBetaAtoms.append(atom)

        # Renumber the atoms
        for i in range(len(CBetaAtoms) - 1):
            a = CBetaAtoms[i]
            spaces = " " * (len(a[6:11]) - len(str(i + 1)))
            a = a[0:6] + spaces + str(i + 1) + a[11:]
            CBetaAtoms[i] = a

        # Determine Connections
        CONECT = []
        ConectChain = []
        atom1 = CBetaAtoms[0]
        chain1 = atom1.split()[4].strip()

        for i, atom in enumerate(CBetaAtoms, 1):
            if "TER" in atom or "END" in atom:
                continue
            chain = atom.split()[4].strip()
            if chain == chain1:
                ConectChain.append(i)
            else:
                CONECT.append(ConectChain)
                ConectChain = []
                ConectChain.append(i)
                chain1 = chain
        CONECT.append(ConectChain)

        # Get the vectors
        # CHANGE HERE # may not be correct change, double check
        vectorf = open(args.modeFilePrefix + mode + ".txt", 'r')
        vectors = vectorf.readlines()
        vectorf.close()

        # Write Trajectories
        w = open("Trajectories/" + structure + '_' + mode + ".pdb", 'w')
        for i in range(0, 100):
            w.writelines(header)
            vIndex = -1
            for atom in CBetaAtoms:
                if "ATOM" in atom:
                    vIndex += 1
                    v = vectors[vIndex].split()
                    vx = float(v[0].strip())
                    vy = float(v[1].strip())
                    vz = float(v[2].strip())
                    x = round(float(atom[30:38].strip()) + (vx * i / 5), 3)
                    xspace = ' ' * (len(atom[30:38]) - len(str(x)))
                    y = round(float(atom[38:46].strip()) + (vy * i / 5), 3)
                    yspace = ' ' * (len(atom[38:46]) - len(str(y)))
                    z = round(float(atom[46:54].strip()) + (vz * i / 5), 3)
                    zspace = ' ' * (len(atom[46:54]) - len(str(z)))
                    atom = atom[0:30] + xspace + \
                        str(x) + yspace + str(y) + zspace + str(z) + atom[54:]
                    w.write(atom)
                else:
                    if "TER " in atom:
                        w.write(atom)
                    else:
                        if "END" in atom:
                            for Con in CONECT:
                                for c in range(len(Con) - 1):
                                    atom1 = str(Con[c])
                                    atom2 = str(Con[c + 1])
                                    w.write("CONECT" + " " * (5 - len(atom1)) +
                                            atom1 + " " * (5 - len(atom2)) + atom2 + "\n")
                            w.write(atom)

        w.close()

    # write arrows

    arrows = ['proc vmd_draw_arrow {mol start end} {\n    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]\n    graphics $mol cylinder $start $middle radius 0.20\n    graphics $mol cone $middle $end radius 0.35\n}\ndraw color green\n']
    steps = [0, 5]
    vIndex = -1
    chainbreaks = [CONECT[0][-1], CONECT[1][-1], CONECT[2][-1], CONECT[3][-1]]
    for atom in CBetaAtoms:
        if "ATOM" in atom:

            vIndex += 1
            v = vectors[vIndex].split()
            vx = float(v[0].strip())
            vy = float(v[1].strip())
            vz = float(v[2].strip())
            x1 = round(float(atom[30:38].strip()) + (vx * steps[0]), 3)
            y1 = round(float(atom[38:46].strip()) + (vy * steps[0]), 3)
            z1 = round(float(atom[46:54].strip()) + (vz * steps[0]), 3)

            x2 = round(float(atom[30:38].strip()) + (vx * steps[1]), 3)
            y2 = round(float(atom[38:46].strip()) + (vy * steps[1]), 3)
            z2 = round(float(atom[46:54].strip()) + (vz * steps[1]), 3)
            arrows.append('draw arrow {' + str(x1) + ' ' + str(y1) + ' ' + str(
                z1) + '} {' + str(x2) + ' ' + str(y2) + ' ' + str(z2) + '}\n')
            if vIndex + 1 == chainbreaks[0]:
                arrows.append('draw color blue2\n')
            elif vIndex + 1 == chainbreaks[1]:
                arrows.append('draw color orange\n')
            elif vIndex + 1 == chainbreaks[2]:
                arrows.append('draw color yellow\n')

    w = open("Trajectories/" + structure + 'ARROWS' + mode + ".txt", 'w')
    w.writelines(arrows)
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
    # '3VBSFull_Aligned.pdb'
    parser.add_argument("--pdb", help="")  # '3VBSProtomer3_SCA.pdb'
    parser.add_argument("--modeF", help="[int]", default=617, type=int)
    parser.add_argument("--modeL", help="[int]", default=617, type=int)
    parser.add_argument("--modeFilePrefix", help="")  # 'ProtomerMode'

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
