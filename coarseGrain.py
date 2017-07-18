#!/usr/bin/env python
# Takes a protomer structure and coarse grains to select a set amount of C-Beta
import os
import sys
import argparse
from datetime import datetime

from utils import *

import math


def main(args):
    f = open(args.pdbFile, 'r')
    lines = f.readlines()
    f.close()

    # Get Info to rewrite PDB
    header = []
    atomLines = []
    end = []

    # get index of first atom
    for i in range(len(lines)):
        if lines[i].startswith("ATOM"):
            index_atom = i
            break
    header = lines[0:index_atom]
    all_atoms = lines[index_atom:]
    number_of_protomers = 0
    macro_molecule = False
    chain_dics = {}
    number_protomer_atoms = 0
    c_beta_atoms = []
    change = False
    for atom in all_atoms:
        cbeta = []
        if atom.startswith("ATOM"):
            info = atom.split()
            chain = info[4]
            res = info[3]
            atom_type = info[2]
            if atom_type == "CB" or (atom_type == "CA" and res == "GLY"):
                if number_protomer_atoms == 0:
                    chain1 = chain
                if chain not in chain_dics:
                    number_protomer_atoms += 1
                if chain != chain1:
                    if chain1 not in chain_dics:
                        chain_dics[chain1] = 1
                    else:
                        chain_dics[chain1] += 1
                        macro_molecule = True
                    chain1 = chain
                c_beta_atoms.append(atom)
        if atom.startswith("END"):
            if chain1 in chain_dics:
                chain_dics[chain1] += 1
            else:
                chain_dics[chain1] = 1

    # print number_protomer_atoms
    # print chain_dics
    number_of_protomers = chain_dics[chain1]

    # Set level of coarsegrain
    c_g = args.cg  # parameter

    # Read in all cbeta atoms for
    cbetas = []
    for line in lines:
        if line.startswith("ATOM"):
            coords = []
            res = line.split()[3].strip()
            atom_type = line.split()[2].strip()
            if atom_type == "CB" or (atom_type == "CA" and res == "GLY"):
                coords.append(float(line[30:38].strip().strip()))
                coords.append(float(line[38:46].strip().strip()))
                coords.append(float(line[46:54].strip().strip()))
                cbetas.append(coords)

    index_of_selected_atoms = []

    starting_atom = args.startingAtom  # residue number of starting atoms
    starting_atom_i = starting_atom - 1  # Index for starting atom
    index_of_selected_atoms.append(starting_atom_i)

    number_protomer_atoms = args.protomerAtoms  # parameter
    protomer_c_betas = cbetas[0:number_protomer_atoms]
    coords_start = protomer_c_betas[starting_atom_i]
    distances_from_start = []
    distance_index = {}  # holds the index of the atoms in order of distance from

    xstart = coords_start[0]
    ystart = coords_start[1]
    zstart = coords_start[2]

    for i in range(len(protomer_c_betas)):
        atom = protomer_c_betas[i]
        if i == starting_atom_i:
            continue
        x = atom[0]
        y = atom[1]
        z = atom[2]
        distance = ((xstart - x) * (xstart - x)) + ((ystart - y)
                                                    * (ystart - y)) + ((zstart - z) * (zstart - z))
        distance = math.sqrt(distance)
        distances_from_start.append(distance)
        distance_index[distance] = i
    distances_from_start.sort()

    # selects atoms which are not within this distance to already selected atoms
    cutoff = distances_from_start[(c_g * (c_g - 1)) - c_g]
    atom_index = distance_index[cutoff]
    index_of_selected_atoms.append(atom_index)

    distribution = {}

    # loops through ordered list and selects all suitable atoms
    for dist in distances_from_start[((c_g * (c_g - 1)) - c_g) + 1:]:
        atom_index = distance_index[dist]
        x = protomer_c_betas[atom_index][0]
        y = protomer_c_betas[atom_index][1]
        z = protomer_c_betas[atom_index][2]
        too_close = False
        local_distribution = []
        for atom in index_of_selected_atoms:
            x1 = protomer_c_betas[atom][0]
            y1 = protomer_c_betas[atom][1]
            z1 = protomer_c_betas[atom][2]
            surrounding_distance = ((x1 - x) * (x1 - x)) + \
                                   ((y1 - y) * (y1 - y)) + ((z1 - z) * (z1 - z))
            surrounding_distance = math.sqrt(surrounding_distance)
            local_distribution.append(surrounding_distance)
            if surrounding_distance < cutoff:
                too_close = True
                break
        if not too_close:
            index_of_selected_atoms.append(atom_index)
            local_distribution.sort()
            distribution[atom_index] = local_distribution

    # print index_of_selected_atoms
    print "No. atoms selected per protomer: " + str(len(index_of_selected_atoms))
    print "No. atoms selected per macro molecule: " + str(len(index_of_selected_atoms) * number_of_protomers)
    # print "No. atoms selected per macro molecule: " + str(len(index_of_selected_atoms) * 60)

    # for a in distribution:
    # print a
    # print distribution[a]

    index_of_selected_atoms.sort()
    # write a pdb file of the Coarse-Grained Capsid with renumbered atoms
    selected_c_beta_lines = []
    # Includes all cbetas of the first pentamer and then coarse grains the rest of the surrounding capsid
    count = 0

    for i in range(0, number_of_protomers):  # parameter
        for j in index_of_selected_atoms:
            index = j + i * number_protomer_atoms
            a = c_beta_atoms[index]
            spaces = " " * (len(a[6:11]) - len(str(count + 1)))
            a = a[0:6] + spaces + str(count + 1) + a[11:]
            selected_c_beta_lines.append(a)
            count += 1
        ter = "TER" + " " * 3 + " " * (5 - len(str(count))) + str(count) + " " * 6 + a.split(
        )[3] + " " + a.split()[4] + " " + " " * (3 - len(a.split()[5])) + a.split()[5] + " \n"
        selected_c_beta_lines.append(ter)

    output_dir = args.outdir

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    w = open(output_dir + str(c_g) + "_SCA.pdb", 'w')
    w.writelines(header)
    w.writelines(selected_c_beta_lines)
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

    # standard arguments for logging and output
    parser.add_argument("--silent", help="Turn off logging",
                        action='store_true', default=False)
    parser.add_argument(
        "--log-file", help="Output log file (default: standard output)", default=None)
    parser.add_argument(
        "--outdir", help="Output directory", default="output")

    # custom arguments
    # parser.add_argument("--output", help="")#output = "3VBSPent"
    parser.add_argument("--pdbFile", help="PDB input file")  # 3VBSPent.pdb
    parser.add_argument(
        "--cg", help="Course grain level [int]", default=4, type=int)
    parser.add_argument(
        "--startingAtom", help="Residue number of starting atoms [int]", default=1, type=int)
    parser.add_argument("--protomerAtoms",
                        help="", default=842, type=int)

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
        print "No arguments provided. Use -h to view help"