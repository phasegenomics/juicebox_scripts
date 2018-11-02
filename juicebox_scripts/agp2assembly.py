#!/usr/bin/env python
# converts agp format to assembly format (useful for juicebox)
# requires exactly 2 inputs, writes to stdout

from __future__ import print_function

import sys
from collections import defaultdict

def printUsage():
    print("\nagp2assembly.py usage:", end='')
    print("\tagp2assembly.py <input_agp_file> <output_assembly_file>")
    return

def read_from_agp(filename):
    lines = []
    clusters = defaultdict(list)
    counter = 0
    order = []
    with open(filename) as agp:
        for line in agp:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[4] != "W":
                continue
            counter += 1
            lines.append(">{0} {1} {2}\n".format(fields[5], counter, fields[7]))
            if fields[8] == "-":
                this_contig = -1 * counter
            else:
                this_contig = counter
            
            if fields[0] not in order:
            	order.append(fields[0])
            	
            clusters[fields[0]].append(str(this_contig))
    return lines, clusters, order

def write_assembly(lines, clusters, order, outfilename):
    with open(outfilename, "w") as outfile:
        for line in lines:
            outfile.write(line)
        for cluster in order:
            outfile.write(" ".join(clusters[cluster]) + "\n")
    return 0

def main():
    if len(sys.argv) != 3:
        printUsage()
        sys.exit()
    lines, clusters, order = read_from_agp(sys.argv[1])
    write_assembly(lines, clusters, order, sys.argv[2])

if __name__ == "__main__":
    main()
