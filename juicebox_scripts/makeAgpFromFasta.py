#!/usr/bin/env python
#requires exactly 2 inputs
#a fasta file to print contig lengths
#writes results to stdout

from __future__ import print_function
import sys

def printUsage():
    print("\nmakeAgpFromFasta.py usage:", end='')
    print("\tmakeAgpFromFasta.py <fasta_file> <agp_out_file>")
    return

#this funtion returns true if line is a contig+bin line
def isContigBinLine(line):
    ret = 0
    if len(line) > 1:
        ret = (line[0:1] == ">")
    return ret

#this function returns the contig+bin portion of a contig+bin line
def getContigBinFromLine(line):
    ret = ""
    if isContigBinLine(line):
        #trim the >, split the line, and grab the first token, in case there is a long description on the line
        ret = line[1:].split()[0].strip()
    return ret


def main():
    if len(sys.argv) != 3:
        printUsage()
        sys.exit()

    fname = sys.argv[1]
    outfile = sys.argv[2]
    contig = ""

    counting = 0
    count = 0
    total = 0

    with open(fname) as file:
        with open(outfile, 'w') as outf:
            outf.write("##agp-version 2.1\n")
            for line in file:
                line = line.strip()
                if isContigBinLine(line):
                    if counting:
                        outf.write("{0}\t0\t{1}\t1\tW\t{0}\t1\t{1}\t+\n".format(contig, count))
                    else:
                        counting = 1
                    contig = getContigBinFromLine(line)
                    count = 0
                elif counting:
                    count += len(line)
                    total += len(line)
                else:
                    continue
    
            outf.write("{0}\t0\t{1}\t1\tW\t{0}\t1\t{1}\t+\n".format(contig, count))
            
if __name__ == "__main__":
    main()
