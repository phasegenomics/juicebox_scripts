#/usr/bin/env python

from __future__ import print_function
import sys

with open(sys.argv[1]) as file:
    gap_index = "not an index"
    for line in file:
        if line.startswith(">"):
            if line.startswith(">hic_gap"):
                gap_index = line.split()[1]
                print("omitting entries with index {} as gaps".format(gap_index), file=sys.stderr)
                continue
            print(line.strip())
        else:
            fields = line.split()
            degapped = [field for field in fields if field != gap_index]
            print(" ".join(degapped))
