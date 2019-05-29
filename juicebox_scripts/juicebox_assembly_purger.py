#!usr/bin/env python
'''
Brad Nelson
May 29, 2019
Phase Genomics

juicebox_scripts/juicebox_assembly_purger.py

This script takes an assembly file and outputs an assembly file
without the specified contigs.

Copyright 2019, Phase Genomics Inc. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see https://www.gnu.org/licenses/agpl-3.0.en.html
'''

from __future__ import print_function
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_assembly")
    parser.add_argument("output_assembly")
    parser.add_argument("--exclude_contigs",
                        default=None,
                        nargs="+",
                        help="Names of contigs to exclude")
    parser.add_argument("--exclude_file",
                        default=None,
                        help="Path to file of contigs to exclude (with contig names in "
                             "first whitespace-delimited column, one per line)")
    args = parser.parse_args()

    if args.exclude_contigs is None and args.exclude_file is None:
        raise ValueError("Error: exclude_contigs or exclude_file must be specified")

    return vars(args)

def get_exclude(contigs, file):
    """Takes a list and/or file of contig names and returns a set of contig names to exclude.
    Assumes contigs and file are not both None.
    """
    exclude = set([])
    if contigs is not None:
        for name in contigs:
            exclude.add(name)
    if file is not None:
        with open(file, 'r') as input:
            for line in input:
                name = line.split()[0]
                exclude.add(name)

    return exclude

def filter_assembly(exclude, input_assembly, output_assembly, **kwargs):
    with open(input_assembly, "r") as input, \
         open(output_assembly, "w") as output:
        index = 0
        purged_names = set([])
        purged_indices = set([])
        purged_from_scaffolds = set([])
        for line in input:
            if line.startswith(">"):
                name, num, length = line.strip().split()
                name = name.replace(">", "")
                if name in exclude:
                    purged_names.add(name)
                    purged_indices.add(num)
                    continue
                else:
                    outstring = ">{} {} {}\n".format(name, index, length)
                    output.write(outstring)
                    index += 1
            else:
                outlist = []
                contigs = line.strip().split()
                for contig_num in contigs:
                    orientation = "-" if contig_num.startswith("-") else "+"
                    if contig_num.replace("-", "") in purged_indices:
                        purged_from_scaffolds.add(contig_num.replace("-", ""))
                        continue
                    else:
                        outlist.append(contig_num)
                if len(outlist) > 0:
                    outstring = " ".join(outlist) + "\n"
                    output.write(outstring)
    all_exclude_contigs_found = exclude.issubset(purged_names) and purged_names.issubset(exclude)
    header_scaffold_match = len(purged_indices) == len(purged_from_scaffolds)

    if not all_exclude_contigs_found:
        not_found = exclude.difference(purged_names)
        mistaken_purge = purged_names.difference(exclude)
        if len(not_found) > 0:
            raise ValueError("Error: contigs specified for exclusion \
                             were not found:\n{}".format(" ".join(not_found)))
        elif len(mistaken_purge) > 0:
            raise ValueError("Error: contigs not specified for exclusion \
                             were mistakenly purged:\n{}".format(" ".join(mistaken_purge)))
    elif not header_scaffold_match:
        missing_from_scaffold = purged_from_scaffolds.difference(purged_indices)
        missing_from_header = purged_indices.difference(purged_from_scaffolds)
        if len(missing_from_scaffold) > 0:
            raise ValueError("Error: contigs found in header were missing from scaffolds: {}".format(" ".join(missing_from_scaffold)))
        elif len(missing_from_header) > 0:
            # This shouldn't happen, but it's included it for completeness
            raise ValueError("Error: contigs found in scaffolds were missing from header: {}".format(" ".join(missing_from_header)))

    print("SUCCESS")
    print("Purged {} contigs from header and scaffolds.".format(len(purged_names), len(purged_from_scaffolds)))
    print("Finished writing output assembly to {}.".format(output_assembly))

def main():
    args = parse_args()
    exclude = get_exclude(args["exclude_contigs"], args["exclude_file"])
    filter_assembly(exclude, **args)

if __name__ == "__main__":
    main()
