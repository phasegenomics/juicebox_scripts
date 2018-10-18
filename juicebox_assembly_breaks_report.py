# some suggestions from code review 4/25/18:
# TODO: add docstrings overall and for each function
# TODO: alter __str__ method from just outputting agp, maybe just use a non-__str__ method.
# TODO: alter __repr__ methods to something more intuitive
# TODO: unit tests!


import argparse
import re
import copy
import os.path

from collections import defaultdict
from natsort import natsorted

class seq_range:
    def __init__ (self, name, index, end, feature_type, last_end=0, last_name=''):
        self.name            = name
        self.index           = int(index)
        self.start           = 0  # was 0 + 1. altered below for writing agp.
        self.end             = int(end)  # was int(end) - 2. altered below for writing agp.
        self.fragment        = False
        self.strand          = '+'
        self.type            = feature_type
        self.scaffold_offset = 0
        self.sid             = 0
        self.part_number     = 0
        self.is_zero_length          = False  # whether or not to use the seq...

        names = self.name.split(':::')
        self.name = names[0]

        if last_name.split(':::')[0] != self.name:
            last_end = 0

        if(len(names) > 1):
            if 'frag' in names[1]:
                self.fragment = True
                #self.name += ':::{0}'.format(names[1])
            if 'debris' in names[1:] and self.start == self.end:
                self.is_zero_length = True  # zero-length contigs are bad.
                #if debris is interesting, do something

        if self.fragment and last_end > 0:
            self.start += last_end
            self.end += last_end

    def setselfid(self, id):
        self.sid = id
    def setpartnumber(self, pn):
        self.part_number = pn

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--assembly', help='juicebox assembly file', required=True)
parser.add_argument('-p', '--prefix', help='the prefix to use for writing outputs', required=True)
args = parser.parse_args()

printed_contigs = defaultdict(int)
#the last bin always has an off-by-one error in its length because that's how Juicebox's .assembly file
#format works, so we want to track the last bins so we can correct it
last_bins_for_contigs = defaultdict(int)
scaffolds   = []
scaffold_id = 1
new_seq_end = 0
new_seq_name = ''
output_file = '{0}.broken_contigs.txt'.format(args.prefix)
seq_items = dict()
out_lines = dict()

with open(args.assembly) as f:
    sline = [x.strip().split() for x in f.readlines()]
    for items in sline:
        if items[0][0] != '>':
            for i in items:
                i = str(abs(int(i)))
                seq = seq_items[str(abs(int(i)))]
                if not seq.fragment:
                    continue
                if seq.index == last_bins_for_contigs[seq.name]:
                    #correct the off-by-one error in the last bin of the contig
                    seq.end = seq.end - 1
                out_lines[i] = '{0}\t{0}_bin{1}\t{2}\t{3}\t{4}'.format(seq.name, printed_contigs[seq.name], seq.start, seq.end, (seq.end-seq.start))
                printed_contigs[seq.name] += 1
        else:
            new_seq = seq_range(items[0][1:], items[1], items[2], 'W', new_seq_end, new_seq_name)
            new_seq_end = new_seq.end
            new_seq_name = new_seq.name
            seq_items[items[1]] = new_seq
            last_bins_for_contigs[new_seq_name] = int(items[1])

#count up the total number of breaks that were introduced
#For every 2 bins beyond the first for a contig, 1 break was introduced.
#It's also impossible to have an even number of bins if using Juicebox
#to make the breaks because it always excises a small sequence at the break
#point. E.g:
# 1 bin    => no breaks
# 2 bins   => impossible
# 3 bins   => 1 break
# 4 bins   => impossible
# 5 bins   => 2 breaks
# 5 bins   => impossible
# 7 bins   => 3 breaks
# etc.
breaks_per_contig = dict()
contigs_with_breaks = set()
for contig in printed_contigs:
    bins = printed_contigs[contig]
    if bins % 2 == 0:
        raise ValueError('Even number of bins detected for {0}: {1} bins'.format(contig, bins))
    breaks_per_contig[contig] = (bins - 1) / 2 #we could rely on integer division instead of subtracting 1, but c'mon, don't get cute
    contigs_with_breaks.add(contig)

total_breaks = sum(breaks_per_contig.values())
total_contigs_broken = len(contigs_with_breaks)

with open(output_file, 'w') as f:
    f.write('#break summary for {0}\n'.format(os.path.basename(args.assembly)))
    f.write('#{0} total breaks in {1} contigs\n'.format(total_breaks, total_contigs_broken))
    f.write('#orig_contig\tfragment\tbreak_start\tbreak_end\tfragment_len\n')
    for line in natsorted(out_lines):
        f.write('{0}\n'.format(out_lines[line]))


