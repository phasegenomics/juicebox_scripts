# some suggestions from code review 4/25/18:
# TODO: add docstrings overall and for each function
# TODO: alter __str__ method from just outputting agp, maybe just use a non-__str__ method.
# TODO: alter __repr__ methods to something more intuitive
# TODO: unit tests!


import argparse
import re
import copy

from pyfaidx import Fasta
from collections import defaultdict


def reverse_compliment(sequence):
    ret = ''
    for b in sequence[::-1]:
        if b == 'N':
            ret += 'N'
        elif b == 'A':
            ret += 'T'
        elif b == 'T':
            ret += 'A'
        elif b == 'G':
            ret += 'C'
        elif b == 'C':
            ret += 'G'
        elif b == 'a':
            ret += 't'
        elif b == 't':
            ret += 'a'
        elif b == 'g':
            ret += 'c'
        elif b == 'c':
            ret += 'g'
        elif b == 'n':
            ret += 'n'
        elif b == 'S':
            ret += 'S'
        elif b == 'K':
            ret += 'K'
        elif b == 'Y':
            ret += 'Y'
        elif b == 'M':
            ret += 'M'
        elif b == 'R':
            ret += 'R'
        elif b == 'W':
            ret += 'W'
        else:
            raise ValueError("Non standard base detected in {0}: {1}".format(__name__, b))
    return ret

# this if for the compliment

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

        if self.fragment:
            self.start += last_end
            self.end += last_end

    def revcomp(self):
        if self.strand == '+':
            self.strand = '-'
        elif self.strand == '-':
            self.strand = '+'
        else:
            raise ValueError("Illegal orientation {1} for contig {0}.".format(self.name, self.strand))


    def setselfid(self, id):
        self.sid = id
    def setpartnumber(self, pn):
        self.part_number = pn

    def __repr__(self):
        ret = self.name
        return ret

    def __str__(self):
        agp_start = self.start + 1  # +1 is right in this case, to increment relative to the last scaffold position.
        agp_end = self.end #- 2  # this was wrong and offset the agp file. not sure why the -2.
        if self.is_zero_length:
            return ''
        # why does scaffold_offset appear twice in the AGP line???
        # Looked at the spec and didn't see anything i could identify with this.
        if (self.type == 'U'):
            ret = 'PGA_scaffold_{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(self.sid, self.scaffold_offset, self.scaffold_offset + (agp_end - agp_start), self.part_number, self.scaffold_offset, self.type, 100,  'paired-ends')
        else:
            ret = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}'.format(self.name, self.scaffold_offset, self.scaffold_offset + (agp_end - agp_start), self.part_number, self.scaffold_offset, self.type, self.name, agp_start, agp_end, self.strand)
        return ret

printed_contigs = defaultdict(int)

class scaffold:
    def __init__ (self, ind, header, sid):

        self.contigs = []
        self.sid = sid

        self.gap_id = 0
        # i think running_offset is just keeping track of what position in the scaffold you're at
        running_offset = 1  # would like to know more about what this is doing.
        pn = 1

        for contig_index in ind:
            contig_index = int(contig_index)
            # in .assembly files, reverse orientation indicated by negative indices
            # be careful to make it positive before using to lookup!

            if contig_index < 0:
#                is_reversed_orientation = True
                is_reversed_orientation = False
                contig_index *= -1
            elif contig_index == 0:  # should never happen
                raise ValueError("AGP indices should be 1-indexed! Index 0 found for contig {0}".format(
                    header[contig_index]
                ))
            else:
                is_reversed_orientation = False

            contig_index -= 1  # translate from 1-based to 0-based system
            #skip tiny stuff
            if header[contig_index].end < 10:
                continue

            if is_reversed_orientation:
                self.contigs.append(copy.deepcopy(header[contig_index]))
#                self.contigs[-1].revcomp()
            else:
                self.contigs.append(copy.deepcopy(header[contig_index]))

            self.contigs[-1].scaffold_offset += running_offset
            self.contigs[-1].setselfid(sid)
            self.contigs[-1].setpartnumber(pn)
            running_offset += self.contigs[-1].end
            pn += 1
#            if (contig_index + 1) != int(ind[-1]):  # here translate back to 1-index for one thing
#                self.contigs.append(seq_range('gap_' + str(self.gap_id), 0, 100, 'U'))
#                self.contigs[-1].scaffold_offset += running_offset
#                running_offset += self.contigs[-1].end
#                self.contigs[-1].setpartnumber(pn)
#                self.contigs[-1].setselfid(sid)
#                self.gap_id += 1
#                pn += 1

    def __repr__(self):
        ret = 'Scaffold {0}'.format(self.sid)
        return ret

    def __str__(self):
        ret = self.__repr__() + '\n'
        ret += '\n'.join([contig for contig in self.contigs])
        return ret

    def dump_fasta(self, fasta, fasta_output_writer):
        if self.contigs.__len__() == 0:
            return
        last_contig = self.contigs[-1]

        sequence = ''
        for contig in self.contigs:
            print "adding contig sequence: {0} to scaffold: {1}".format(contig.name, self.sid)
            if contig.name not in fasta and contig.type != 'U':
                raise ValueError('FATAL: missing contig: {0}'.format(contig.name))
            if contig.type == 'U':
                sequence += 'N' * 100
                raise Exception("Failed!")
            else:
                if contig.is_zero_length:
                    print "Skipping zero-length contig: contig {0} in scaffold {1}".format(contig.name, self.sid)
                    continue
                tmpseq = fasta[contig.name][contig.start:contig.end].seq
                if contig.strand == '-123456':
                    print "revcomping sequence: {0}".format(contig.name)
                    tmpseq = reverse_compliment(tmpseq)
                sequence += tmpseq

        total_scaffold_length = len(sequence)
        print "about to print scaffold"

        printed = 0
        line_size = 80
        if len(sequence) == 0:  # to handle zero-length scaffolds- just don't print em!
            return

        fasta_output_writer.write('>{0}_bin{1}\n'.format(self.contigs[0].name, printed_contigs[self.contigs[0].name]))
        printed_contigs[self.contigs[0].name] = printed_contigs[self.contigs[0].name] + 1
        while len(sequence) - printed > line_size:
            fasta_output_writer.write('{0}\n'.format(sequence[printed:printed+line_size]))
            printed += line_size
        if len(sequence) - printed > 0:
            fasta_output_writer.write('{0}\n'.format(sequence[printed:]))
            printed += len(sequence[printed:])



parser = argparse.ArgumentParser()
parser.add_argument('-a', '--assembly', help='juicebox assembly file', required=True)
parser.add_argument('-f', '--fasta', help='the fasta file', required=True)
parser.add_argument('-p', '--prefix', help='the prefix to use for writing outputs', required=True)
args = parser.parse_args()

header      = []
scaffolds   = []
scaffold_id = 1
new_seq_end = 0
new_seq_name = ''
agp_output_file = '{0}.agp'.format(args.prefix)
fasta_output_file = '{0}.fasta'.format(args.prefix)

with open(args.assembly) as f:
    lines = f.readlines()
    sline = [x.strip().split() for x in lines]
    for items in sline:
        if items[0][0] != '>':
            for i in range(len(items)):
                print "Adding contig", items[i], "as scaffold", scaffold_id
                scaffolds.append(scaffold([items[i]], header, scaffold_id))
                scaffold_id += 1
        else:
            new_seq = seq_range(items[0][1:], items[1], items[2], 'W', new_seq_end, new_seq_name)
            new_seq_end = new_seq.end
            new_seq_name = new_seq.name
            header.append(new_seq)

with open(agp_output_file, 'w') as f:
    f.write('##agp-version 2.0\n# This file was generated by converting juicebox assembly format\n')
    for i in range(0, scaffold_id-1):
        for contig in scaffolds[i].contigs:
            if str(contig) == '':
                print "omitting debris contig {0} from AGP, zero-length.".format(contig.name)
                continue
            f.write('{0}\n'.format(contig))

input_fasta = Fasta(args.fasta)

with open(fasta_output_file, 'w') as f:
    f.write('')

fasta_output_writer = open(fasta_output_file, 'a')
for i in range(0, scaffold_id - 1):
    scaffolds[i].dump_fasta(input_fasta, fasta_output_writer)
