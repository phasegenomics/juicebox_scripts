#!/usr/bin/env python
'''
Shawn Sullivan
October 31, 2018
Phase Genomics

juicebox_scripts/test_juicebox_assembly_converter.py

This file contains unit tests for functions of the juicebox_assembly_converter.py script.

Copyright 2018, Phase Genomics Inc. All rights reserved.

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
from _collections import defaultdict

class JuiceboxConverter:
    def __init__(self):
        pass
    
    def process(self, fasta, assembly, contig_mode=False):
        sequences = self._read_fasta(fasta)
        assembly_map, scaffolds = self._read_assembly(assembly, contig_mode=contig_mode)
        sequences = self._add_breaks(sequences, assembly_map)
        return ProcessedAssembly(sequences, assembly_map, scaffolds)
    
    def _read_fasta(self, fasta):
        sequences = dict()
        active_seq = None
        with open(fasta) as f:
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    continue
                if line[0] == '>':
                    active_seq = line[1:]
                    if active_seq in sequences.keys():
                        raise InvalidFastaError('Fasta {0} contains multiple contigs named {1}'.format(fasta, active_seq))
                    sequences[active_seq] = ''
                elif active_seq is not None:
                    sequences[active_seq] += line
                else:
                     raise InvalidFastaError('Fasta {0} does not begin with a contig name'.format(fasta))
        return sequences
    
    def _read_assembly(self, assembly, contig_mode=False):
        assembly_map = list()
        scaffolds = list()
        unscaffolded_contigs = list()
        with open(assembly) as f:
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    continue
                if line[0] == '>':
                    # >cname index len
                    tokens = line[1:].split()
                    if int(tokens[1]) != len(assembly_map) + 1:
                        raise MissingFragmentError('Assembly {0} is missing the sequence for index {1}'.format(assembly, len(assembly_map) + 1))
                    if int(tokens[2]) == 0:
                        raise ZeroLengthContigError('Assembly {0} lists contig {1} as zero length'.format(assembly, tokens[0]))
                    assembly_map.append((tokens[0], tokens[2]))
                    unscaffolded_contigs.append(tokens[0])
                else:
                    if contig_mode:
                        for contig in assembly_map:
                            scaffolds.append([(contig[0], contig[1], '+', contig_mode)])
                            unscaffolded_contigs.remove(contig[0])
                        break
                    else:
                        scaffold = list()
                        contigs = line.split()
                        for contig in contigs:
                            index = abs(int(contig)) - 1
                            strand = '+' if int(contig) > 0 else '-'
                            scaffold.append((assembly_map[index][0], assembly_map[index][1], strand, contig_mode))
                            unscaffolded_contigs.remove(assembly_map[index][0])
                        scaffolds.append(scaffold)
        if len(unscaffolded_contigs) != 0:
            raise UnscaffoldedContigError('Contigs are not included in scaffolding output: {0}'.format(unscaffolded_contigs))
        return assembly_map, scaffolds
    
    def _add_breaks(self, sequences, assembly_map):
        sequence_offsets = defaultdict(int)
        new_sequences = dict() 
        for fragment in assembly_map:
            fragment_name = fragment[0]
            fragment_size = int(fragment[1])
            if '::fragment' in fragment_name:
                orig_contig = fragment_name.split('::fragment')[0]
                new_sequences[fragment_name] = sequences[orig_contig][sequence_offsets[orig_contig]:sequence_offsets[orig_contig]+fragment_size]
                sequence_offsets[orig_contig] += fragment_size
            else:
                new_sequences[fragment_name] = sequences[fragment_name]
        return new_sequences

class ProcessedAssembly:
    def __init__(self, sequences, assembly_map, scaffolds):
        self.sequences = sequences
        self.assembly_map = assembly_map
        self.scaffolds = scaffolds
        self.contig_mode = scaffolds[0][0][3]
        self.gap_size = 100
    
    def write_fasta(self, outfile):
        self._write_file(outfile, self.fasta())
    
    def write_agp(self, outfile):
        self._write_file(outfile, self.agp())
    
    def write_bed(self, outfile):
        self._write_file(outfile, self.bed())
    
    def write_break_report(self, outfile):
        self._write_file(outfile, self.break_report())
    
    def fasta(self):
        ret = list()
        for index, scaffold in enumerate(self.scaffolds):
            scaffold_name = self._make_scaffold_name(index+1, scaffold)
            ret.append('>' + scaffold_name + '\n')
            seq = ''
            for contig in scaffold:
                seq += self.sequences[contig[0]] if contig[2] == '+' else self._reverse_complement(self.sequences[contig[0]]) 
                if contig != scaffold[-1]:
                    seq += 'n' * self.gap_size
            ret += self._chunk_sequence(seq)
        ret[-1] = ret[-1].strip()
        return ret
    
    def agp(self):
        ret = list()
        ret.append('##agp-version 2.0\n')
        ret.append('# This file was generated by converting juicebox assembly format\n')
        for index, scaffold in enumerate(self.scaffolds):
            scaffold_name = self._make_scaffold_name(index+1, scaffold)
            if not self.contig_mode:
                scaffold_name = scaffold_name.split()[0]
            offset_coord = 1
            part_number = 1
            for contig in scaffold:
                ret.append(self._make_agp_line(scaffold_name, contig, offset_coord, part_number))
                offset_coord += int(contig[1])
                part_number += 1
                if contig != scaffold[-1]:
                    ret.append(self._make_agp_gap_line(scaffold_name, offset_coord, part_number))
                    offset_coord += self.gap_size
                    part_number += 1
        ret[-1] = ret[-1].strip()
        return ret 
    
    def bed(self):
        ret = list()
        ret.append('##bed file\n')
        ret.append('# This file was generated by converting juicebox assembly format\n')
        gap_number = 1
        for index, scaffold in enumerate(self.scaffolds):
            scaffold_name = self._make_scaffold_name(index+1, scaffold)
            if not self.contig_mode:
                scaffold_name = scaffold_name.split()[0]
            offset_coord = 0
            for contig in scaffold:
                ret.append(self._make_bed_line(scaffold_name, contig, offset_coord))
                offset_coord += int(contig[1])
                if contig != scaffold[-1]:
                    ret.append(self._make_bed_gap_line(scaffold_name, offset_coord, gap_number))
                    offset_coord += self.gap_size
                    gap_number += 1
        ret[-1] = ret[-1].strip()
        return ret
    
    def break_report(self):
        ret = list()
        break_count = 0
        broken_orig_contigs = set()
        break_offsets = defaultdict(int)
        for contig in self.assembly_map:
            fragment_name = contig[0]
            fragment_size = contig[1]
            if '::fragment' in fragment_name:
                orig_contig = fragment_name.split('::fragment')[0]
                break_start = break_offsets[orig_contig]
                break_end = break_start + int(fragment_size)
                line = '\t'.join([
                                    orig_contig,
                                    fragment_name,
                                    str(break_start),
                                    str(break_end),
                                    fragment_size
                                ])
                ret.append(line + '\n')
                if '::debris' in fragment_name:
                    break_count += 1
                broken_orig_contigs.add(orig_contig)
                break_offsets[orig_contig] += int(fragment_size)
        ret.insert(0, '#orig_contig\tfragment\tbreak_start\tbreak_end\tfragment_len\n')
        ret.insert(0, '#{0} total breaks in {1} contigs\n'.format(break_count, len(broken_orig_contigs)))
        ret[-1] = ret[-1].strip()
        return ret
    
    def _write_file(self, outfile, contents):
        with open(outfile, 'w') as f:
            f.writelines(contents)
    
    def _make_scaffold_name(self, index, scaffold):
        if self.contig_mode:
            scaffold_name = '{0}'.format(scaffold[0][0])
        else:
            contig_count = len(scaffold)
            scaffold_length = 0
            for contig in scaffold:
                scaffold_length += int(contig[1])
                if contig != scaffold[-1]:
                    scaffold_length += self.gap_size
            scaffold_name = 'PGA_scaffold_{0} {1}_contigs length_{2}'.format(index,
                                                                             contig_count,
                                                                             scaffold_length)
        return scaffold_name
    
    def _chunk_sequence(self, sequence, line_len=80):
        chunked_sequence = list()
        len_added = 0
        while len(sequence) - len_added > line_len:
            chunked_sequence.append(sequence[len_added:len_added+line_len] + '\n')
            len_added += line_len
        if len(sequence) - len_added > 0: 
            chunked_sequence.append(sequence[len_added:] + '\n')
        return chunked_sequence

    def _make_agp_line(self, scaffold_name, contig, offset_coord, part_number):
        scaff_start_coord = str(offset_coord)
        scaff_end_coord = str(offset_coord + int(contig[1]) - 1)
        elt_type = 'W'
        contig_name = contig[0]
        contig_start_coord = '1'
        contig_end_coord = contig[1]
        line = '\t'.join([
                            scaffold_name,
                            scaff_start_coord,
                            scaff_end_coord,
                            str(part_number),
                            elt_type,
                            contig_name,
                            contig_start_coord,
                            contig_end_coord
                        ])
        return line + '\n'
    
    def _make_agp_gap_line(self, scaffold_name, offset_coord, part_number):
        gap_start_coord = str(offset_coord)
        gap_end_coord = str(offset_coord + self.gap_size - 1)
        elt_type = 'U'
        gap_size = str(self.gap_size)
        gap_type = 'paired-ends'
        line = '\t'.join([
                        scaffold_name,
                        gap_start_coord,
                        gap_end_coord,
                        str(part_number),
                        elt_type,
                        gap_size,
                        gap_type
                    ])
        return line + '\n'     
    
    def _make_bed_line(self, scaffold_name, contig, offset_coord):
        start_coord = str(offset_coord)
        end_coord = str(offset_coord + int(contig[1]))
        contig_name = contig[0]
        contig_len = contig[1]
        strand = contig[2]
        line = '\t'.join([
                            scaffold_name,
                            start_coord,
                            end_coord,
                            contig_name,
                            contig_len,
                            strand
                        ])
        return line + '\n'
    
    def _make_bed_gap_line(self, scaffold_name, offset_coord, gap_number):
        start_coord = str(offset_coord)
        end_coord = str(offset_coord + self.gap_size)
        gap_name = 'pg_gap_{0}'.format(gap_number)
        line = '\t'.join([
                            scaffold_name,
                            start_coord,
                            end_coord,
                            gap_name
                        ])
        return line + '\n'
    
    def _reverse_complement(self, sequence):
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
    


class ZeroLengthContigError(ValueError):
    pass

class UnscaffoldedContigError(ValueError):
    pass

class MissingFragmentError(ValueError):
    pass

class InvalidFastaError(ValueError):
    pass




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--assembly', help='juicebox assembly file', required=True)
    parser.add_argument('-f', '--fasta', help='the fasta file', required=True)
    parser.add_argument('-p', '--prefix', help='the prefix to use for writing outputs. '\
                        'Default: the assembly file, minus the file extension', default=None)
    parser.add_argument('-c', '--contig_mode', action='store_true', help='ignore scaffold '\
                        'specification and just output contigs. useful when only trying to '\
                        'obtain a fasta reflecting juicebox breaks', default=False)
    args = parser.parse_args()
    
    assembly = args.assembly
    fasta = args.fasta
    prefix = args.prefix
    if prefix is None:
        import os.path.splitext
        splitext(args.assembly)[0]
    contig_mode = args.contig_mode
    
    processed_assembly = JuiceboxConverter().process(fasta, assembly, contig_mode)
    processed_assembly.write_agp(prefix + '.agp')
    processed_assembly.write_bed(prefix + '.bed')
    processed_assembly.write_break_report(outfile)(prefix + '.break_report.txt')
    processed_assembly.write_fasta(outfile)(prefix + '.fasta')

