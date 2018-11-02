#!/usr/bin/env python
'''
Shawn Sullivan
October 31, 2018
Phase Genomics

juicebox_scripts/test_juicebox_assembly_converter.py

This file contains unit tests for functions of the
juicebox_assembly_converter.py script.

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
    '''The JuiceboxConverter class offers methods to read in a Juicebox
    .assembly file and the accompanying .fasta file, and generates a
    flexible data structure (ProcessedAssembly) which can be used to
    output additional file formats. Via the contig_mode flag, it
    supports simply reading in a .assembly file and reflecting the
    scaffold structures present in that file, as well as reading in the
    .assembly file and only outputting the contigs it describes, which
    is useful when only desiring contigs which have been broken using
    Juicebox.
    
    Attributes:
        None
    '''
    def __init__(self):
        pass
    
    def process(self, fasta, assembly, contig_mode=False):
        '''Read in a .assembly file and .fasta file, generating a
        ProcessedAssembly reflecting them.
        
        Args:
            fasta (str): path to the fasta corresponding to the assembly
            assembly (str): path to the .assembly file generated by
                Juicebox containing the results of using it
            contig_mode (bool) [ptional]: instead of generating a
                ProcessedAssembly reflecting the scaffolds in the bottom
                portion of the .assembly file, only reflect the contigs
                described in the top portion of the file. Useful when
                you have broken contigs in Juicebox and only want to get
                an assembly reflecting those breaks, without
                scaffolding. Default: False
        
        Returns:
            ProcessedAssembly: a ProcessedAssembly object reflecting the
                inputs
        '''
        
        sequences = self._read_fasta(fasta)
        assembly_map, scaffolds = self._read_assembly(assembly, contig_mode=contig_mode)
        sequences = self._add_breaks(sequences, assembly_map)
        return ProcessedAssembly(sequences, assembly_map, scaffolds)
    
    def _read_fasta(self, fasta):
        '''Read in a .fasta file and return a dictionary mapping the
        sequence names to their sequences.
        
        Args:
            fasta (str): path to the fasta containing the sequences you
                wish to read
        
        Returns:
            dict[str:str]: dict mapping contig/sequence names present in
                the .fasta to their sequences
        '''
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
        '''Read in a .assembly file and return two lists reflecting its
        contents, one for the contig list in the top of the file and the
        other the scaffolds listed at the bottom of the file.
        
        Args:
            assembly (str): path to the Juicebox .assembly file to
                process
            contig_mode (bool): whether to process the assembly in a
                manner which reflects the scaffolds in the bottom of the
                file, or just lists the contigs without scaffolding in
                the scaffold list. Default: False
        
        Returns:
            (list, list): tuple of two lists reflecting the .assembly
                file
                list[(str, str), (str, str), ...]: information about the
                    contigs shown in the top section of the file. The
                    list includes a tuple for each contig listed in the
                    .assembly file, in the order they are listed in the
                    top section of the assembly file.
                    str: the contig name
                    str: the contig length
                list[list[(str, str, str, bool), (str, str, str, bool), ...]:
                    information about the scaffolds shown in the bottom
                    section of the file, or about the contigs in the top
                    section of the file if contig_mode is True. Each
                    individual scaffold is a list containing the contigs
                    on that scaffold. Ordering within these lists matches
                    the ordering of contigs shown in the scaffolds
                    section at the bottom of the .assembly file if
                    contig_mode is False. If contig_mode is True, each
                    scaffold only contains a single contig, and they
                    are ordered in the same order they are present
                    in the top section of the .assembly file.
                    str: the contig name
                    str: the contig length
                    str: the contig strand (+ or -)
                    bool: whether the contig was placed using contig_mode
                        or not
        '''
            
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
        '''Introduce breaks into an assembly_map, as generated by
        _read_assembly, where indicated in the .assembly file. Does not
        mutate sequences - instead it returns a brand new sequences
        dict.
        
        Args:
            sequences (dict[str:str]): dictionary mapping contig/sequence
                names to their sequences
            assembly_map (list[(str, str)]: list containing contig
                information from the .assembly file, as generated by
                _read_assembly
        
        Returns:
            dict[str:str]: a brand new dictionary mapping contig/
                sequence names to their sequences, after performing
                breaks suggested by the contig names and information in
                assembly_map
        '''
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
    '''The ProcessedAssembly class represents the results of
    JuiceboxConverter parsing a .fasta file and a .assembly file. It
    provides methods to output those processed results in several formats
    (.fasta, .agp, .bed) and also can generate a report describing any
    breaks that were introduced into the assembly by the .assembly file.
    
    Attributes:
        sequences (dict[str:str]): dictionary mapping contig names to
            their sequences
        assembly_map (list[(str, str)]: list containing contig
            information from the .assembly file, as generated by
            JuiceboxConverter._read_assembly
        scaffolds (list[list[(str, str, str, bool)): list containing
            scaffold information from the .assembly file, as generated by
            JuiceboxConverter._read_assembly
    '''
    def __init__(self, sequences, assembly_map, scaffolds):
        self.sequences = sequences
        self.assembly_map = assembly_map
        self.scaffolds = scaffolds
        self.contig_mode = scaffolds[0][0][3]
        self.gap_size = 100
    
    def write_fasta(self, outfile):
        '''Write FASTA output to a specified file
        
        Args:
            outfile (str): path to the file to write to
        '''
        self._write_file(outfile, self.fasta())
    
    def write_agp(self, outfile):
        '''Write AGP output to a specified file
        
        Args:
            outfile (str): path to the file to write to
        '''
        self._write_file(outfile, self.agp())
    
    def write_bed(self, outfile):
        '''Write BED output to a specified file
        
        Args:
            outfile (str): path to the file to write to
        '''
        self._write_file(outfile, self.bed())
    
    def write_break_report(self, outfile):
        '''Write a report about breaks shown in the .assembly file to a
        specified file
        
        Args:
            outfile (str): path to the file to write to
        '''
        self._write_file(outfile, self.break_report())
    
    def fasta(self):
        '''Generate a FASTA format representation of the
        ProcessedAssembly
        
        Returns:
            list[str]: list of strings representing the ProcessedAssembly
                in FASTA format
        '''
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
        '''Generate an AGP format representation of the
        ProcessedAssembly
        
        Returns:
            list[str]: list of strings representing the ProcessedAssembly
                in AGP format
        '''
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
        '''Generate a BED format representation of the
        ProcessedAssembly
        
        Returns:
            list[str]: list of strings representing the ProcessedAssembly
                in BED format
        '''
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
        '''Generate a report about breaks shown in the .assembly file 
        
        Returns:
            list[str]: list of strings summarizing the breaks present in
                the ProcessedAssembly
        '''
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
        '''Write contents to a specified file
        
        Args:
            outfile (str): path to the file to write to
            contents (list[str]): list of strings to write to a file
        '''
        with open(outfile, 'w') as f:
            f.writelines(contents)
    
    def _make_scaffold_name(self, index, scaffold):
        '''Construct a string that shows the proper name of a scaffold.
        If running in contig_mode, just use the name of the contig
        instead.
        
        Args:
            index (int): the index of the scaffold
            scaffold list[(str, str, str, bool)]: a list of tuples
                reflecting the scaffold. Each tuple conforms to the
                schema used in JuiceboxConverter._read_assembly
                str: the contig name
                str: the contig length
                str: the contig strand (+ or -)
                bool: whether the contig was placed using contig_mode
                    or not
        '''
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
        '''Process a long string representation of a sequence into a list
        of chunks that are line_len long (the last string may be shorter)
        
        Args:
            sequence (str): the sequence to chunk
            line_len (int) [optional]: the length of the lines to break
                the sequence into. The last line may be shorter as it
                will be the remainder after chunking the rest of the
                sequence. Default: 80
        
        Returns:
            list(str): a list of strings reflecting the input sequence
                broken into chunks
        '''
        chunked_sequence = list()
        len_added = 0
        while len(sequence) - len_added > line_len:
            chunked_sequence.append(sequence[len_added:len_added+line_len] + '\n')
            len_added += line_len
        if len(sequence) - len_added > 0: 
            chunked_sequence.append(sequence[len_added:] + '\n')
        return chunked_sequence

    def _make_agp_line(self, scaffold_name, contig, offset_coord, part_number):
        '''Make a line for an AGP file reflecting the positioning of a
        single contig.
        
        Args:
            scaffold_name (str): the name of the scaffold the contig is
                placed on
            contig ((str, str)): tuple containing the name of the contig
                and its length as strings. Matches the representation
                of contigs in an assembly_map as generated by
                JuiceboxConverter._read_assembly
            offset_coord (int): the offset of the contig in the scaffold
            part_number (int): the part number of the contig in the
                scaffold
        
        Returns:
            str: an AGP-format line reflecting the contig and other
                inputs
        '''
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
        '''Make a line for an AGP file reflecting the positioning of a
        gap in the scaffold.
        
        Args:
            scaffold_name (str): the name of the scaffold the gap is
                placed on
            offset_coord (int): the offset of the gap in the scaffold
            part_number (int): the part number of the gap in the
                scaffold
        
        Returns:
            str: an AGP-format line reflecting the gap and other
                inputs
        '''
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
        '''Make a line for a BED file reflecting the positioning of a
        single contig.
        
        Args:
            scaffold_name (str): the name of the scaffold the contig is
                placed on
            contig ((str, str)): tuple containing the name of the contig
                and its length as strings. Matches the representation
                of contigs in an assembly_map as generated by
                JuiceboxConverter._read_assembly
            offset_coord (int): the offset of the contig in the scaffold
        
        Returns:
            str: a BED-format line reflecting the contig and other
                inputs
        '''
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
        '''Make a line for a BED file reflecting the positioning of a
        gap in the scaffold.
        
        Args:
            scaffold_name (str): the name of the scaffold the gap is
                placed on
            offset_coord (int): the offset of the gap in the scaffold
            gap_number (int): the index of the gap in the scaffold
        
        Returns:
            str: a BED-format line reflecting the gap and other
                inputs
        '''
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
        '''Reverse complement a sequence
        
        Args:
            sequence (str): the sequence to reverse complement
        
        Returns:
            str: the reverse complement of the sequence
        '''
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
    '''An Error caused when a contig is listed as having zero length in
    a .assembly file'''
    pass

class UnscaffoldedContigError(ValueError):
    '''An Error caused when a contig is listed in a .assembly file, but
    is not placed anywhere in a scaffold.'''
    pass

class MissingFragmentError(ValueError):
    '''An Error caused when a contig or fragment is missing from the list
    of contigs at the start of the file'''
    pass

class InvalidFastaError(ValueError):
    '''An Error caused when an invalid .fasta file is attempted to be
    read.'''
    pass


if __name__ == '__main__':
    import argparse
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
