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

import filecmp
import os
import unittest
#uncomment the below lines if you need to have the ability to run this
#file directly either from the the root directory, but you can avoid
#that by running `python -m unittest discover` within the project root
#directory instead of trying to call this file directly
#import sys
# sys.path.append('src/')
# sys.path.append('../src/')
from juicebox_scripts.juicebox_assembly_converter import JuiceboxConverter, ProcessedAssembly
from juicebox_scripts.juicebox_assembly_converter import InvalidFastaError, MissingFragmentError, \
    UnscaffoldedContigError, ZeroLengthContigError, FragmentAndContigLengthDiscrepancy

class JuiceboxConverterTestCase(unittest.TestCase):
    def setUp(self):
        self.converter = JuiceboxConverter()
        self.test_output_fasta = 'test_output.fasta'
        self.test_output_agp = 'test_output.agp'
        self.test_output_bed = 'test_output.bed'
        self.test_output_break_report = 'test_output.break_report.txt'
        
        self.collateral_dir = os.path.dirname(__file__) + '/collateral/'
        self.test_file_dir = self.collateral_dir + 'test_inputs/'
        self.test_bad_fasta_1 = self.test_file_dir + 'test_bad_1.fasta'
        self.test_bad_fasta_2 = self.test_file_dir + 'test_bad_2.fasta'
        self.test_breaks_bad_assembly = self.test_file_dir + 'test_breaks_bad.assembly'
        self.test_breaks_assembly = self.test_file_dir + 'test_breaks.assembly'
        #self.test_breaks_oversize_fragment = self.test_file_dir + 'test_breaks_oversize_fragment.assembly'
        self.test_reordered_breaks_assembly = self.test_file_dir + 'test_breaks_reordered.assembly'
        self.test_contigs_bad_assembly = self.test_file_dir + 'test_contigs_bad.assembly'
        self.test_contigs_assembly = self.test_file_dir + 'test_contigs.assembly'
        self.test_scaffolds_bad_assembly = self.test_file_dir + 'test_scaffolds_bad.assembly'
        self.test_scaffolds_assembly = self.test_file_dir + 'test_scaffolds.assembly'
        self.test_fasta = self.test_file_dir + 'test.fasta'
        
        self.expected_outputs_dir = self.collateral_dir + 'expected_test_outputs/'
        self.expected_result_breaks_straight_agp = self.expected_outputs_dir + 'expected_result_breaks_straight.agp'
        self.expected_result_breaks_straight_bed = self.expected_outputs_dir + 'expected_result_breaks_straight.bed'
        self.expected_result_breaks_straight_fasta = self.expected_outputs_dir + 'expected_result_breaks_straight.fasta'
        self.expected_result_breaks_agp = self.expected_outputs_dir + 'expected_result_breaks.agp'
        self.expected_result_breaks_bed = self.expected_outputs_dir + 'expected_result_breaks.bed'
        self.expected_result_breaks_broken_contigs_txt = self.expected_outputs_dir + 'expected_result_breaks.broken_contigs.txt'
        self.expected_result_breaks_fasta = self.expected_outputs_dir + 'expected_result_breaks.fasta'
        self.expected_result_contigs_straight_agp = self.expected_outputs_dir + 'expected_result_contigs_straight.agp'
        self.expected_result_contigs_straight_bed = self.expected_outputs_dir + 'expected_result_contigs_straight.bed'
        self.expected_result_contigs_straight_fasta = self.expected_outputs_dir + 'expected_result_contigs_straight.fasta'
        self.expected_result_contigs_agp = self.expected_outputs_dir + 'expected_result_contigs.agp'
        self.expected_result_contigs_bed = self.expected_outputs_dir + 'expected_result_contigs.bed'
        self.expected_result_contigs_broken_contigs_txt = self.expected_outputs_dir + 'expected_result_contigs.broken_contigs.txt'
        self.expected_result_contigs_fasta = self.expected_outputs_dir + 'expected_result_contigs.fasta'
        self.expected_result_scaffolds_agp = self.expected_outputs_dir + 'expected_result_scaffolds.agp'
        self.expected_result_scaffolds_bed= self.expected_outputs_dir + 'expected_result_scaffolds.bed'
        self.expected_result_scaffolds_fasta = self.expected_outputs_dir + 'expected_result_scaffolds.fasta'
        

    def tearDown(self):
        if os.path.exists(self.test_output_fasta):
            os.remove(self.test_output_fasta)
        if os.path.exists(self.test_output_agp):
            os.remove(self.test_output_agp)
        if os.path.exists(self.test_output_bed):
            os.remove(self.test_output_bed)
        if os.path.exists(self.test_output_break_report):
            os.remove(self.test_output_break_report)
    
    def read_file_lines(self, fname):
        with open(fname) as f:
            return f.readlines()

    def test_make_contigs(self):
        contigs = self.converter.process(self.test_fasta, self.test_contigs_assembly)
        expected_contigs_fasta = self.read_file_lines(self.expected_result_contigs_fasta)
        expected_contigs_agp = self.read_file_lines(self.expected_result_contigs_agp)
        expected_contigs_bed = self.read_file_lines(self.expected_result_contigs_bed)
        self.assertEqual(contigs.fasta(), expected_contigs_fasta)
        self.assertEqual(contigs.agp(), expected_contigs_agp)
        self.assertEqual(contigs.bed(), expected_contigs_bed)
    
    def test_make_scaffolds(self):
        scaffolds = self.converter.process(self.test_fasta, self.test_scaffolds_assembly)
        expected_scaffolds_fasta = self.read_file_lines(self.expected_result_scaffolds_fasta)
        expected_scaffolds_agp = self.read_file_lines(self.expected_result_scaffolds_agp)
        expected_scaffolds_bed = self.read_file_lines(self.expected_result_scaffolds_bed)
        self.assertEqual(scaffolds.fasta(), expected_scaffolds_fasta)
        self.assertEqual(scaffolds.agp(), expected_scaffolds_agp)
        self.assertEqual(scaffolds.bed(), expected_scaffolds_bed)
    
    def test_make_breaks(self):
        breaks = self.converter.process(self.test_fasta, self.test_breaks_assembly)
        expected_breaks_fasta = self.read_file_lines(self.expected_result_breaks_fasta)
        expected_breaks_agp = self.read_file_lines(self.expected_result_breaks_agp)
        expected_breaks_bed = self.read_file_lines(self.expected_result_breaks_bed)
        self.assertEqual(breaks.fasta(), expected_breaks_fasta)
        self.assertEqual(breaks.agp(), expected_breaks_agp)
        self.assertEqual(breaks.bed(), expected_breaks_bed)
    
    def test_make_contigs_in_contig_mode(self):
        contigs = self.converter.process(self.test_fasta, self.test_contigs_assembly, contig_mode=True)
        expected_contigs_fasta = self.read_file_lines(self.expected_result_contigs_straight_fasta)
        expected_contigs_agp = self.read_file_lines(self.expected_result_contigs_straight_agp)
        expected_contigs_bed = self.read_file_lines(self.expected_result_contigs_straight_bed)
        self.assertEqual(contigs.fasta(), expected_contigs_fasta)
        self.assertEqual(contigs.agp(), expected_contigs_agp)
        self.assertEqual(contigs.bed(), expected_contigs_bed)

    def test_make_scaffolds_in_contig_mode(self):
        contigs = self.converter.process(self.test_fasta, self.test_scaffolds_assembly, contig_mode=True)
        #the scaffolds assembly with no breaks should do the same thing as the contigs assembly
        expected_contigs_fasta = self.read_file_lines(self.expected_result_contigs_straight_fasta)
        expected_contigs_agp = self.read_file_lines(self.expected_result_contigs_straight_agp)
        expected_contigs_bed = self.read_file_lines(self.expected_result_contigs_straight_bed)
        self.assertEqual(contigs.fasta(), expected_contigs_fasta)
        self.assertEqual(contigs.agp(), expected_contigs_agp)
        self.assertEqual(contigs.bed(), expected_contigs_bed)
    
    def test_make_breaks_in_contig_mode(self):
        contigs = self.converter.process(self.test_fasta, self.test_breaks_assembly, contig_mode=True)
        expected_contigs_fasta = self.read_file_lines(self.expected_result_breaks_straight_fasta)
        expected_contigs_agp = self.read_file_lines(self.expected_result_breaks_straight_agp)
        expected_contigs_bed = self.read_file_lines(self.expected_result_breaks_straight_bed)
        self.assertEqual(contigs.fasta(), expected_contigs_fasta)
        self.assertEqual(contigs.agp(), expected_contigs_agp)
        self.assertEqual(contigs.bed(), expected_contigs_bed)
    
    def test_make_reordered_contigs_in_contig_mode(self):
        contigs = self.converter.process(self.test_fasta, self.test_reordered_breaks_assembly, contig_mode=True)
        expected_contigs_fasta = self.read_file_lines(self.expected_result_breaks_straight_fasta)
        expected_contigs_agp = self.read_file_lines(self.expected_result_breaks_straight_agp)
        expected_contigs_bed = self.read_file_lines(self.expected_result_breaks_straight_bed)
        self.assertEqual(contigs.fasta(), expected_contigs_fasta)
        self.assertEqual(contigs.agp(), expected_contigs_agp)
        self.assertEqual(contigs.bed(), expected_contigs_bed)

    def test_no_breaks_break_report(self):
        contigs = self.converter.process(self.test_fasta, self.test_contigs_assembly)
        expected_no_break_report = self.read_file_lines(self.expected_result_contigs_broken_contigs_txt)
        self.assertEqual(contigs.break_report(), expected_no_break_report)
    
    def test_breaks_break_report(self):
        contigs = self.converter.process(self.test_fasta, self.test_breaks_assembly)
        expected_break_report = self.read_file_lines(self.expected_result_breaks_broken_contigs_txt)
        self.assertEqual(contigs.break_report(), expected_break_report)
    
    def test_write_fasta(self):
        contigs = self.converter.process(self.test_fasta, self.test_contigs_assembly)
        contigs.write_fasta(self.test_output_fasta)
        self.assertTrue(filecmp.cmp(self.expected_result_contigs_fasta, self.test_output_fasta))
    
    def test_write_agp(self):
        contigs = self.converter.process(self.test_fasta, self.test_contigs_assembly)
        contigs.write_agp(self.test_output_agp)
        self.assertTrue(filecmp.cmp(self.expected_result_contigs_agp, self.test_output_agp))
    
    def test_write_bed(self):
        contigs = self.converter.process(self.test_fasta, self.test_contigs_assembly)
        contigs.write_bed(self.test_output_bed)
        self.assertTrue(filecmp.cmp(self.expected_result_contigs_bed, self.test_output_bed))
    
    def test_write_break_report(self):
        contigs = self.converter.process(self.test_fasta, self.test_breaks_assembly)
        contigs.write_break_report(self.test_output_break_report)
        self.assertTrue(filecmp.cmp(self.expected_result_breaks_broken_contigs_txt, self.test_output_break_report))
    
    def test_bad_contigs(self):
        with self.assertRaises(ZeroLengthContigError):
            self.converter.process(self.test_fasta, self.test_contigs_bad_assembly)
    
    def test_bad_scaffolds(self):
        with self.assertRaises(UnscaffoldedContigError):
            self.converter.process(self.test_fasta, self.test_scaffolds_bad_assembly)
    
    def test_bad_breaks(self):
        with self.assertRaises(MissingFragmentError):
            self.converter.process(self.test_fasta, self.test_breaks_bad_assembly)

    #def test_oversize_fragment(self):
    #    with self.assertRaises(FragmentAndContigLengthDiscrepancy):
    #        self.converter.process(self.test_fasta, self.test_breaks_oversize_fragment)

    def test_bad_fasta(self):
        with self.assertRaises(InvalidFastaError):
            self.converter.process(self.test_bad_fasta_1, self.test_contigs_assembly)
        with self.assertRaises(InvalidFastaError):
            self.converter.process(self.test_bad_fasta_2, self.test_contigs_assembly)    


if __name__ == '__main__':
    unittest.main()
