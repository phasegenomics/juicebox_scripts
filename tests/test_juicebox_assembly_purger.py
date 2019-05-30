#!/usr/bin/env python
'''
Brad Nelson
May 30, 2019
Phase Genomics

tests/test_juicebox_assembly_purger.py

This file contains unit tests for functions of the juicebox_assembly_purger.py script.

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

import os
import filecmp
import unittest

from juicebox_scripts.juicebox_assembly_purger import filter_assembly

class JuiceboxPurgerTestCase(unittest.TestCase):
    def setUp(self):
        self.collateral = os.path.dirname(__file__) + "/collateral/"
        self.input_assembly = self.collateral + "test_inputs/test_scaffolds.assembly"
        self.output_dir = self.collateral + "test_outputs/"
        self.expected_output_dir = self.collateral + "expected_test_outputs/"

    def tearDown(self):
        outfiles = ["test_filter_assembly_null.assembly",
                    "test_filter_assembly_minus_one.assembly",
                    "test_filter_assembly_minus_three.assembly",
                    "test_filter_assembly_missing_contig.assembly"]
        for outfile in outfiles:
            outfile = self.output_dir + outfile
            if os.path.exists(outfile):
                os.remove(outfile)

    def test_filter_assembly_null(self):
        exclude = set([])
        outfile = self.output_dir + "test_filter_assembly_null.assembly"
        filter_assembly(exclude, self.input_assembly, outfile, logging="silent")
        self.assertTrue(filecmp.cmp(outfile, self.input_assembly))

    def test_filter_assembly_minus_one(self):
        exclude = set(["contig_1_len_29"])
        outfile = self.output_dir + "test_filter_assembly_minus_one.assembly"
        expected_output = self.expected_output_dir + "test_filter_assembly_minus_one.assembly"
        filter_assembly(exclude, self.input_assembly, outfile, logging="silent")
        self.assertTrue(filecmp.cmp(outfile, expected_output))

    def test_filter_assembly_minus_three(self):
        exclude = set(["contig_1_len_29", "contig_5_len_25", "contig_8_len_16"])
        outfile = self.output_dir + "test_filter_assembly_minus_three.assembly"
        expected_output = self.expected_output_dir + "test_filter_assembly_minus_three.assembly"
        filter_assembly(exclude, self.input_assembly, outfile, logging="silent")
        self.assertTrue(filecmp.cmp(outfile, expected_output))

    def test_filter_assembly_missing_contig(self):
        exclude = set(["missing"])
        outfile = self.output_dir + "test_filter_assembly_missing_contig.assembly"
        with self.assertRaises(ValueError):
            filter_assembly(exclude, self.input_assembly, outfile, logging="silent")
