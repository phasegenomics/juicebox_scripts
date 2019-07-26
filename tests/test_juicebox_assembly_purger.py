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

from juicebox_scripts.juicebox_assembly_purger import filter_assembly, get_exclude

class JuiceboxPurgerTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.collateral = os.path.dirname(__file__) + "/collateral/"
        self.input_assembly = self.collateral + "test_inputs/test_scaffolds.assembly"
        self.input_exclude_file = self.collateral + "test_inputs/test_exclude.txt"
        self.output_dir = self.collateral + "test_outputs/"
        self.expected_output_dir = self.collateral + "expected_test_outputs/"

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    @classmethod
    def tearDownClass(self):
        if os.path.exists(self.output_dir):
            os.removedirs(self.output_dir)

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

    def test_get_exclude_from_list(self):
        exclude_list = ["1", "2", "3"]
        expected_exclude = set(["1", "2", "3"])
        self.assertEquals(get_exclude(exclude_list, None), expected_exclude)

    def test_get_exclude_with_duplicates(self):
        exclude_list = ["1", "2", "3", "1"]
        expected_exclude = set(["1", "2", "3"])
        self.assertEquals(get_exclude(exclude_list, None), expected_exclude)

    def test_get_exclude_from_file(self):
        expected_exclude = set(["1", "2", "3"])
        self.assertEquals(get_exclude(None, self.input_exclude_file), expected_exclude)

    def test_get_exclude_from_file_and_list(self):
        exclude_list = ["3", "4", "5"]
        expected_exclude = set(["1", "2", "3", "4", "5"])
        self.assertEquals(get_exclude(exclude_list, self.input_exclude_file), expected_exclude)

if __name__ == "__main__":
    unittest.main()
