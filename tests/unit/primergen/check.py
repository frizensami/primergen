#!/usr/bin/env python3

from primergen.check import is_right_length, is_gc_valid, are_primers_valid
from Bio.Seq import Seq
import unittest


class TestCheckMethods(unittest.TestCase):
    def test_gc_content(self):
        s = Seq("ATGC")
        self.assertTrue(is_gc_valid(s))

    def test_gc_content_toolittle(self):
        s = Seq("ATTC")
        self.assertFalse(is_gc_valid(s))

    def test_gc_content_toomuch(self):
        s = Seq("AGGC")
        self.assertFalse(is_gc_valid(s))

    def test_length_correct(self):
        s = Seq("ATGCA" * 4)
        self.assertTrue(is_right_length(s))

    def test_length_toolong(self):
        s = Seq("ATGCA" * 5)
        self.assertFalse(is_right_length(s))

    def test_length_tooshort(self):
        s = Seq("ATGCA" * 3)
        self.assertFalse(is_right_length(s))

    def test_all_primers_ok(self):
        s1 = Seq("ATGC" * 5)
        s2 = Seq("GCTA" * 5)
        s3 = Seq("CGTA" * 5)
        self.assertTrue(are_primers_valid([s1, s2, s3]))

    def test_all_primers_ok_two_identical_fail(self):
        s1 = Seq("ATGC" * 5)
        s2 = Seq("ATGC" * 5)
        self.assertFalse(are_primers_valid([s1, s2]))


if __name__ == "__main__":
    unittest.main()
