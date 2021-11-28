#!/usr/bin/env python3

from primergen.check import is_gc_valid
from Bio.Seq import Seq
import unittest


class TestCheckMethods(unittest.TestCase):
    def test_gc_content(self):
        s = Seq("ATGC")
        self.assertTrue(is_gc_valid(s))


if __name__ == "__main__":
    unittest.main()
