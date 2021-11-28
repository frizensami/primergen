#!/usr/bin/env python3

"""
Checks whether a list of primers satisfies these constraints:
1. All same length L = 20
2. Balanced CG content for each primer, between 45% and 55%: 0.45 <= (G + C) / (A + G + T + C) <= 0.55
3. Edit distance between any two primers is at least 0.4L (8)
"""

PRIMER_LENGTH = 20
MIN_EDIT_DISTANCE = 0.4 * PRIMER_LENGTH
MIN_CG_CONTENT = 45
MAX_CG_CONTENT = 55

from Bio.SeqUtils import GC


def are_primers_valid(primers: [str]) -> bool:
    pass


def is_gc_valid(primer) -> bool:
    gc_percentage = GC(primer)
    return MIN_CG_CONTENT <= gc_percentage and gc_percentage <= MAX_CG_CONTENT
