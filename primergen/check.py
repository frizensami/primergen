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
import editdistance
import itertools
import functools
from polyleven import levenshtein


def are_primers_valid(primers):
    """
    Returns a boolean indicating of a list of primers all mutually pass our conditions to be a valid primer library.
    Exits early if any of them fail.
    """
    # Run the easy checks first
    easy_valids = map(is_len_gc_valid, primers)
    is_all_valid_easy = all(easy_valids)
    if not is_all_valid_easy:
        print(
            f"PRIMER LIBRARY FAILURE: >= 1 primers have wrong length or wrong GC proportion"
        )
        return False

    # Run the n-choose-2 valids now
    for (p1, p2) in itertools.combinations(primers, 2):
        if not is_primer_pair_valid(p1, p2, MIN_EDIT_DISTANCE):
            print(
                f"PRIMER LIBRARY FAILURE: {p1} and {p2} are not {MIN_EDIT_DISTANCE} edit distance apart!"
            )
            return False

    return True


def get_valid_primers(primers):
    """
    Cuts down an input list of primers into a list that only has valid primers.

    Problem to solve: how do we deal with the combinations? Ideally we would return the maximum subset of primers satisfying the edit distance, but
    this smells of an NP-hard problem. May have to just do greedy in some way? No idea.
    """
    pass


def is_right_length(primer) -> bool:
    return len(primer) == 20


def is_gc_valid(primer) -> bool:
    gc_percentage = GC(primer)
    return MIN_CG_CONTENT <= gc_percentage and gc_percentage <= MAX_CG_CONTENT


def is_len_gc_valid(primer) -> bool:
    return is_right_length(primer) and is_gc_valid(primer)


@functools.lru_cache(maxsize=100000)
def is_primer_pair_valid(p1, p2, limit=None):
    """
    Memoized to avoid computing edit distance multiple times for identical string pairs
    Add a limit to use the possibly faster polyleven library (only care about edit distances up to a point, then return early I expect)
    TODO: is this memoized on value or object ID itself?
    """
    dist = -1
    if limit:
        dist = levenshtein(p1, p2)
    else:
        dist = editdistance.eval(p1, p2)
    return dist > MIN_EDIT_DISTANCE
