#!/usr/bin/env python3

import time
import os
import random

DATA_FOLDER = "data/"
ALPHABET = ["A", "T", "C", "G"]
GC = ["G", "C"]
AT = ["A", "T"]

def get_filestamp(suffix=None):
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if suffix:
        timestr += f"-{suffix}"
    return timestr


def write_primers(primers, total_time_sec=-1, strategy=""):
    # Write primer
    os.makedirs(DATA_FOLDER, exist_ok=True)
    with open(
        os.path.join(
            DATA_FOLDER, get_filestamp(suffix=f"{strategy}-{len(primers)}primers.txt")
        ),
        "w",
    ) as f:
        f.write(f"Total time (seconds):\t{total_time_sec}\n")
        f.write("\n".join(primers))


def random_primer(length=20):
    return "".join(random.choices(ALPHABET, k=length))


def random_primer_with_balanced_gc(length=20, min_gc=9, max_gc=11):
    """
    Generate primer of length with at least min_g GCs and at most max_gc GCs.
    """
    # Generate string of the right length with only As and Ts
    primer = random.choices(AT, k=length)
    # Choose how many Gs or Cs to have
    num_gc = random.randint(min_gc, max_gc)
    # For each of those nts, choose if it's a G or C
    gc_str = random.choices(GC, k=num_gc)
    # Pick the indices we should swap to G or C
    indices = range(0, length)
    gc_indices = random.sample(indices, num_gc)
    # Swap all necessary indices to G/C
    for idx, g_or_c in enumerate(gc_str):
        idx_to_replace = gc_indices[idx]
        primer[idx_to_replace] = g_or_c
    return "".join(primer)
