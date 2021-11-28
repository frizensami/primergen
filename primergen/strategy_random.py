#!/usr/bin/env python3
from Bio.Seq import Seq, MutableSeq
import numpy as np
from check import *
from util import write_primers
import time


def generate_random(target=2000):
    """
    This is a baseline primer generation strategy based on pure randomness.
    1. Generate random 20-mer.
    2. See if it satisfies the GC content requirements. If not, regenerate. (TODO: just generate it with the GC already)
    3. If satisfies, chcek its edit distance against all other primers in the library. If it's OK, add it.
    4. Continue until the target number of primers is hit.

    target is the number of primers we expect at the end
    """
    return []
    pass


if __name__ == "__main__":
    start_time = time.process_time()
    # Actually execute
    primers = generate_random()
    # Finished executing
    end_time = time.process_time()
    total_time = end_time - start_time
    print(primers)
    print(f"Total time (sec): {total_time}")
    write_primers(primers, total_time_sec=total_time, strategy="random")
