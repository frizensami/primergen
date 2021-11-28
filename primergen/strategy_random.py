#!/usr/bin/env python3
from Bio.Seq import Seq, MutableSeq
import numpy as np
from check import *
from util import write_primers, random_primer
import time
import random


def generate_random(target=2000):
    """
    This is a baseline primer generation strategy based on pure randomness.
    1. Generate random 20-mer.
    2. See if it satisfies the GC content requirements. If not, regenerate. (TODO: just generate it with the GC already)
    3. If satisfies, chcek its edit distance against all other primers in the library. If it's OK, add it.
    4. Continue until the target number of primers is hit.

    target is the number of primers we expect at the end
    """
    primers = []
    iters = 0
    gc_invalid = 0
    edit_invalid = 0

    try:
        while len(primers) < 2000:
            # Stats
            iters += 1
            print(
                f"Iteration: {iters}\tPrimers:{len(primers)}\tGC invalid: {gc_invalid}\tEdit invalid\t{edit_invalid}"
            )
            # Create new primer
            primer = random_primer()
            print(primer)
            # Check for GC content
            if not is_gc_valid(primer):
                gc_invalid += 1
                continue
            # Check for edit distance
            for other_primer in primers:
                if not is_primer_pair_valid(primer, other_primer):
                    edit_invalid += 1
                    break
            else:
                # We didn't break out of previous loop, all good, add to primer library
                primers.append(primer)
    except KeyboardInterrupt:
        print(f"Exited early with {len(primers)} primers!")
    return primers


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
