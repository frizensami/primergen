#!/usr/bin/env python3
import numpy as np
from check import *
from util import write_primers, random_primer_with_weights_and_balanced_gc
import time
import random
from strategy_generic import BasePrimerGenerator

"""
TODO: New strings are generated based on the inverse of the current frequencies for each position. Fix GC after.
The hope is that we get a string that is as different as possible from the rest.
"""

class RandomBalancedGCFrequenciesPrimerGenerator(BasePrimerGenerator):
    def __init__(self, target=TARGET_PRIMERS):
        super().__init__(target=target, strategy="random-balanced-gc")
        self.counts = [{'A': 0, 'T': 0, 'G': 0, 'C': 0}] * PRIMER_LENGTH

    def generate(self):
        while len(self.primers) < self.target:
            # Stats
            super().new_iteration()
            super().print_metrics()

            # Create new primer
            primer = random_primer_with_weights_and_balanced_gc()
            print(primer)

            # Update frequency counts for each position from the new primer
            self.update_counts_from_new_primer(primer)

            # Check for GC content
            if not is_gc_valid(primer):
                super().new_gc_error()
                continue
            # Check for edit distance
            for other_primer in self.primers:
                if not is_primer_pair_valid(primer, other_primer, limit=(MIN_EDIT_DISTANCE + 1)):
                    super().new_edit_error()
                    break
            else:
                super().found_new_primer(primer)

    def update_counts_from_new_primer(primer):
        # Increase the count for the nucleotide in each position as necessary
        for idx, nt in enumerate(primer):
            self.counts[idx][nt] = self.counts[idx][nt] + 1


if __name__ == "__main__":
    RandomBalancedGCFrequenciesPrimerGenerator().execute()
