#!/usr/bin/env python3
import numpy as np
from check import *
from util import (
    write_primers,
    generate_primer_from_frequencies_and_balanced_gc,
    generate_primer_from_frequencies_and_balanced_gc_by_rerolling,
)
import time
import random
from strategy_generic import BasePrimerGenerator

"""
TODO: New strings are generated based on the inverse of the current frequencies for each position. Fix GC after.
The hope is that we get a string that is as different as possible from the rest.

E.g., if we current have primers ATGC.., TGCA, GCAT.., the next primer will probably be CATG... (least popular choices for all those positions)

I guess this is inspired by Cantor's diagonalization argument where we generate new numbers by putting a new value in each position.
Well, not really because we don't have an infinite number of nts.

TODO: make new strategy where the GC-imbalance is corrected deterministically to prevent re-rolls.
TODO: make new strategy where sometimes we randomize completely
TODO: store primers found over time / iteraitons as a graph
TODO: store all the primer to primer distances in a file so we can plot a force directed graph of how far apart everything is
"""


class RandomBalancedGCFrequenciesPrimerGenerator(BasePrimerGenerator):
    def __init__(self, target=TARGET_PRIMERS):
        super().__init__(target=target, strategy="random-balanced-gc-frequencies")
        # Frequencies for each position in the primer
        self.counts = [{"A": 1, "T": 1, "G": 1, "C": 1} for i in range(PRIMER_LENGTH)]

    def generate(self):
        while len(self.primers) < self.target:
            # Stats
            super().new_iteration()
            super().print_metrics()

            # Get an inverse of the current frequencies so we can generate an "as different as possible" primer
            inv_frequencies = self.counts_to_inverse_frequencies()
            # print(self.counts)
            # print(inv_frequencies)

            # Create new primer
            primer = generate_primer_from_frequencies_and_balanced_gc_by_rerolling(
                frequencies=inv_frequencies
            )
            print(primer)

            # Check for GC content
            if not is_gc_valid(primer):
                super().new_gc_error()
                continue
            # Check for edit distance
            for other_primer in self.primers:
                if not is_primer_pair_valid(
                    primer, other_primer, limit=MIN_EDIT_DISTANCE
                ):
                    super().new_edit_error()
                    break
            else:
                super().found_new_primer(primer)
                # Update frequency counts for each position from the new primer
                self.update_counts_from_new_primer(primer)

    def update_counts_from_new_primer(self, primer):
        # Increase the count for the nucleotide in each position as necessary
        for idx, nt in enumerate(primer):
            # print(f"Updating counts idx {idx}, nt {nt} by 1")
            self.counts[idx][nt] = self.counts[idx][nt] + 1
        # print(f"Finished updating, final result {self.counts}")

    def counts_to_inverse_frequencies(self):
        """
        Given the raw frequency counts, generate weights for random.choice that prefers the /least/ frequent elements.
        Imagine if one position's frequency count is {A: 2, T: 2, G: 7, C: 5}.
        We can just inverse each frequency to get the relative weights to "restore" equilibrium.
        For e.g., the result should be {A: 1/2, T: 1/2, G: 1/7, C: 1/5}
        """
        inverse_frequencies = []
        for count in self.counts:
            inverse_count = [
                1 / count["A"],
                1 / count["T"],
                1 / count["G"],
                1 / count["C"],
            ]
            inverse_frequencies.append(inverse_count)
        return inverse_frequencies


if __name__ == "__main__":
    RandomBalancedGCFrequenciesPrimerGenerator().execute()
