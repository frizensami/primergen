#!/usr/bin/env python3
from check import *
from util import write_primers, random_primer
import time
import random
from strategy_generic import BasePrimerGenerator


class RandomPrimerGenerator(BasePrimerGenerator):
    def __init__(self, target=TARGET_PRIMERS):
        super().__init__(target=target, strategy="random")

    def generate(self):
        while len(self.primers) < self.target:
            # Stats
            self.iterations += 1
            super().print_metrics()

            # Create new primer
            primer = random_primer()
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


if __name__ == "__main__":
    RandomPrimerGenerator().execute()
