#!/usr/bin/env python3

from extractor_generic import BasePrimerExtractor
from check import *
from extractor import get_primers


class GreedyPrimerExtractor(BasePrimerExtractor):
    def __init__(self, initial_primers, target=TARGET_PRIMERS, strategy="greedy"):
        super().__init__(initial_primers, target, strategy)

    def generate(self):
        """
        Just run through the list of initial primers and grab ones that don't clash with the current list
        """
        for init_primer in self.initial_primers:
            for other_primer in self.primers:
                if not is_primer_pair_valid(
                    init_primer, other_primer, limit=MIN_EDIT_DISTANCE
                ):
                    super().new_edit_error()
                    break
            else:
                super().found_new_primer(init_primer)


if __name__ == "__main__":
    GreedyPrimerExtractor(initial_primers=get_primers()).execute()
