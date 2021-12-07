#!/usr/bin/env python3
import random
import math
import time
import sys
import re

import networkx as nx
from primergen.common.check import *
from primergen.common.util import random_primer_with_balanced_gc
from primergen.common.graph_utils import *

from .base import BasePrimerExtractor

PRINT_EVERY_NTH_ITERATION_N = 20000


class GreedyPrimerExtractor(BasePrimerExtractor):
    def __init__(self, initial_primers=None, target=TARGET_PRIMERS, strategy="greedy"):
        super().__init__(initial_primers, target, strategy)

    def generate(self):
        """
        Just run through the list of initial primers and grab ones that don't clash with the current list
        """
        start_compute_time = time.process_time()
        for init_primer in self.initial_primers:
            self.iterations += 1
            # Don't let printing slow us down
            if self.iterations % PRINT_EVERY_NTH_ITERATION_N == 0:
                elapsed = time.process_time() - start_compute_time
                progress_fraction = self.iterations / len(self.initial_primers)
                eta_min = int(
                    (elapsed / progress_fraction) * (1 - progress_fraction) / 60
                )
                eta_sec = int(
                    (elapsed / progress_fraction) * (1 - progress_fraction) % 60
                )
                progress_percent = round(progress_fraction * 100, 2)
                print(
                    f"Iter {self.iterations}\t({progress_percent}%)\tETA: {eta_min} m {eta_sec} s"
                )

            for other_primer in self.primers:
                if not is_primer_pair_valid(
                    init_primer, other_primer, limit=MIN_EDIT_DISTANCE
                ):
                    super().new_edit_error()
                    break
            else:
                super().found_new_primer(init_primer)


if __name__ == "__main__":
    GreedyPrimerExtractor().execute()
