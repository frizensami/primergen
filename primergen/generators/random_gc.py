#!/usr/bin/env python3#!
import random
import math
import time
import sys
import re

import networkx as nx
from primergen.common.check import *
from primergen.common.util import (
    random_primer_with_balanced_gc,
)
from primergen.common.graph_utils import *

from .base import BasePrimerGenerator


class RandomBalancedGCPrimerGenerator(BasePrimerGenerator):
    def __init__(self, target=TARGET_PRIMERS):
        super().__init__(target=target, strategy="random-balanced-gc")

    def generate(self):
        while len(self.primers) < self.target:
            # Stats
            super().new_iteration()
            self.iterations += 1
            super().print_metrics()

            # Create new primer
            primer = random_primer_with_balanced_gc()
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
    RandomBalancedGCPrimerGenerator().execute()
