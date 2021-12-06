#!/usr/bin/env python3
from strategy_generic import BasePrimerGenerator
from check import *
import random
import numpy

# OG Seed
RANDOM_SEED = 246
NP_RANDOM_SEED = 4812


# RANDOM_SEED = 111
# NP_RANDOM_SEED = 7777


class BasePrimerExtractor(BasePrimerGenerator):
    """
    Is passed a list of GC-valid primers (pre-generated), goal is to find the maximum size valid primer library within it.
    """

    def __init__(self, initial_primers, target=TARGET_PRIMERS, strategy="base"):
        super().__init__(target=target, strategy=strategy)
        self.initial_primers = initial_primers
        self.num_starting_primers = len(initial_primers)

        # Seed random same for all extractors
        random.seed(RANDOM_SEED)
        numpy.random.seed(NP_RANDOM_SEED)

    def generate(self):
        """
        NOTE: Inheriting classes will implement this
        """
        raise NotImplementedError
