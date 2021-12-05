#!/usr/bin/env python3
from strategy_generic import BasePrimerGenerator
from check import *


class BasePrimerExtractor(BasePrimerGenerator):
    """
    Is passed a list of GC-valid primers (pre-generated), goal is to find the maximum size valid primer library within it.
    """

    def __init__(self, initial_primers, target=TARGET_PRIMERS, strategy="base"):
        super().__init__(target=target, strategy=strategy)
        self.initial_primers = initial_primers
        self.num_starting_primers = len(initial_primers)

    def generate(self):
        """
        NOTE: Inheriting classes will implement this
        """
        raise NotImplementedError
