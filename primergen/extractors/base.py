#!/usr/bin/env python3
import random
import numpy
import re
import pickle

from importlib import resources
from primergen.generators.generic import BasePrimerGenerator
from primergen.common.check import *


# OG Seed
RANDOM_SEED = 246
NP_RANDOM_SEED = 4812

# Alt seed
# RANDOM_SEED = 111
# NP_RANDOM_SEED = 7777

# PRIMER_FILE = "20211205-141133-100000-primers.txt"
# PRIMER_FILE = "20211204-141808-10000-primers.txt"
# PRIMER_FILE = "20211204-174137-5000-primers.txt"
# PRIMER_FILE = "20211205-182302-2000-primers.txt"
# PRIMER_FILE = "20211206-100844-1000-primers.txt"
# PRIMER_FILE = "20211206-101224-500-primers.txt"
# PRIMER_FILE = "20211206-102017-200-primers.txt"
PRIMER_FILE = "20211206-093547-100-primers.txt"


class BasePrimerExtractor(BasePrimerGenerator):
    """
    Is passed a list of GC-valid primers (pre-generated), goal is to find the maximum size valid primer library within it.
    """

    def __init__(self, initial_primers=None, target=TARGET_PRIMERS, strategy="base"):
        super().__init__(target=target, strategy=strategy)
        if initial_primers:
            self.initial_primers = initial_primers
        else:
            self.initial_primers = self.get_primers()
        self.num_starting_primers = len(self.initial_primers)

        # Seed random same for all extractors
        random.seed(RANDOM_SEED)
        numpy.random.seed(NP_RANDOM_SEED)

    def generate(self):
        """
        NOTE: Inheriting classes will implement this
        """
        raise NotImplementedError

    def get_primers(self):
        initial_primers = []
        with resources.open_text("primergen.input", PRIMER_FILE) as f:
            # with open(PRIMER_FILE, "r") as f:
            initial_primers = f.read().splitlines()
            return initial_primers

    def get_edges_file_txt(self):
        edges_filename = PRIMER_FILE.split(".")[0] + "-edges.txt"
        with resources.open_text("primergen.input", edges_filename) as f:
            edges_str = f.read()
            edges_ints = list(map(int, re.sub("[ [\\]\(\)]", "", edges_str).split(",")))
            edges = [
                (edges_ints[i], edges_ints[i + 1]) for i in range(0, len(edges_ints), 2)
            ]
            # print(edges)
            del edges_str
            del edges_ints
            # print(edges)
            return edges

    def get_edges_file_pkl(self):
        edges_filename = PRIMER_FILE.split(".")[0] + "-edges.pkl"
        with resources.open_binary("primergen.input", edges_filename) as f:
            return pickle.load(f)
