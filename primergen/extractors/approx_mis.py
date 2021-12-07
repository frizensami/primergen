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

sys.setrecursionlimit(150000)


class ApproxMisPrimerExtractor(BasePrimerExtractor):
    def __init__(
        self, initial_primers=None, target=TARGET_PRIMERS, strategy="approx-clique"
    ):
        super().__init__(initial_primers, target, strategy)

    def generate(self):
        # How many primer combos total (n choose 2)
        combinations = math.comb(self.num_starting_primers, 2)

        # Initialize graph with n primers
        g = nx.Graph()
        edges = []

        if USE_EDGES_FILE:
            edges = self.get_edges_file_pkl()
        else:

            # Compute weights between all edges
            # Get all pairs of primers + their indices - this returns [((index of item, pair1_item1), (index of item, pair1_item2)), ...]
            pairs_with_idx = itertools.combinations(enumerate(self.initial_primers), 2)
            # Start a timer for when we start to compute edit dists
            start_edit_dist_compute_time = time.process_time()
            for ((idx1, primer1), (idx2, primer2)) in pairs_with_idx:
                self.iterations += 1
                # Don't let printing slow us down
                if self.iterations % PRINT_EVERY_NTH_ITERATION_N == 0:
                    elapsed = time.process_time() - start_edit_dist_compute_time
                    progress_fraction = (self.iterations) / combinations
                    eta_min = int(
                        (elapsed / progress_fraction) * (1 - progress_fraction) / 60
                    )
                    eta_sec = int(
                        (elapsed / progress_fraction) * (1 - progress_fraction) % 60
                    )
                    progress_percent = round(progress_fraction * 100, 2)
                    print(
                        f"Iter {self.iterations}\tComputing edit distance between primer {idx1} and {idx2}\t({progress_percent}%)\tETA: {eta_min} m {eta_sec} s"
                    )
                # Compute edit distance
                dist = get_edit_distance_with_limit(
                    primer1, primer2, limit=MIN_EDIT_DISTANCE
                )
                # DeLOB: only if the nodes are too close, add them to the graph
                if dist < MIN_EDIT_DISTANCE:
                    # Don't let printing slow us down
                    if self.iterations % PRINT_EVERY_NTH_ITERATION_N == 0:
                        print(
                            f"Iter {self.iterations}\tAdding edge between {idx1} and {idx2}, distance is {dist} >= {MIN_EDIT_DISTANCE}"
                        )
                    # g.add_edge(idx1, idx2, weight=dist)
                    edges.append((idx1, idx2))

        # Add all the edges we computed as tuples of (node1, node2)
        print(f"Adding {len(edges)} edges...")
        g.add_edges_from(edges)
        print(f"Done adding edges")

        # Add primers that weren't added to the graph to the final list (no conflicts)
        self.found_new_primers(get_no_conflict_primers(g, self.initial_primers))

        """
        Approximation algorithm to find maximum independent set (built-in from networkx)
        """

        print(f"Computing APPROXIMATE (V/(logV)^2) largest maximum independent set...")
        largest_clique = nx.algorithms.approximation.clique.maximum_independent_set(g)
        print(f"Done computing largest cliques...")
        print(largest_clique)

        final_primers = [self.initial_primers[idx] for idx in largest_clique]
        self.found_new_primers(final_primers)
        print(f"Number of primers found: {len(final_primers)}")


if __name__ == "__main__":
    ApproxMisPrimerExtractor().execute()
