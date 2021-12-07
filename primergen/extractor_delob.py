#!/usr/bin/env python3
from extractor_generic import BasePrimerExtractor
import networkx as nx
from check import *
from util import random_primer_with_balanced_gc
import random
import math
import time
from extractor import get_primers, PRIMER_FILE
import re
from graph_utils import *

import pickle

PRINT_EVERY_NTH_ITERATION_N = 20000


class DelobPrimerExtractor(BasePrimerExtractor):
    def __init__(self, initial_primers, target=TARGET_PRIMERS, strategy="delob"):
        super().__init__(initial_primers, target, strategy)

    def generate(self):
        # How many primer combos total (n choose 2)
        combinations = math.comb(self.num_starting_primers, 2)

        # Initialize graph with n primers
        g = nx.Graph()
        edges = []

        if USE_EDGES_FILE:
            edges_filename = PRIMER_FILE.split(".")[0] + "-edges.txt"
            with open(edges_filename, "r") as f:
                edges_str = f.read()
                edges_ints = list(
                    map(int, re.sub("[ [\\]\(\)]", "", edges_str).split(","))
                )
                edges = [
                    (edges_ints[i], edges_ints[i + 1])
                    for i in range(0, len(edges_ints), 2)
                ]
                # print(edges)
                del edges_str
                del edges_ints
                # print(edges)
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
                            f"Iter {self.iterations}\tAdding edge between {idx1} and {idx2}, distance is {dist}, less than {MIN_EDIT_DISTANCE}"
                        )
                    # g.add_edge(idx1, idx2, weight=dist)
                    edges.append((idx1, idx2))

            # # Write to edges
            # with open("edges-9ed.pkl", "wb") as f:
            #     pickle.dump(edges, f)

        # Add all the edges we computed as tuples of (node1, node2)
        print(f"Adding {len(edges)} edges...")
        g.add_edges_from(edges)
        print(f"Done adding edges")

        # Add primers that weren't added to the graph to the final list (no conflicts)
        self.found_new_primers(get_no_conflict_primers(g, self.initial_primers))

        """
        DeLOB network algorithm
        1. Pick a random node in the graph
        2. Take note of all its neighbours
        3. Remove the random node from the graph and insert it into the list of valid primers
        4. Remove all its neighbors from the graph
        5. Repeat steps 1 -- 4 until there are no more edges
        """

        number_of_edges = g.number_of_edges()
        initial_number_of_nodes = g.number_of_nodes()
        # Start a timer for when we start to compute edit dists
        start_node_remove_time = time.process_time()
        while number_of_edges != 0:
            # 1. Pick random node
            new_valid_primer = random.choice(list(g.nodes()))
            print(f"New primer: {new_valid_primer}")
            self.found_new_primer(self.initial_primers[new_valid_primer])
            # 2. Get neighbors
            neighbors = g.neighbors(new_valid_primer)
            # 3 & 4: Remove node and neighbors from graph
            g.remove_node(new_valid_primer)
            g.remove_nodes_from(neighbors)
            # 5. Recompute current state, print, and cycle back
            number_of_edges = g.number_of_edges()
            number_of_nodes = g.number_of_nodes()
            # Stats
            elapsed = time.process_time() - start_node_remove_time
            progress_fraction = (
                initial_number_of_nodes - number_of_nodes
            ) / initial_number_of_nodes
            eta_min = int((elapsed / progress_fraction) * (1 - progress_fraction) / 60)
            eta_sec = int((elapsed / progress_fraction) * (1 - progress_fraction) % 60)
            progress_percent = round(progress_fraction * 100, 2)
            print(
                f"Edges remaining: {number_of_edges}, Nodes remaining: {number_of_nodes}\t({progress_percent}%)\tETA: {eta_min} m {eta_sec} s"
            )
            self.print_metrics()


if __name__ == "__main__":
    DelobPrimerExtractor(initial_primers=get_primers()).execute()
