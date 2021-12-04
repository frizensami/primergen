#!/usr/bin/env python3

"""
1. Generate some N = large number of valid GC-content strings, where N > target (probably)
2. Compute pair-wise edit distance between all pairs, generate edge between two nodes if edit distance < 8
3. Pick a random node in the graph and remove it
4. Delete all of its neighbors
5. Repeat 3 and 4 until there are no more nodes
"""

import networkx as nx
from check import *
from strategy_generic import BasePrimerGenerator
from util import random_primer_with_balanced_gc
import random
import math

PRINT_EVERY_NTH_ITERATION_N = 20000


class CliquePrimerGenerator(BasePrimerGenerator):
    def __init__(self, target=TARGET_PRIMERS):
        super().__init__(target=target, strategy="delob")
        self.init_primers = []

    def generate(self):
        # Generating large number N of initial primers
        N = int(RANDOM_PRIMERS_TO_GENERATE)
        primers = self.generate_n_primers(n=N)
        self.num_starting_primers = N  # For writing to file later

        # How many
        combinations = math.comb(N, 2)

        # Initialize graph with n primers
        g = nx.Graph()
        edges = []
        weights = []

        # Compute weights between all edges
        # Get all pairs of primers + their indices - this returns [((index of item, pair1_item1), (index of item, pair1_item2)), ...]
        pairs_with_idx = itertools.combinations(enumerate(primers), 2)
        for ((idx1, primer1), (idx2, primer2)) in pairs_with_idx:
            self.new_iteration()
            # Don't let printing slow us down
            if self.iterations % PRINT_EVERY_NTH_ITERATION_N == 0:
                print(
                    f"Iter {self.iterations}\tComputing edit distance between primer {idx1} and {idx2}\t({round(self.iterations / combinations * 100, 2)}%)"
                )
            # Compute edit distance
            dist = get_edit_distance(primer1, primer2, limit=MIN_EDIT_DISTANCE)
            # DeLOB: only if the nodes are too close, add them to the graph
            if dist < MIN_EDIT_DISTANCE:
                # Don't let printing slow us down
                if self.iterations % PRINT_EVERY_NTH_ITERATION_N == 0:
                    print(
                        f"Iter {self.iterations}\tAdding edge between {idx1} and {idx2}, distance is {dist}, less than {MIN_EDIT_DISTANCE}"
                    )
                # g.add_edge(idx1, idx2, weight=dist)
                edges.append((idx1, idx2))
                weights.append(dist)

        # Add all the edges we computed as tuples of (node1, node2)
        print(f"Adding {len(edges)} edges...")
        g.add_edges_from(edges)
        print(f"Done adding edges")

        """
        DeLOB network algorithm
        1. Pick a random node in the graph
        2. Take note of all its neighbours
        3. Remove the random node from the graph and insert it into the list of valid primers
        4. Remove all its neighbors from the graph
        5. Repeat steps 1 -- 4 until there are no more edges
        """

        number_of_edges = g.number_of_edges()
        while number_of_edges != 0:
            # 1. Pick random node
            new_valid_primer = random.choice(list(g.nodes()))
            print(f"New primer: {new_valid_primer}")
            self.found_new_primer(primers[new_valid_primer])
            # 2. Get neighbors
            neighbors = g.neighbors(new_valid_primer)
            # 3 & 4: Remove node and neighbors from graph
            g.remove_node(new_valid_primer)
            g.remove_nodes_from(neighbors)
            # 5. Recompute current state, print, and cycle back
            number_of_edges = g.number_of_edges()
            number_of_nodes = g.number_of_nodes()
            print(
                f"Edges remaining: {number_of_edges}, Nodes remaining: {number_of_nodes}"
            )
            self.print_metrics()

    def generate_n_primers(self, n):
        print(f"Generating {n} primers to begin...")
        primers = []
        for i in range(n):
            primers.append(random_primer_with_balanced_gc())
        print(f"Done generating {n} primers")
        return primers


if __name__ == "__main__":
    CliquePrimerGenerator().execute()
