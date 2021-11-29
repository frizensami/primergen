#!/usr/bin/env python3

"""
1. Generate some N = large number of valid GC-content strings, where N > target (probably)
2. Compute pair-wise edit distance between all pairs, generate edge between two nodes if edit distance >= 8
3. Find the largest clique in the graph and see its size
4. Repeat 1 - 3 until we're happy? Or from 3, do a random strategy to try to add more nodes?
"""

import networkx as nx
from check import *
from strategy_generic import BasePrimerGenerator
from util import random_primer_with_balanced_gc

class CliqueMaxWeightPrimerGenerator(BasePrimerGenerator):
    def __init__(self, target=TARGET_PRIMERS):
        super().__init__(target=target, strategy="max-weight-clique")
        self.init_primers = []

    def generate(self):
        # Generating large number N of initial primers
        N = int(self.target * 0.1)
        primers = self.generate_n_primers(n=N)

        # Initialize graph with n primers
        g = nx.Graph()

        edges = []
        # Compute weights between all edges
        # Get all pairs of primers + their indices - this returns [((index of item, pair1_item1), (index of item, pair1_item2)), ...]
        pairs_with_idx = itertools.combinations(enumerate(primers), 2)
        for ((idx1, primer1), (idx2, primer2)) in pairs_with_idx:
            print(f"Computing edit distance between primer {idx1} and {idx2}")
            # Compute edit distance
            dist = get_edit_distance(primer1, primer2)
            # If these are a valid pair, add them as an edge
            if dist >= MIN_EDIT_DISTANCE:
                print(f"Adding edge between {idx1} and {idx2}, distance is {dist}")
                #g.add_edge(idx1, idx2, weight=dist)
                edges.append((idx1, idx2, dist))

        print(f"Adding {len(edges)} edges...")
        g.add_weighted_edges_from(edges)
        print(f"Done adding edges")

        print(f"Computing node weights...")
        for node in g.nodes:
            edges_with_weights = g.edges(node, data='weight')
            weight_sum = sum(map(lambda x: x[2], edges_with_weights))
            g.nodes[node]['weight'] = weight_sum
        print(f"Done computing node weights")

        print(f"Computing largest and max_weight clique...")
        largest_clique, weight = nx.algorithms.clique.max_weight_clique(g)
        print(f"Done computing largest clique with total weight {weight}...")
        print(largest_clique)

        final_primers = [primers[idx] for idx in largest_clique]
        self.found_new_primers(final_primers)
        print(f"Number of primers found: {len(final_primers)}")




    def generate_n_primers(self, n):
        print(f"Generating {n} primers to begin...")
        primers = []
        for i in range(n):
            primers.append(random_primer_with_balanced_gc())
        print(f"Done generating {n} primers")
        return primers



if __name__ == "__main__":
    CliqueMaxWeightPrimerGenerator().execute()
