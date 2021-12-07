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


class DelobMinDegreeNeighborsPrimerExtractor(BasePrimerExtractor):
    def __init__(
        self,
        initial_primers=None,
        target=TARGET_PRIMERS,
        strategy="delob-mindegree-neighbors",
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
                            f"Iter {self.iterations}\tAdding edge between {idx1} and {idx2}, distance is {dist}, less than {MIN_EDIT_DISTANCE}"
                        )
                    # g.add_edge(idx1, idx2, weight=dist)
                    edges.append((idx1, idx2))

        # Add all the edges we computed as tuples of (node1, node2)
        print(f"Adding {len(edges)} edges...")
        g.add_edges_from(edges)
        del edges
        print(f"Done adding edges")

        # Add primers that weren't added to the graph to the final list (no conflicts)
        self.found_new_primers(get_no_conflict_primers(g, self.initial_primers))

        # with open("edges.txt", "w") as f:
        #     f.write(str(edges))

        """
        DeLOB-mindegree-neighbors network algorithm
        1. Find all the nodes in the graph with /minimum degree/
        2. Compute the minimum degree of all the neighbors of neighbours of each such node
        3. Pick the node N with the smallest average degree of neighbours(neighbours(N)) -- removal of this node will "free" as many nodes as possible
        4. Remove the node N from the graph and insert it into the list of valid primers
        5. Delete all its neighbors from the graph (min degree --> fewer neighbors deleted)
        6. Repeat steps 1 -- 5 until there are no more edges
        """

        number_of_edges = g.number_of_edges()
        initial_number_of_nodes = g.number_of_nodes()
        # Start a timer for when we start to compute edit dists
        start_node_remove_time = time.process_time()
        while number_of_edges != 0:
            # MINDEGREE-NEIGHBORS STARTS HERE
            sorted_nodes_by_degree = sorted(g.degree(), key=lambda x: x[1])
            ## FIND THE NODE WITH THE LOWEST AVERAGE SECOND-LEVEL-NEIGHBOR DEGREE
            # The minimum degree - ignore nodes if degree is higher than this
            min_degree_in_graph = sorted_nodes_by_degree[0][1]
            best_node = None
            lowest_nbr_degree = float("inf")
            for node, degree in sorted_nodes_by_degree:
                # print(
                #     f"Processing node {node} with degree {degree} (min is: {min_degree_in_graph})"
                # )
                if degree > min_degree_in_graph:
                    break
                # For each valid node, we need to find the average degree of all of its neighbors of neighbors
                snd_nbr_avg_degree = self.avg_degree_second_neighbors(g, node)
                if snd_nbr_avg_degree < lowest_nbr_degree:
                    lowest_nbr_degree = snd_nbr_avg_degree
                    best_node = node
                # print(
                #     f"Processed node {node} with degree {degree} (min is: {min_degree_in_graph}), second neighbor degree: {snd_nbr_avg_degree}, lowest so far is {lowest_nbr_degree} (node: {best_node})\n"
                # )

            # MINDEGREE-NEIGHBORS ENDS HERE
            # 1. Pick random node
            new_valid_primer = best_node
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

    def avg_degree_second_neighbors(self, graph, node):
        """
        Compute the average of the degrees of 2nd level neighbors from a given node
        """
        second_nbrs = []
        second_nbr_degrees = []
        # Find all neighbors of current node
        for n in graph.neighbors(node):
            # TODO: this doesn't account for repeated neighbors, but not sure if it's worth it
            # Add the degree of all their neighbors to a list
            second_nbrs.append(n)
            second_nbr_degrees.append(graph.degree(n))
        # print(
        #     f"Second neighbor degrees: {second_nbr_degrees} (seconds nbrs: {second_nbrs})"
        # )
        if len(second_nbr_degrees) == 0:
            # No second neigbors - return 0 to avoid divide by 0 error
            return 0
        else:
            # Average of second neighbor degrees
            return sum(second_nbr_degrees) / len(second_nbr_degrees)


if __name__ == "__main__":
    DelobMinDegreeNeighborsPrimerExtractor().execute()
