#!/usr/bin/env python3


def get_no_conflict_primers(g, initial_primers):
    # We need a step where all primers that aren't added to the graph are added as valid primers
    # Since they don't have even 1 edge, they confict with no other primers in the set at all
    """
    1. Convert the set of node IDs in the graph into a set
    2. Take R = (set of all possible node IDs) - (set of node IDs in graph)
    3. Add each primer in R into the valid primer set
    """
    print("Finding primers with 0 conflicts with rest of set.")
    all_nodes_in_graph = set(g.nodes)
    all_possible_nodes = set(range(0, len(initial_primers)))
    valid_nodes_with_no_conflicts = list(
        map(
            lambda idx: initial_primers[idx],
            all_possible_nodes - all_nodes_in_graph,
        )
    )
    print(f"{len(valid_nodes_with_no_conflicts)} no-conflict primers found.")
    return valid_nodes_with_no_conflicts
