#!/usr/bin/env python3

import time
import os
import random

DATA_FOLDER = "data/"
ALPHABET = ["A", "T", "C", "G"]
GC = ["G", "C"]
AT = ["A", "T"]


def get_filestamp(suffix=None):
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if suffix:
        timestr += f"-{suffix}"
    return timestr


def write_primers(primers, primer_found_times, total_time_sec=-1, strategy=""):
    # Write primer
    os.makedirs(DATA_FOLDER, exist_ok=True)
    path_without_ext = os.path.join(
        DATA_FOLDER, get_filestamp(suffix=f"{strategy}-{len(primers)}primers")
    )
    with open(
        path_without_ext + ".txt",
        "w",
    ) as f:
        f.write(f"Total time (seconds):\t{total_time_sec}\n")
        f.write("\n".join(primers))

    with open(
        path_without_ext + "-primertimes.txt",
        "w",
    ) as f:
        for (time, primers) in primer_found_times:
            f.write(f"{time},{primers}\n")


def random_primer(length=20):
    return "".join(random.choices(ALPHABET, k=length))


def random_primer_with_balanced_gc(length=20, min_gc=9, max_gc=11):
    """
    Generate primer of length with at least min_g GCs and at most max_gc GCs.
    """
    # Generate string of the right length with only As and Ts
    primer = random.choices(AT, k=length)
    # Choose how many Gs or Cs to have
    num_gc = random.randint(min_gc, max_gc)
    # For each of those nts, choose if it's a G or C
    gc_str = random.choices(GC, k=num_gc)
    # Pick the indices we should swap to G or C
    indices = range(0, length)
    gc_indices = random.sample(indices, num_gc)
    # Swap all necessary indices to G/C
    for idx, g_or_c in enumerate(gc_str):
        idx_to_replace = gc_indices[idx]
        primer[idx_to_replace] = g_or_c
    return "".join(primer)


def generate_primer_from_frequencies_and_balanced_gc(
    frequencies, length=20, min_gc=9, max_gc=11
):
    # Generate /length/ number of nts with relative frequencies
    primer = [random.choices(ALPHABET, frequencies[i], k=1)[0] for i in range(length)]
    # This primer may not have a balanced GC count, so count GCs
    num_gcs = sum(map(lambda nt: 1 if nt == "G" or nt == "C" else 0, primer))

    # print(f"Initial primer is {primer}")
    # print(f"We have {num_gcs} G and Cs")

    if num_gcs < min_gc:
        # Do we need MORE GCs?
        # How many more GCs we have to add
        diff = min_gc - num_gcs
        # Where we could add GCs them (where there's only As and Ts)
        non_gc_indices = [
            idx for idx, elem in enumerate(primer) if elem == "A" or elem == "T"
        ]
        # Choose indices where we should insert G/C from this list
        indices_to_insert_at = random.sample(non_gc_indices, diff)
        # Insert the G/Cs at those positions
        for idx in indices_to_insert_at:
            primer[idx] = random.choice(GC)
        # print(f"Need {diff} more GCs")
        # print(f"Non GCs are at indices {non_gc_indices}")
        # print(f"We will insert GCs at {indices_to_insert_at}")

    elif num_gcs > max_gc:
        # Or maybe we have too many GCs
        diff = num_gcs - max_gc
        # Where we could add ATs (where there's only Gs and Cs)
        gc_indices = [
            idx for idx, elem in enumerate(primer) if elem == "G" or elem == "C"
        ]
        # Choose indices where we should insert G/C from this list
        indices_to_insert_at = random.sample(gc_indices, diff)
        # Insert the A/Ts at those positions
        for idx in indices_to_insert_at:
            primer[idx] = random.choice(AT)
        # print(f"Need {diff} fewer GCs")
        # print(f"GCs are at indices {gc_indices}")
        # print(f"We will insert ATs at {indices_to_insert_at}")

    # print(f"Final primer: {primer}")
    return "".join(primer)


def generate_primer_from_frequencies_and_balanced_gc_by_rerolling(
    frequencies, length=20, min_gc=9, max_gc=11
):
    """
    Generate primers based on frequency metrics but re-rolls if we don't get the right GC proportion
    """
    while True:
        # Generate /length/ number of nts with relative frequencies
        primer = [
            random.choices(ALPHABET, frequencies[i], k=1)[0] for i in range(length)
        ]
        # This primer may not have a balanced GC count, so count GCs
        num_gcs = sum(map(lambda nt: 1 if nt == "G" or nt == "C" else 0, primer))

        # print(f"Initial primer is {primer}")
        # print(f"We have {num_gcs} G and Cs")

        if num_gcs >= min_gc and num_gcs <= max_gc:
            # print(f"Got a balanced GC primer {primer}, returning")
            return "".join(primer)
