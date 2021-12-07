# PrimerGen: A primer generation / extraction testbed with benchmarks

This is an **open-source** set of **techniques** with **reproducible benchmarks** to evaluate primer generation strategies. 

## Project goals and metrics
Generate as many primers (or extract as many primers from a fixed list) as possible with these constraints:

1. Each primer has L=20 nucleotides
1. Each primer must have a ratio of ((G + C) / (A + T + G + C)) between 0.45 and 0.55 ("balanced GC content")
1. Any two pairs of primers must have an edit distance >= 8 (0.4 * L)

Metrics:

1. **Primer Extraction Ratio (PER%)**: of a given list of valid (except edit distance) primers, filter it down to the largest possible set of edit-distance valid primers.
2. **Primer Generation Speed**: any methods must remain tractable.

## Setup
1. `git clone` this repo.
1. If you use conda, install dependencies with `conda env create -f environment.yml`. 
1. If you use pip, see the following section.

### Pip packages to install (no requirements.txt yet)
1. `biopython`
2. `editdistance` https://github.com/roy-ht/editdistance
3. `polyleven` https://github.com/fujimotos/polyleven (faster levenshtein than `editdistance`)
4. `networkx` - for graph algorithms

## Running
- (Extractors) Run an extractor with e.g., `python -m primergen.extractors.delob_mindegree`
- (Generators) Run a generator with e.g., `python -m primergen.generators.random_gc`


## High-level code organization and concepts 
- Modules in `extractors`: Set of techniques that take in a valid (except edit distance) list of primers, and filters it down to an edit-distance valid set of primers.
- Modules in `generators`: Set of techniques that continually generate new valid primers (no fixed input list).
- Files in `input` represent fixed lists of primers that are used to test `extractors` consistently. `*-edges.txt` or `*.pkl` are files that cache all pairs of edit distance calculations to make testing graph-based methods faster (no need to recompute all-pairs edit distance each time).
- Files in `output` represent the result of each run, which are timestamped and include the name of technique that generated them.

## Techniques implemented
### Extractors
1. `naive_clique`: Optimal algo (assuming initial number of primers is large enough), find maximum size clique in graph, where vertices are primers and edges represent two primers that are *far enough away* in edit distance
1. `approx-mis`: Approximation of maximum independent set in graph, where vertices are primers and edges represent two primers that are *too close* in edit distance.
1. `delob`: Same representation as `approx-mis`. We remove 1 vertex (primer) randomly from graph each time, then remove all of its neighbors, repeat until no more edges.
1. `delob-mindegree`: Only change: we remove vertices with the minimum degree first.
1. `delob-mindegree-neighbors`: Only change: we remove vertices that have neighbors of neighbors with the lowest average degree. 
1. `greedy`: Choose a random primer from input set, evaluate edit distance to all edit-distance valid primers so far, add it to the edit-distance valid primers set if all distances are more than the threshold.

### Generators
1. `random`: Completely random primer generation, check GC, check edit dist against all accepted primers so far, add to pool if so.
1. `random_gc`: Generate primers that already have balanced GC content, rest of steps are identical
1. `random_gc_frequencies`: Generate primers that are likely to be as different as possible from the average primer in the pool for each position. Reroll primer generation until we get balanced GC content.
1. `random_gc_frequencies_noreroll`: Same as before but we don't reroll for valid GC content, we fix it in constant time.




