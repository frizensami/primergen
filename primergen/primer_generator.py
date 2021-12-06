#!/usr/bin/env python3

from check import *
from util import *

import random
import numpy

if __name__ == "__main__":
    # Seed random same for all extractors
    random.seed(246)
    numpy.random.seed(4812)
    # Generating large number N of initial primers
    N = int(RANDOM_PRIMERS_TO_GENERATE)
    primers = generate_n_primers(n=N)
    filename = get_filestamp(suffix=f"{len(primers)}-primers.txt")
    with open(filename, "w") as f:
        f.write("\n".join(primers))
