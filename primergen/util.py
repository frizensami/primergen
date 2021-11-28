#!/usr/bin/env python3

import time
import os

DATA_FOLDER = "data/"


def get_filestamp(suffix=None):
    timestr = time.strftime("%Y%m%d-%H%M%S")
    if suffix:
        timestr += f"-{suffix}"
    return timestr


def write_primers(primers, total_time_sec=-1, strategy=""):
    # Write primer
    os.makedirs(DATA_FOLDER, exist_ok=True)
    with open(
        os.path.join(DATA_FOLDER, get_filestamp(suffix=f"{strategy}.txt")), "w"
    ) as f:
        f.write(f"Total time (seconds):\t{total_time_sec}\n")
        f.write("\n".join(primers))
