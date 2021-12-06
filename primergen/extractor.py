#!/usr/bin/env python3

# PRIMER_FILE = "20211205-141133-100000-primers.txt"
# PRIMER_FILE = "20211204-141808-10000-primers.txt"
# PRIMER_FILE = "20211204-174137-5000-primers.txt"
PRIMER_FILE = "20211205-182302-2000-primers.txt"


def get_primers():
    initial_primers = []
    with open(PRIMER_FILE, "r") as f:
        initial_primers = f.read().splitlines()
    return initial_primers
