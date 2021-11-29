#!/usr/bin/env python3
from util import write_primers
from check import *
from util import write_primers, random_primer
import time
import random

PRIMERS_PER_SECOND_PERIOD_SEC = 1

class BasePrimerGenerator:
    """
    Generic class for all primer generators.
    Helps to unify metrics to compare different primer generation methods, saving to file, etc.
    """
    def __init__(self, target=2000, strategy="base"):
        self.target = target
        self.primers = []
        self.iterations = 0
        self.gc_errors = 0
        self.edit_errors = 0
        self.strategy = strategy
        # Tracking primers per second
        self.prev_primers_time = 0
        self.prev_primers_count = 0
        self.primers_per_second = 0

    def execute(self):
        """
        Method to call by any runners (e.g., main methods).
        Starts stats, runs generate(), then finishes up and writes to file
        """
        self.start()
        try:
            self.generate()
        except KeyboardInterrupt:
            print(f"Exited early with {len(self.primers)} primers!")
        self.finish()

    def start(self):
        """
        Call just before starting primer generation
        """
        self.start_time = time.process_time()
        self.prev_primers_time = self.start_time

    def generate(self):
        """
        NOTE: Inheriting classes will implement this
        """
        raise NotImplementedError


    def finish(self):
        """
        Call after filling up primer list
        """
        # Timing stats
        end_time = time.process_time()
        total_time = end_time - self.start_time

        # Print stats
        print(self.primers)
        print(f"Total time (sec): {total_time}")

        # Write to file
        write_primers(self.primers, total_time_sec=total_time, strategy=self.strategy)

        # Check the primers for errors (after writing, so we can check it ourselves later)
        print(f"Validating primers after saving...")
        valid = are_primers_valid(self.primers)
        if not valid:
            print("ERROR: Found invalid primers in final library! Error with algorithm.")
        else:
            print("Primers validated (all good) before writing!")


    def found_new_primer(self, primer):
        """
        Called when we confirm a new primer is found
        """
        self.primers.append(primer)

    def found_new_primers(self, primers):
        """
        Called when we confirm a new primer is found
        """
        self.primers.extend(primers)


    def new_iteration(self):
        self.iterations += 1

    def new_gc_error(self):
        self.gc_errors += 1

    def new_edit_error(self):
        self.edit_errors += 1

    def print_metrics(self):
        # Elapsed
        cur_time =  time.process_time()
        elapsed_sec = cur_time - self.start_time
        elapsed_min = int(elapsed_sec / 60)
        elapsed_min_sec = int(elapsed_sec % 60)
        # Count primers per second
        if cur_time - self.prev_primers_time > PRIMERS_PER_SECOND_PERIOD_SEC:
            # We've passed the period we want to compute PPS for
            current_num_primers = len(self.primers)
            self.primers_per_second = round((current_num_primers - self.prev_primers_count) / (cur_time - self.prev_primers_time), 2)
            self.prev_primers_time = cur_time
            self.prev_primers_count = current_num_primers
        print(
            f"CPU Elapsed: {elapsed_min} min {elapsed_min_sec} sec\tIteration: {self.iterations}\tPrimers: {len(self.primers)}\tPPS: {self.primers_per_second}\tGC invalid: {self.gc_errors}\tEdit invalid\t{self.edit_errors}\n"
        )
