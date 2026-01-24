#!/usr/bin/env python3
"""
Extract timed out PrimeWitness entries from SHA order 3 data files.

Goes through all files in sha_order3_data and extracts only the lines
where PrimeWitness computation timed out, outputting them into a single file.
"""

import os
import glob

INPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "sha_order3_data")
OUTPUT_FILE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "sha3_timeouts.processed.txt")


def extract_timeouts(input_path):
    """Extract lines with PrimeWitness timeout from a single file."""
    timeout_lines = []
    with open(input_path, 'r') as f:
        for line in f:
            # Keep only data lines (no '[' for headers) that have "(timeout)"
            if line.strip() and '[' not in line and '(timeout)' in line:
                # Remove the "PrimeWitness: 0 (timeout)" suffix
                cleaned = line.split('PrimeWitness:')[0].rstrip() + '\n'
                timeout_lines.append(cleaned)
    return timeout_lines


def main():
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    input_files = sorted(glob.glob(os.path.join(INPUT_DIR, "sha3eqns.*.txt")))

    if not input_files:
        print("No input files found!")
        return

    all_timeouts = []
    for input_path in input_files:
        basename = os.path.basename(input_path)
        timeouts = extract_timeouts(input_path)
        all_timeouts.extend(timeouts)
        print(f"Processed {basename}: {len(timeouts)} timeouts")

    with open(OUTPUT_FILE, 'w') as f:
        f.writelines(all_timeouts)

    print(f"\nDone! Found {len(all_timeouts)} total timeouts")
    print(f"Output written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
