#!/usr/bin/env python3
"""
Merge PrimeWitness data from processed files back into the original input files.

This script:
1. Copies sha_order3_input to sha_order3_data
2. For each file, reads the corresponding witnessed file to get PrimeWitness values
3. Appends PrimeWitness values to each data row in the copied file
"""

import os
import shutil
import re

PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))
INPUT_DIR = os.path.join(PROJECT_ROOT, "data", "sha_order3_input")
WITNESSED_DIR = os.path.join(PROJECT_ROOT, "data", "sha_order3_processed_witnessed")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "sha_order3_data")


def is_header_line(line):
    """Check if a line is a header line (contains bracket notation for curve)."""
    return '[' in line and ']' in line


def is_data_row(line):
    """Check if a line is a data row (non-empty, not a header)."""
    stripped = line.strip()
    return stripped and not is_header_line(stripped)


def extract_primewitness(line):
    """Extract the PrimeWitness portion from a witnessed line."""
    match = re.search(r'(PrimeWitness: \S+(?:\s+\(timeout\))?)', line)
    if match:
        return match.group(1)
    return None


def get_witnessed_filename(input_filename):
    """Convert input filename to witnessed filename."""
    # sha3eqns.00000-09999.txt -> witnessed.sha3eqns.00000-09999.processed.txt
    base = input_filename.replace('.txt', '')
    return f"witnessed.{base}.processed.txt"


def process_file(input_path, witnessed_path, output_path):
    """Merge PrimeWitness data from witnessed file into input file."""

    # Read all PrimeWitness values from witnessed file
    primewitness_values = []
    with open(witnessed_path, 'r') as f:
        for line in f:
            pw = extract_primewitness(line)
            if pw:
                primewitness_values.append(pw)

    # Read input file and merge
    output_lines = []
    pw_index = 0

    with open(input_path, 'r') as f:
        for line in f:
            if is_data_row(line):
                # This is a data row - append PrimeWitness
                if pw_index < len(primewitness_values):
                    # Remove trailing whitespace/newline, add PrimeWitness, add newline
                    output_lines.append(f"{line.rstrip()} {primewitness_values[pw_index]}\n")
                    pw_index += 1
                else:
                    print(f"Warning: Ran out of PrimeWitness values at row {pw_index}")
                    output_lines.append(line)
            else:
                # Header or empty line - keep as is
                output_lines.append(line)

    # Check if we used all PrimeWitness values
    if pw_index != len(primewitness_values):
        print(f"Warning: Used {pw_index} of {len(primewitness_values)} PrimeWitness values")

    # Write output file
    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    return pw_index


def main():
    # Create output directory (fresh copy)
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)

    # Get list of input files
    input_files = sorted([f for f in os.listdir(INPUT_DIR) if f.endswith('.txt')])

    total_rows = 0
    for input_filename in input_files:
        input_path = os.path.join(INPUT_DIR, input_filename)
        witnessed_filename = get_witnessed_filename(input_filename)
        witnessed_path = os.path.join(WITNESSED_DIR, witnessed_filename)
        output_path = os.path.join(OUTPUT_DIR, input_filename)

        if not os.path.exists(witnessed_path):
            print(f"Warning: No witnessed file for {input_filename}")
            shutil.copy(input_path, output_path)
            continue

        rows = process_file(input_path, witnessed_path, output_path)
        total_rows += rows
        print(f"Processed {input_filename}: {rows} data rows")

    print(f"\nTotal: {total_rows} data rows processed")
    print(f"Output written to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
