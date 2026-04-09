#!/usr/bin/env python3
"""
Merge PrimeWitness data from witnessed files back into the filtered input files.

This script:
1. For each file in sha_order3_no3_iso, reads the corresponding witnessed file
   to get PrimeWitness values
2. Appends PrimeWitness values to each torsor data row
3. Writes the result to sha_order3_no3_iso_witnessed
"""

import os
import re

PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))
INPUT_DIR = os.path.join(PROJECT_ROOT, "data", "sha_order3_no3_iso")
WITNESSED_DIR = os.path.join(PROJECT_ROOT, "data", "sha_order3_processed_witnessed")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "sha_order3_no3_iso_witnessed")


def is_header_line(line):
    """Check if a line is a header line (contains bracket notation for curve)."""
    return '[' in line and ']' in line


def is_data_row(line):
    """Check if a line is a data row (non-empty, not a header)."""
    stripped = line.strip()
    return stripped and not is_header_line(stripped)


def extract_primewitness(line):
    """Extract the PrimeWitness portion from a witnessed line."""
    match = re.search(r'(PrimeWitness: \S+(?:\s+\(no witness found\))?)', line)
    if match:
        return match.group(1)
    return None


def get_witnessed_filename(input_filename):
    """Convert input filename to witnessed filename."""
    # sha3eqns.00000-09999.txt -> witnessed.notimeout.sha3eqns.00000-09999.processed.txt
    base = input_filename.replace('.txt', '')
    return f"witnessed.notimeout.{base}.processed.txt"


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
                    output_lines.append(f"{line.rstrip()} {primewitness_values[pw_index]}\n")
                    pw_index += 1
                else:
                    print(f"Warning: Ran out of PrimeWitness values at row {pw_index}")
                    output_lines.append(line)
            else:
                # Header or empty line - keep as is
                output_lines.append(line)

    if pw_index != len(primewitness_values):
        print(f"Warning: Used {pw_index} of {len(primewitness_values)} PrimeWitness values")

    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    return pw_index


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    input_files = sorted([f for f in os.listdir(INPUT_DIR) if f.endswith('.txt')])

    total_rows = 0
    for input_filename in input_files:
        input_path = os.path.join(INPUT_DIR, input_filename)
        witnessed_filename = get_witnessed_filename(input_filename)
        witnessed_path = os.path.join(WITNESSED_DIR, witnessed_filename)
        output_path = os.path.join(OUTPUT_DIR, input_filename)

        if not os.path.exists(witnessed_path):
            print(f"Warning: No witnessed file for {input_filename}")
            continue

        rows = process_file(input_path, witnessed_path, output_path)
        total_rows += rows
        print(f"Processed {input_filename}: {rows} data rows")

    print(f"\nTotal: {total_rows} data rows processed")
    print(f"Output written to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
