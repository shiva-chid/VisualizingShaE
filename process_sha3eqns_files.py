#!/usr/bin/env python3
"""
Process SHA order 3 equation files.

Removes elliptic curve header lines (containing '[') and blank lines,
keeping only the ternary cubic coefficient lines.
"""

import os
import glob

INPUT_DIR = os.path.join(os.path.dirname(__file__), "data", "sha_order3_input")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "data", "sha_order3_processed")


def process_file(input_path, output_path):
    """Process a single file, removing curve headers and blank lines."""
    with open(input_path, 'r') as f:
        lines = f.readlines()

    # Keep only lines that:
    # - Are not blank (after stripping whitespace)
    # - Do not contain '[' (which indicates elliptic curve header lines)
    processed_lines = [
        line for line in lines
        if line.strip() and '[' not in line
    ]

    with open(output_path, 'w') as f:
        f.writelines(processed_lines)

    return len(processed_lines)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    input_files = sorted(glob.glob(os.path.join(INPUT_DIR, "sha3eqns.*.txt")))

    if not input_files:
        print("No input files found!")
        return

    total_lines = 0
    for input_path in input_files:
        basename = os.path.basename(input_path)
        # sha3eqns.00000-09999.txt -> sha3eqns.00000-09999.processed.txt
        output_name = basename.replace(".txt", ".processed.txt")
        output_path = os.path.join(OUTPUT_DIR, output_name)

        num_lines = process_file(input_path, output_path)
        total_lines += num_lines
        print(f"Processed {basename} -> {output_name} ({num_lines} lines)")

    print(f"\nDone! Processed {len(input_files)} files ({total_lines} total lines)")


if __name__ == "__main__":
    main()
