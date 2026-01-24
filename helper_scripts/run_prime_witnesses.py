#!/usr/bin/env python3
"""
Wrapper script to run PrimeWitness computation with per-case timeouts.

Usage:
  python3 run_prime_witnesses.py <input_filename>

Where input_filename is just the filename (e.g., sha3eqns.00000-09999.processed.txt)
located in data/sha_order3_processed/

To run in parallel:
  ls data/sha_order3_processed | parallel -j30 "python3 run_prime_witnesses.py {}"
"""

import os
import sys
import subprocess

# Timeout in seconds for each case
TIMEOUT_SECONDS = 5

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(PROJECT_ROOT, "data", "sha_order3_processed")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "data", "sha_order3_processed_witnessed")

# Magma script to compute prime witness for a single case
MAGMA_SINGLE_CASE = '''
AttachSpec("magma/spec");
L := {coeffs};
try
    C := CMcurveForIndex3Torsor(L);
    p := PrimeWitnessForMinimality(C);
    print p;
catch e
    print 0;
end try;
quit;
'''


def compute_prime_witness(coeffs_str):
    """
    Run Magma to compute prime witness for a single case.
    Returns (prime, timed_out) where prime is int (0 on failure) and timed_out is bool.
    """
    magma_code = MAGMA_SINGLE_CASE.format(coeffs=coeffs_str)

    try:
        result = subprocess.run(
            ["magma", "-b"],
            input=magma_code,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SECONDS,
            cwd=PROJECT_ROOT
        )
        # Parse the output - should be just the prime number
        output = result.stdout.strip()
        # Get last non-empty line (in case there's other output)
        lines = [l.strip() for l in output.split('\n') if l.strip()]
        if lines:
            return int(lines[-1]), False
        return 0, False
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT for input {coeffs_str}", file=sys.stderr)
        return 0, True
    except Exception as e:
        print(f"FAILURE for input {coeffs_str}: {e}", file=sys.stderr)
        return 0, False


def process_file(input_filename):
    """Process all lines in an input file."""
    input_path = os.path.join(INPUT_DIR, input_filename)
    output_filename = "witnessed." + input_filename
    output_path = os.path.join(OUTPUT_DIR, output_filename)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with open(input_path, 'r') as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    with open(output_path, 'w') as out:
        for i, line in enumerate(lines, 1):
            # Convert "1 2 3 4 5 6 7 8 9 10" to "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]"
            parts = line.split()
            coeffs_str = "[" + ", ".join(parts) + "]"

            print(f"Processing line {i}: {coeffs_str}")
            p, timed_out = compute_prime_witness(coeffs_str)

            if timed_out:
                out.write(f"{line} PrimeWitness: {p} (timeout)\n")
            else:
                out.write(f"{line} PrimeWitness: {p}\n")
            out.flush()  # Flush after each line so progress is visible

    print(f"Done processing {input_filename}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_filename>", file=sys.stderr)
        sys.exit(1)

    process_file(sys.argv[1])
