#!/usr/bin/env python3
"""
Extract Cremona labels from the raw SHA3 equation files.

Each block starts with a header line like:
    182 b 3 [1,0,0,-15663,-755809] 0 1 9

The first three tokens (e.g. '182', 'b', '3') form the Cremona label '182b3'.
"""

import os
import glob

RAW_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "sha_order3_raw")
OUTPUT_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "cremona_labels.txt")


def extract_labels():
    labels = []
    for filepath in sorted(glob.glob(os.path.join(RAW_DIR, "*.txt"))):
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Header lines contain '[' (the Weierstrass coefficients)
                if '[' in line:
                    parts = line.split()
                    label = parts[0] + parts[1] + parts[2]
                    labels.append(label)

    with open(OUTPUT_PATH, 'w') as f:
        for label in labels:
            f.write(label + '\n')

    print(f"Extracted {len(labels)} Cremona labels to {OUTPUT_PATH}")


if __name__ == "__main__":
    extract_labels()
