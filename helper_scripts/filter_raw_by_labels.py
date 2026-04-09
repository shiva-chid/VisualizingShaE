#!/usr/bin/env python3
"""
Filter the raw SHA3 equation files, keeping only curves whose Cremona labels
appear in cremona_labels_no_3_iso.txt.

Reads blocks from data/sha_order3_raw/ and writes the filtered blocks
to data/sha_order3_no3_iso/, preserving the same file structure.
"""

import os
import glob

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
RAW_DIR = os.path.join(BASE_DIR, "data", "sha_order3_raw")
LABELS_PATH = os.path.join(BASE_DIR, "data", "cremona_labels_no_3_iso.txt")
OUTPUT_DIR = os.path.join(BASE_DIR, "data", "sha_order3_no3_iso")


def load_labels():
    with open(LABELS_PATH) as f:
        return set(line.strip() for line in f if line.strip())


def parse_blocks(filepath):
    """Yield (cremona_label, block_text) for each block in a raw file."""
    with open(filepath) as f:
        lines = f.readlines()

    block_lines = []
    label = None
    for line in lines:
        stripped = line.strip()
        if stripped == "":
            if label is not None:
                yield label, "".join(block_lines)
                block_lines = []
                label = None
        else:
            if label is None and "[" in stripped:
                parts = stripped.split()
                label = parts[0] + parts[1] + parts[2]
            block_lines.append(line)

    # Handle last block if file doesn't end with a blank line
    if label is not None:
        yield label, "".join(block_lines)


def filter_files():
    labels = load_labels()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    total_kept = 0
    total_skipped = 0

    for filepath in sorted(glob.glob(os.path.join(RAW_DIR, "*.txt"))):
        filename = os.path.basename(filepath)
        kept_blocks = []

        for label, block_text in parse_blocks(filepath):
            if label in labels:
                kept_blocks.append(block_text)
                total_kept += 1
            else:
                total_skipped += 1

        output_path = os.path.join(OUTPUT_DIR, filename)
        with open(output_path, "w") as f:
            f.write("\n".join(kept_blocks))
            if kept_blocks:
                f.write("\n")

        print(f"{filename}: kept {len(kept_blocks)} blocks")

    print(f"\nTotal: kept {total_kept}, skipped {total_skipped}")


if __name__ == "__main__":
    filter_files()
