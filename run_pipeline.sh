#!/bin/bash
#
# Master script to reproduce the full data pipeline.
#
# A reviewer can run this script to recreate all derived data from scratch.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Step 1: Download raw SHA3 data from Fisher's website ==="
python3 helper_scripts/download_sha3_data.py

echo ""
echo "=== Step 2: Filter to curves without 3-isogenies ==="
python3 helper_scripts/filter_raw_by_labels.py

echo ""
echo "=== Step 3: Process filtered files (extract ternary cubic coefficients) ==="
python3 helper_scripts/process_sha3eqns_files.py

echo ""
echo "=== Step 4: Compute prime witnesses (requires Magma) ==="
echo "Running PrimeWitness on all processed files in parallel..."
ls data/sha_order3_processed | parallel -j30 "python3 helper_scripts/run_prime_witnesses.py {}"
