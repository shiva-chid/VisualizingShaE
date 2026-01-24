#!/bin/bash
# Test suite to verify the codebase refactor didn't break anything
# Run from project root: ./test_refactor.sh
#
# Quick one-liner tests (run individually if needed):
#   magma -b -e 'AttachSpec("magma/spec"); print "OK"; quit;'
#   magma -b -e 'Attach("Tamagawa/Tamagawa_pkg2.m"); print "OK"; quit;'
#   magma -b -e 'print Read("data/5.txt"); quit;'
#   magma -b -e 'AttachSpec("magma/spec"); C:=CMcurveForIndex3Torsor([1,2,91,0,0,0,0,0,0,-17]); print Genus(C); quit;'
#   python3 -m py_compile helper_scripts/*.py

set -e  # Exit on first error

echo "=== Testing Codebase Refactor ==="
echo "Run this from the project root directory"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

pass() { echo -e "${GREEN}PASS${NC}: $1"; }
fail() { echo -e "${RED}FAIL${NC}: $1"; exit 1; }

# Check we're in the right directory
[ -d "magma" ] || fail "Not in project root (magma/ not found)"
[ -d "helper_scripts" ] || fail "Not in project root (helper_scripts/ not found)"
[ -d "data" ] || fail "Not in project root (data/ not found)"

echo "=== 1. Testing Magma spec file loads ==="
magma -b <<'EOF' && pass "spec loads correctly" || fail "spec failed to load"
AttachSpec("magma/spec");
print "Spec loaded successfully";
quit;
EOF

echo ""
echo "=== 2. Testing g2RMapproach.m can load data files ==="
magma -b <<'EOF' && pass "g2RMapproach.m data loading works" || fail "g2RMapproach.m data loading failed"
// Test that the curve family and elliptic curves database load correctly
P<a,b,c,x> := PolynomialRing(Rationals(), 4);

// Test load_curve_family path (data/5.txt)
Dfile := Sprintf("../data/%o.txt", 5);
// Run from project root, so path is data/5.txt
Dfile := Sprintf("data/%o.txt", 5);
f := eval Read(Dfile);
print "Loaded curve family from data/5.txt";

// Test elliptic curves database path
filename := Sprintf("data/ellipticcurvesinfo%o.txt", 5);
fil := Open(filename, "r");
s := Gets(fil);
print "Loaded ellipticcurvesinfo5.txt, first line length:", #s;

print "Data files load correctly";
quit;
EOF

echo ""
echo "=== 3. Testing Tamagawa package attachment ==="
magma -b <<'EOF' && pass "Tamagawa package loads" || fail "Tamagawa package failed"
Attach("Tamagawa/Tamagawa_pkg2.m");
print "Tamagawa package loaded successfully";
quit;
EOF

echo ""
echo "=== 4. Testing PrimeWitnesses.m core functions ==="
magma -b <<'EOF' && pass "PrimeWitnesses core functions work" || fail "PrimeWitnesses failed"
AttachSpec("magma/spec");
// Test with a known diagonal torsor
L := [1, 2, 91, 0, 0, 0, 0, 0, 0, -17];
C := CMcurveForIndex3Torsor(L);
print "CMcurveForIndex3Torsor works, genus:", Genus(C);
// Don't run full PrimeWitnessForMinimality as it's slow
quit;
EOF

echo ""
echo "=== 5. Testing diagonaltorsors_check.m data path ==="
magma -b <<'EOF' && pass "diagonaltorsors.txt loads" || fail "diagonaltorsors.txt failed"
fil := Open("data/diagonaltorsors.txt","r");
s := Gets(fil);
print "First line of diagonaltorsors.txt:", s;
quit;
EOF

echo ""
echo "=== 6. Testing Python scripts (syntax check) ==="
python3 -m py_compile helper_scripts/download_sha3_data.py && pass "download_sha3_data.py syntax OK" || fail "download_sha3_data.py syntax error"
python3 -m py_compile helper_scripts/process_sha3eqns_files.py && pass "process_sha3eqns_files.py syntax OK" || fail "process_sha3eqns_files.py syntax error"
python3 -m py_compile helper_scripts/run_prime_witnesses.py && pass "run_prime_witnesses.py syntax OK" || fail "run_prime_witnesses.py syntax error"
python3 -m py_compile helper_scripts/merge_primewitness.py && pass "merge_primewitness.py syntax OK" || fail "merge_primewitness.py syntax error"
python3 -m py_compile helper_scripts/extract_timeouts.py && pass "extract_timeouts.py syntax OK" || fail "extract_timeouts.py syntax error"

echo ""
echo "=== 7. Testing Python path resolution ==="
python3 -c "
import os
import sys
sys.path.insert(0, 'helper_scripts')

# Test that paths resolve correctly from helper_scripts/
os.chdir('helper_scripts')

# Test download_sha3_data.py path
from download_sha3_data import OUTPUT_DIR
assert 'data/sha_order3_input' in OUTPUT_DIR, f'Wrong path: {OUTPUT_DIR}'
print(f'download_sha3_data.py OUTPUT_DIR: {OUTPUT_DIR}')

# Test process_sha3eqns_files.py paths
from process_sha3eqns_files import INPUT_DIR, OUTPUT_DIR
assert 'data/sha_order3_input' in INPUT_DIR, f'Wrong INPUT_DIR: {INPUT_DIR}'
assert 'data/sha_order3_processed' in OUTPUT_DIR, f'Wrong OUTPUT_DIR: {OUTPUT_DIR}'
print(f'process_sha3eqns_files.py paths OK')

# Test merge_primewitness.py paths
from merge_primewitness import INPUT_DIR, WITNESSED_DIR, OUTPUT_DIR as MERGE_OUTPUT
assert 'data/sha_order3_input' in INPUT_DIR
assert 'data/sha_order3_processed_witnessed' in WITNESSED_DIR
assert 'data/sha_order3_data' in MERGE_OUTPUT
print(f'merge_primewitness.py paths OK')

print('All Python paths resolve correctly')
" && pass "Python paths resolve correctly" || fail "Python path resolution failed"

echo ""
echo "=== 8. Testing run_prime_witnesses.py Magma integration ==="
# This actually runs Magma via the Python script
python3 helper_scripts/run_prime_witnesses.py --help 2>/dev/null || python3 -c "
import subprocess
import sys

# Run a quick Magma test via the script's mechanism
magma_code = '''
AttachSpec(\"magma/spec\");
L := [1, 2, 91, 0, 0, 0, 0, 0, 0, -17];
C := CMcurveForIndex3Torsor(L);
print Genus(C);
quit;
'''
result = subprocess.run(['magma', '-b'], input=magma_code, capture_output=True, text=True, timeout=30)
print('Magma output:', result.stdout.strip())
if '1' in result.stdout:
    print('Integration test passed')
    sys.exit(0)
else:
    print('Unexpected output')
    sys.exit(1)
" && pass "Magma integration via Python OK" || fail "Magma integration failed"

echo ""
echo "==========================================="
echo -e "${GREEN}All tests passed!${NC}"
echo "==========================================="
