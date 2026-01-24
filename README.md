# On the visibility category of the Shafarevich-Tate group

This is code to accompany the paper of the same name by Barinder S. Banwait, Jerson Caro, and Shiva Chidambaram.

## Directory Structure

```
.
├── magma/                      # Core Magma scripts
│   ├── spec                    # Package specification
│   ├── PrimeWitnesses.m
│   ├── g2RMapproach.m
│   └── ...
├── helper_scripts/             # Python/Sage data processing
│   ├── download_sha3_data.py
│   ├── process_sha3eqns_files.py
│   ├── run_prime_witnesses.py
│   ├── merge_primewitness.py
│   └── ...
├── data/
│   ├── 5.txt, 13.txt           # Parametric curve families
│   ├── ellipticcurvesinfo*.txt # Elliptic curve databases
│   ├── diagonaltorsors.txt
│   └── sha_order3_data/        # ~124k SHA order 3 entries
├── Tamagawa/                   # Van Bommel's Tamagawa package
├── magma_scripts/              # Standalone examples
└── test_refactor.sh            # Test suite
```

### `magma/`
Magma scripts for computing Cremona-Mazur curves, prime witnesses for minimality, and verifying mod-ell congruences between genus 2 Jacobians and elliptic curves. The `spec` file defines the core package (FractionalLinearTransformation, KummerElementAndTransformation, Minimality_CremonaMazur_index3).

### `helper_scripts/`
Python and Sage scripts for data processing. Includes a pipeline for downloading Tom Fisher's SHA order 3 data, computing prime witnesses in parallel, and merging results. Also contains LMFDB query scripts for finding candidate curves.

### `data/`
Contains parametric curve families (`5.txt`, `13.txt`), precomputed elliptic curve databases (`ellipticcurvesinfo*.txt`), diagonal torsor data, and the SHA order 3 dataset with computed prime witnesses (~124k entries across 30 files).

### `Tamagawa/`
Raymond van Bommel's Magma package for computing Tamagawa numbers of genus 2 curves via regular models.

### `magma_scripts/`
Standalone Magma examples for gluing constructions and Fisher's algebras.

## Section 3 verification

The main code for Section 3 is `magma/g2RMapproach.m`. When running it, look for the following, indicating success:

=========================================================
SUCCESS: Tamagawa number computation
=========================================================
Curve coefficients: [ 1, 6, 5, 12, 2, 0, 1 ]
Rank bounds: [ 4, 4 ]
Congruent elliptic curves (a-invariants):[
[ 0, -1, 0, -169321, -28379327 ]]
Tamagawa numbers: [ <2, 9>, <3, 4>, <29, 1> ]

## Testing

Run the full test suite from the project root:

```bash
./test_refactor.sh
```

Or run individual tests manually:

```bash
# Test spec file loads
magma -b -e 'AttachSpec("magma/spec"); print "OK"; quit;'

# Test Tamagawa package
magma -b -e 'Attach("Tamagawa/Tamagawa_pkg2.m"); print "OK"; quit;'

# Test data file paths
magma -b -e 'print Read("data/5.txt"); quit;'

# Test core function
magma -b -e 'AttachSpec("magma/spec"); C:=CMcurveForIndex3Torsor([1,2,91,0,0,0,0,0,0,-17]); print Genus(C); quit;'

# Test Python script syntax
python3 -m py_compile helper_scripts/*.py
```

