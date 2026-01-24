# On the visibility category of the Shafarevich-Tate group

## Directory Structure

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

