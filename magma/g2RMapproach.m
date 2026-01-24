/*
    g2RMapproach.m

    This script searches for genus 2 curves C whose Jacobian J(C) has
    ell-torsion that is congruent (mod ell) to an elliptic curve E.

    The approach:
    1. Load a database of elliptic curves with precomputed Frobenius traces
    2. Load a parametric family of genus 2 curves (polynomial f in a,b,c,x)
    3. Search through low-height rational specializations (a,b,c)
    4. For each specialization, filter candidate elliptic curves by matching
       Frobenius traces modulo ell
    5. Further filter by conductor exponent at ell (must be 0 for good reduction)
    6. Verify congruences rigorously and compute Mordell-Weil rank bounds
    7. Compute Tamagawa numbers for the resulting curves
*/

//=============================================================================
// GLOBAL SETUP
//=============================================================================

ZZ := Integers();
P<a,b,c,x> := PolynomialRing(Rationals(), 4);
P1<t> := PolynomialRing(Rationals());

//=============================================================================
// We need Raymond van Bommel's Tamagawa package for Tamagawa numbers
//=============================================================================

Attach("../Tamagawa/Tamagawa_pkg2.m");

//=============================================================================
// CONFIGURATION PARAMETERS
//=============================================================================

// Height bound for generating low-height rationals
HEIGHT_BOUND := 10;

// Prime threshold for checking bad primes of discriminant
PRIMES_THRESHOLD := 10^4;

// Bound for computing Frobenius traces
PRIMES_BOUND := 10^2;

// Maximum number of (a,b,c) triples to search
LOW_HT_RATIONALS_BOUND := 25000;

// The prime ell for the mod-ell representation
ELL := 5;

// Parameter D for the curve family (loads from D.txt)
D_PARAM := 5;

// Precision for p-adic computations
PADIC_PRECISION := 100;

// Higher prime bound for rigorous congruence verification
HIGH_PRIMES_BOUND := 10^3;

//=============================================================================
// UTILITY FUNCTIONS
//=============================================================================

/*
    Generate a sorted list of low-height rationals a/b where:
    - |a| <= htbd, 1 <= b <= htbd
    - gcd(a,b) = 1
    Sorted by string length (simpler rationals first)
*/
function generate_low_height_rationals(htbd)
    rats := Setseq({a/b : a in [-htbd..htbd], b in [1..htbd] | GCD(a,b) eq 1});
    rats := Sort(rats, func<a,b | #Sprint(a) - #Sprint(b)>);
    return rats;
end function;

/*
    Load elliptic curves from a precomputed file.
    File format: ainvariants:conductor:label:ap_values (colon-separated)

    Returns:
        Es     - list of elliptic curves
        apEs   - list of associative arrays mapping p -> a_p(E)
*/
function load_elliptic_curves(filename, primesbound)
    fil := Open(filename, "r");
    s := Gets(fil);
    Es := [];
    apEs := [];
    primes := PrimesUpTo(primesbound);

    while not IsEof(s) do
        s_split := Split(s, ":");
        ainvs := [StringToInteger(x) : x in Split(s_split[1], "[], ")];
        // cond := StringToInteger(s_split[2]);  // conductor (unused here)
        // label := s_split[3];                   // curve label (unused here)
        aps := [StringToInteger(x) : x in Split(s_split[4], "[], ")];

        // Build associative array for Frobenius traces
        apsdict := AssociativeArray();
        for i := 1 to #aps do
            apsdict[primes[i]] := aps[i];
        end for;

        Append(~Es, EllipticCurve(ainvs));
        Append(~apEs, apsdict);
        s := Gets(fil);
    end while;

    return Es, apEs;
end function;

/*
    Load the parametric polynomial f from a file.
    The file should contain a polynomial in variables a,b,c,t.
*/
function load_curve_family(D)
    Dfile := Sprintf("../data/%o.txt", D);
    f := eval Read(Dfile);
    return f;
end function;

//=============================================================================
// DISCRIMINANT AND VALIDITY CHECKS
//=============================================================================

/*
    Compute the discriminant of a hyperelliptic curve as an integer.
    Returns the product of numerator and denominator (for coprimality checks).
*/
function discriminant_as_integer(C)
    discC := Discriminant(C);
    return ZZ!(Numerator(discC) * Denominator(discC));
end function;

/*
    Check if the discriminant is valid for our search:
    - Not divisible by ell
    - All prime factors are below primesthreshold

    Returns: true if valid, false otherwise
*/
function is_discriminant_valid(discCmod, ell, primesthreshold)
    // Must have good reduction at ell
    if discCmod mod ell eq 0 then
        return false;
    end if;

    // Check that all bad primes are below threshold
    badprimes := [p : p in PrimesUpTo(primesthreshold) | discCmod mod p eq 0];
    product_of_prime_powers := &*([1] cat [p^Valuation(discCmod, p) : p in badprimes]);

    if product_of_prime_powers ne Abs(discCmod) then
        return false;
    end if;

    return true;
end function;

//=============================================================================
// FROBENIUS TRACE COMPUTATIONS
//=============================================================================

/*
    Compute point counts on a hyperelliptic curve over F_p and F_{p^2}.

    Returns:
        n1 - number of points over F_p
        n2 - number of points over F_{p^2}
        or false if computation fails
*/
function compute_point_counts(C, p)
    try
        Cp := ChangeRing(C, GF(p));
        Cp2 := ChangeRing(C, GF(p,2));
        return #Cp, #Cp2;
    catch e;
        return false, false;
    end try;
end function;

/*
    Compute the Frobenius trace tp and norm np for a genus 2 curve.

    For a genus 2 curve, the characteristic polynomial of Frobenius is:
        X^4 - t_p X^3 + (n_p + 2p) X^2 - p t_p X + p^2

    where t_p = p + 1 - #C(F_p) and n_p is computed from #C(F_p) and #C(F_{p^2}).
*/
function compute_frobenius_invariants(n1, n2, p)
    tp := p + 1 - n1;
    np := (n1^2 + n2)/2 - (p + 1)*n1 - p;
    return tp, np;
end function;

/*
    Compute the set of possible a_p(E) mod ell values that would be
    consistent with a congruence J(C)[ell] ~ E[ell].

    If such a congruence exists, then a_p(E)^2 + t_p * a_p(E) + n_p = 0 (mod ell)
    or a_p(E)^2 - t_p * a_p(E) + n_p = 0 (mod ell).

    Returns the set of roots (and their negatives) modulo ell.
*/
function compute_valid_ap_values_mod_ell(tp, np, ell)
    F_ell := GF(ell);
    // Roots of X^2 + t_p X + n_p = 0 (mod ell)
    roo := {r[1] : r in Roots(Polynomial([np, tp, 1]), F_ell)};
    // Include negatives (accounts for the other factor)
    roo := roo join {-r : r in roo};
    return roo;
end function;

//=============================================================================
// CANDIDATE FILTERING
//=============================================================================

/*
    Filter elliptic curve indices by Frobenius trace matching.

    For a hyperelliptic curve C, find which elliptic curves E (by index)
    have a_p(E) mod ell consistent with a mod-ell congruence to J(C).

    Returns: list of surviving indices (into Es/apEs arrays)
*/
function filter_indices_by_frobenius(C, discCmod, apEs, ell, primesbound)
    indsleft := [1..#apEs];
    F_ell := GF(ell);

    for p in PrimesUpTo(primesbound) do
        // Skip primes of bad reduction
        if discCmod mod p eq 0 then
            continue;
        end if;

        // Compute point counts
        n1, n2 := compute_point_counts(C, p);
        if n1 cmpeq false then
            continue;
        end if;

        // Compute Frobenius invariants
        tp, np := compute_frobenius_invariants(n1, n2, p);

        // Find valid a_p values mod ell
        valid_ap_mod_ell := compute_valid_ap_values_mod_ell(tp, np, ell);

        // Filter indices: keep only those with matching a_p mod ell
        indsleft := [i : i in indsleft |
                     F_ell!(apEs[i][p]) in valid_ap_mod_ell];

        if #indsleft eq 0 then
            break;
        end if;
    end for;

    return indsleft;
end function;

//=============================================================================
// MAIN SEARCH FUNCTIONS
//=============================================================================

/*
    Search for "good pairs" - specializations (a,b,c) of the parametric
    family f that yield hyperelliptic curves with potential mod-ell
    congruences to elliptic curves in our database.

    Returns: list of <indices, [a,b,c]> pairs
*/
function find_good_pairs(f, lowhtrationals, apEs, ell, primesthreshold,
                         primesbound, lowhtrationalsBound)
    goodpairs := [**];
    count := 0;

    for aa, bb, cc in lowhtrationals do
        // Respect the search bound
        if count gt lowhtrationalsBound then
            break;
        end if;
        count +:= 1;

        // Progress indicator
        if count mod 1000 eq 0 then
            print count;
        end if;

        // Specialize the parametric family
        g := Evaluate(f, [aa, bb, cc, t]);

        // Try to construct the hyperelliptic curve
        try
            C := HyperellipticCurve(g);
        catch e;
            continue;
        end try;

        // Check discriminant validity
        discCmod := discriminant_as_integer(C);
        if not is_discriminant_valid(discCmod, ell, primesthreshold) then
            continue;
        end if;

        // Filter elliptic curves by Frobenius matching
        indsleft := filter_indices_by_frobenius(C, discCmod, apEs, ell, primesbound);

        if #indsleft eq 0 then
            continue;
        end if;

        // Found a good pair!
        print goodpairs;
        Append(~goodpairs, <indsleft, [aa, bb, cc]>);
    end for;

    return goodpairs;
end function;

/*
    Filter good pairs by conductor exponent at ell.
    Keep only those with good reduction at ell (conductor exponent = 0).

    Returns: list of pairs with good reduction at ell
*/
function filter_by_conductor_at_ell(goodpairs, f, ell, prec)
    verygoodpairs := [**];
    Qell := pAdicField(ell, prec);

    for x in goodpairs do
        g := Evaluate(f, x[2] cat [t]);
        C := HyperellipticCurve(g);

        // Check conductor exponent at ell
        Cell := ChangeRing(C, Qell);
        nell := ConductorExponent(Cell);

        if nell ne 0 then
            continue;
        end if;

        Append(~verygoodpairs, x);
    end for;

    return verygoodpairs;
end function;

//=============================================================================
// CONGRUENCE VERIFICATION
//=============================================================================

/*
    Verify that a congruence E[ell] ~ J(C)[ell] holds by checking
    Frobenius traces up to a given prime bound.

    For a valid congruence, we need:
        a_p(E)^2 + t_p * a_p(E) + n_p = 0 (mod ell)  OR
        a_p(E)^2 - t_p * a_p(E) + n_p = 0 (mod ell)

    Returns: true if congruence is verified, false otherwise
*/
function verify_congruence(E, C, ell : primesbound := 1000)
    F_ell := GF(ell);

    // Compute discriminant for bad prime check
    discC := Discriminant(C);
    N := ZZ!(Numerator(discC) * Denominator(discC));

    for p in PrimesUpTo(primesbound) do
        // Skip bad primes
        if N mod p eq 0 then
            continue;
        end if;

        // Get a_p(E)
        ap := TraceOfFrobenius(E, p);

        // Compute point counts on C
        n1, n2 := compute_point_counts(C, p);
        if n1 cmpeq false then
            continue;
        end if;

        // Compute Frobenius invariants
        tp, np := compute_frobenius_invariants(n1, n2, p);

        // Check both possible congruence conditions
        val1 := F_ell!(ap^2 + tp*ap + np);
        val2 := F_ell!(ap^2 - tp*ap + np);

        if val1 ne 0 and val2 ne 0 then
            return false;
        end if;
    end for;

    return true;
end function;

//=============================================================================
// RANK BOUNDS AND FINAL VERIFICATION
//=============================================================================

/*
    For each very good pair, verify the congruence rigorously and
    compute Mordell-Weil rank bounds for the Jacobian.

    Returns: list of <coefficients, rank_bounds, verified_curves>
*/
function compute_verified_pairs_with_ranks(verygoodpairs, Es, f, ell, highprimesbound)
    verygoodpairs1 := [**];
    mwranks := [**];

    for x in verygoodpairs do
        g := Evaluate(f, x[2] cat [t]);
        C := HyperellipticCurveOfGenus(2, g);

        // Verify congruence for each candidate elliptic curve
        Es_for_C := [];
        for y in x[1] do
            E := Es[y];
            if verify_congruence(E, C, ell : primesbound := highprimesbound) then
                printf "Verified congruence for primes upto %o\n", highprimesbound;
                Append(~Es_for_C, aInvariants(E));
            end if;
        end for;

        // Compute rank bounds for the Jacobian
        J := Jacobian(C);
        rkbdlow, rkbdupp := RankBounds(J);
        rkbds := [rkbdlow, rkbdupp];
        printf "Rank Bounds for Jac(C): %o\n", rkbds;

        Append(~mwranks, rkbds);
        Append(~verygoodpairs1, <Coefficients(g), rkbds, Es_for_C>);
    end for;

    return verygoodpairs1, mwranks;
end function;

//=============================================================================
// TAMAGAWA NUMBER COMPUTATION
//=============================================================================

/*
    Compute Tamagawa numbers for each curve in the final list.
    Uses the Tamagawa package for regular models.
*/
procedure compute_tamagawa_numbers(verygoodpairs1)
    for x in verygoodpairs1 do
        C := HyperellipticCurveOfGenus(2, P1!(x[1]));
        discC := Discriminant(C);
        N := ZZ!(Numerator(discC) * Denominator(discC));

        try
            Tamagawa := [<p, TamagawaNumber(RegularModel(C, p))> : p in PrimeFactors(N)];
            printf "Tamagawa numbers for C:\n%o\n", Tamagawa;
        catch e;
            print e;
        end try;
    end for;
end procedure;

//=============================================================================
// MAIN EXECUTION
//=============================================================================

/*
    Main procedure that orchestrates the entire computation:
    1. Generate low-height rationals
    2. Load elliptic curves database
    3. Load parametric curve family
    4. Search for good pairs
    5. Filter by conductor
    6. Verify congruences and compute rank bounds
    7. Compute Tamagawa numbers
*/
procedure main()
    print "=== Searching for mod-ell congruences between genus 2 Jacobians and elliptic curves ===";
    printf "Parameters: ell = %o, D = %o\n", ELL, D_PARAM;

    // Step 1: Generate low-height rationals
    print "\n--- Generating low-height rationals ---";
    lowhtrationals := generate_low_height_rationals(HEIGHT_BOUND);
    printf "Generated %o low-height rationals\n", #lowhtrationals;
    lowhtrationals;

    // Step 2: Load elliptic curves
    print "\n--- Loading elliptic curves ---";
    filename := Sprintf("../data/ellipticcurvesinfo%o.txt", ELL);
    Es, apEs := load_elliptic_curves(filename, PRIMES_BOUND);
    printf "Loaded %o elliptic curves with Frobenius data\n", #Es;

    // Step 3: Load parametric family
    print "\n--- Loading curve family ---";
    f := load_curve_family(D_PARAM);

    // Step 4: Search for good pairs
    print "\n--- Searching for good pairs ---";
    goodpairs := find_good_pairs(f, lowhtrationals, apEs, ELL,
                                  PRIMES_THRESHOLD, PRIMES_BOUND,
                                  LOW_HT_RATIONALS_BOUND);
    printf "Found %o good pairs\n", #goodpairs;

    // Step 5: Filter by conductor at ell
    print "\n--- Filtering by conductor exponent at ell ---";
    verygoodpairs := filter_by_conductor_at_ell(goodpairs, f, ELL, PADIC_PRECISION);
    printf "After conductor filtering: %o pairs remain\n", #verygoodpairs;

    // Step 6: Verify congruences and compute rank bounds
    print "\n--- Verifying congruences and computing rank bounds ---";
    SetClassGroupBounds("GRH");
    verygoodpairs1, mwranks := compute_verified_pairs_with_ranks(
        verygoodpairs, Es, f, ELL, HIGH_PRIMES_BOUND);

    // Step 7: Compute Tamagawa numbers
    print "\n--- Computing Tamagawa numbers ---";
    Attach("../Tamagawa/Tamagawa_pkg2.m");
    compute_tamagawa_numbers(verygoodpairs1);

    print "\n=== Computation complete ===";
end procedure;

//=============================================================================
// RUN MAIN
//=============================================================================

main();

/*
    Sample output from a previous run:

    Verified congruence for primes upto 1000
    Rank Bounds for Jac(C): [ 4, 4 ]
    Verified congruence for primes upto 1000
    Rank Bounds for Jac(C): [ 0, 2 ]
    Verified congruence for primes upto 1000
    Rank Bounds for Jac(C): [ 2, 2 ]
    ...

    Sample good pairs found:
    [*
    <[ 1401 ], [ 0, 0, -5/2 ]>,
    <[ 2035 ], [ 0, 4, -1/3 ]>,
    <[ 1078 ], [ 0, 3, -9/8 ]>,
    <[ 4092, 4093 ], [ 0, 5, -1/4 ]>,
    <[ 310 ], [ 0, 5, -2/19 ]>,
    <[ 3191, 3192 ], [ 0, 2, 17/4 ]>,
    <[ 1401 ], [ 0, -3, 13/2 ]>,
    <[ 516 ], [ 0, 1/2, -7/8 ]>,
    <[ 947 ], [ 0, 9/2, -5/8 ]>,
    <[ 2369, 2370 ], [ 0, 9/4, -1/4 ]>
    *]
*/
