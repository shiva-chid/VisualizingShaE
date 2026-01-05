/*
    findCongruentCurvev2.m

    This script finds an elliptic curve F that is ell-congruent to a given
    elliptic curve E (i.e., E[ell] ~ F[ell] as Galois modules).

    The approach uses Fisher's Hesse polynomials to parametrize curves with
    the same mod-ell representation. We search through low-height rational
    parameters to find curves F satisfying:
    - Good reduction at ell
    - All bad primes below a threshold
    - Tamagawa numbers coprime to ell
    - Mordell-Weil rank at least 2

    This is used in conjunction with g2RMapproach.m to find curves for
    visualizing elements of Sha.
*/

//=============================================================================
// GLOBAL SETUP
//=============================================================================

ZZ := Integers();

//=============================================================================
// CONFIGURATION PARAMETERS
//=============================================================================

// The prime ell for the mod-ell representation
ELL := 5;

// Prime threshold for checking bad primes of discriminant
PRIMES_THRESHOLD := 5*10^4;

// Height bound for searching rational parameters
HEIGHT_BOUND := 5000;

// Minimum required rank for the partner curve
MIN_RANK := 2;

// Prime bound for Frobenius trace verification
FROBENIUS_VERIFICATION_BOUND := 10^3;

//=============================================================================
// INPUT CURVE
//=============================================================================

// The elliptic curve E for which we seek a congruent partner
// (This curve typically comes from g2RMapproach.m output)
INPUT_CURVE_AINVS := [0, -1, 0, -169321, -28379327];

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

//=============================================================================
// CURVE PREPARATION
//=============================================================================

/*
    Prepare the input curve for the congruence search.
    Converts to short Weierstrass form and extracts coefficients.

    Returns:
        E    - the curve in Weierstrass form
        a, b - the Weierstrass coefficients (y^2 = x^3 + ax + b)
        J    - the normalized j-invariant (j/1728)
*/
function prepare_input_curve(ainvs)
    E := EllipticCurve(ainvs);
    printf "Factorisation of conductor of E: %o\n", Factorisation(Conductor(E));

    // Convert to short Weierstrass form: y^2 = x^3 + ax + b
    E := WeierstrassModel(E);
    a := Coefficients(E)[4];
    b := Coefficients(E)[5];
    J := jInvariant(E) / 1728;

    return E, a, b, J;
end function;

/*
    Compute Fisher's Hesse polynomial data for parametrizing
    ell-congruent curves.

    The Hesse polynomials give a rational parametrization of curves
    with isomorphic mod-ell Galois representations.

    Returns:
        scal1, scal2 - scaling factors for c4, c6
        cc4, cc6     - the Hesse polynomials
*/
function compute_hesse_data(a, b, ell)
    // Standard scaling for Hesse polynomials
    scal1 := -27;
    scal2 := -54;

    c4 := a / scal1;
    c6 := b / scal2;
    cE := [c4, c6];

    // Get Fisher's Hesse polynomials for the given ell
    DD, cc4, cc6 := HessePolynomials(ell, 1, cE);

    return scal1, scal2, cc4, cc6;
end function;

/*
    Compute Rubin-Silverberg polynomials (alternative parametrization).
    These give another way to parametrize ell-congruent curves.

    Returns:
        alpha, beta - the Rubin-Silverberg polynomials in R<t>
        R           - the parent polynomial ring
*/
function compute_rubin_silverberg_data(J, ell)
    alpha, beta := RubinSilverbergPolynomials(ell, J);
    R<t> := Parent(alpha);
    beta := R!beta;
    return alpha, beta, R;
end function;

//=============================================================================
// DISCRIMINANT AND VALIDITY CHECKS
//=============================================================================

/*
    Compute the discriminant of an elliptic curve as an integer.
*/
function discriminant_as_integer(F)
    disc := Discriminant(F);
    return ZZ!(Numerator(disc) * Denominator(disc));
end function;

/*
    Check if the discriminant satisfies our requirements:
    - Not divisible by ell
    - All prime factors are below primesthreshold

    Returns: true if valid, false otherwise
*/
function is_discriminant_valid(disc, ell, primesthreshold)
    // Must have good reduction at ell
    if disc mod ell eq 0 then
        return false;
    end if;

    // Check that all bad primes are below threshold
    badprimes := [p : p in PrimesUpTo(primesthreshold) | disc mod p eq 0];
    product_of_prime_powers := &*([1] cat [p^Valuation(disc, p) : p in badprimes]);

    if product_of_prime_powers ne Abs(disc) then
        return false;
    end if;

    return badprimes;
end function;

//=============================================================================
// PARTNER CURVE CONSTRUCTION
//=============================================================================

/*
    Construct a partner curve using Fisher's Hesse polynomials.
    Given a rational parameter n, evaluates the Hesse polynomials
    to get a curve with the same mod-ell representation.

    Returns: the constructed elliptic curve, or false if construction fails
*/
function construct_partner_curve_hesse(n, scal1, scal2, cc4, cc6)
    try
        numdenn := [Numerator(n), Denominator(n)];
        F := EllipticCurve([scal1 * Evaluate(cc4, numdenn),
                            scal2 * Evaluate(cc6, numdenn)]);
        return F, numdenn;
    catch e;
        return false, false;
    end try;
end function;

/*
    Alternative: Construct partner curve using Rubin-Silverberg polynomials.
    (Currently commented out in main search, but available for use)
*/
function construct_partner_curve_rubin_silverberg(n, a, b, alpha, beta)
    try
        F := EllipticCurve([a * Evaluate(alpha, n), b * Evaluate(beta, n)]);
        return F;
    catch e;
        return false;
    end try;
end function;

//=============================================================================
// FILTERING FUNCTIONS
//=============================================================================

/*
    Check if a curve passes the bad primes filter.
    The curve must have good reduction at ell.

    Returns: true if passes, false otherwise
*/
function passes_bad_primes_filter(F, ell)
    badPrimes := BadPrimes(F);
    return not (ell in badPrimes);
end function;

/*
    Check if a curve passes the Tamagawa number filter.
    All Tamagawa numbers must be coprime to ell.

    Returns: true if passes, false otherwise
*/
function passes_tamagawa_filter(F, ell)
    tamagawaNumbers := TamagawaNumbers(F);
    for x in tamagawaNumbers do
        if x mod ell eq 0 then
            return false;
        end if;
    end for;
    return true;
end function;

/*
    Check if a curve passes the rank filter.
    The Mordell-Weil rank must be at least min_rank.

    Returns: rank_bounds if passes, false otherwise
*/
function passes_rank_filter(F, min_rank)
    rkF, gens, shainfo := MordellWeilShaInformation(F : RankOnly := true);
    printf "Mordell-Weil rank bounds for F = %o\n", rkF;

    if rkF[2] lt min_rank then
        return false;
    end if;

    return rkF;
end function;

//=============================================================================
// MAIN SEARCH
//=============================================================================

/*
    Search for partner curves that are ell-congruent to E.

    Iterates through low-height rational parameters, constructs candidate
    curves using Hesse polynomials, and filters by various criteria.

    Returns: list of partner curves satisfying all criteria
*/
function find_partner_curves(scal1, scal2, cc4, cc6, ell,
                              primesthreshold, lowhtrationals, min_rank)
    partnerCurves := [];
    count := 0;

    for n in lowhtrationals do
        count +:= 1;
        if count mod 1000 eq 0 then
            printf "%o ", count;
        end if;

        // Construct candidate curve using Hesse polynomials
        F, numdenn := construct_partner_curve_hesse(n, scal1, scal2, cc4, cc6);
        if F cmpeq false then
            continue;
        end if;

        // Get minimal model
        F := MinimalModel(F);

        // Check discriminant validity
        disc := discriminant_as_integer(F);
        badprimes := is_discriminant_valid(disc, ell, primesthreshold);
        if badprimes cmpeq false then
            continue;
        end if;
        printf "parameter = %o\nbad primes: %o\n", numdenn, badprimes;

        // Check bad primes at ell
        printf "Computing bad primes...\n";
        time passes_bp := passes_bad_primes_filter(F, ell);
        if not passes_bp then
            continue;
        end if;

        // Check Tamagawa numbers
        printf "Computing Tamagawa numbers...\n";
        time passes_tam := passes_tamagawa_filter(F, ell);
        if not passes_tam then
            continue;
        end if;

        // Check rank
        printf "Computing rank...\n";
        time rkF := passes_rank_filter(F, min_rank);
        if rkF cmpeq false then
            continue;
        end if;

        // If we reach here, we have a valid partner curve
        print F;
        Append(~partnerCurves, F);
    end for;

    return partnerCurves;
end function;

//=============================================================================
// CONGRUENCE VERIFICATION
//=============================================================================

/*
    Verify that two elliptic curves E and F are ell-congruent by
    comparing their Frobenius traces modulo ell.

    For congruent curves, a_p(E) = a_p(F) (mod ell) for all primes p
    of good reduction.

    Returns:
        diff_multiset - multiset of (a_p(E) - a_p(F)) mod ell values
        bad_primes    - set of primes where traces differ mod ell
*/
function verify_frobenius_congruence(E, F, ell, primesbound)
    primes := PrimesUpTo(primesbound);
    apsE := [TraceOfFrobenius(E, p) : p in primes];
    apsF := [TraceOfFrobenius(F, p) : p in primes];

    // Multiset of differences mod ell (should be mostly 0s)
    diff_multiset := {* GF(ell)!(apsE[i] - apsF[i]) : i in [1..#apsE] *};

    // Primes where traces differ mod ell
    bad_primes := { primes[i] : i in [1..#apsE] | GF(ell)!(apsE[i] - apsF[i]) ne 0 };

    return diff_multiset, bad_primes;
end function;

//=============================================================================
// MAIN EXECUTION
//=============================================================================

/*
    Main procedure that orchestrates the search for congruent curves:
    1. Prepare the input curve
    2. Compute Hesse polynomial data
    3. Generate low-height rationals
    4. Search for partner curves
    5. Analyze and verify the results
*/
procedure main()
    print "=== Searching for ell-congruent partner curves ===";
    printf "Parameters: ell = %o\n", ELL;

    // Step 1: Prepare the input curve
    print "\n--- Preparing input curve ---";
    E, a, b, J := prepare_input_curve(INPUT_CURVE_AINVS);

    // Step 2: Compute Hesse polynomial data
    print "\n--- Computing Hesse polynomials ---";
    scal1, scal2, cc4, cc6 := compute_hesse_data(a, b, ELL);

    // Also compute Rubin-Silverberg data (for reference/alternative use)
    alpha, beta, R := compute_rubin_silverberg_data(J, ELL);

    // Step 3: Generate low-height rationals
    print "\n--- Generating search parameters ---";
    lowhtrationals := generate_low_height_rationals(HEIGHT_BOUND);
    printf "Generated %o low-height rationals\n", #lowhtrationals;

    // Step 4: Search for partner curves
    print "\n--- Searching for partner curves ---";
    SetClassGroupBounds("GRH");
    partnerCurves := find_partner_curves(scal1, scal2, cc4, cc6, ELL,
                                          PRIMES_THRESHOLD, lowhtrationals,
                                          MIN_RANK);
    printf "\nFound %o partner curves\n", #partnerCurves;

    // Step 5: Analyze results
    if #partnerCurves gt 0 then
        print "\n--- Analyzing first partner curve ---";
        F := partnerCurves[1];
        printf "a-invariants: %o\n", aInvariants(F);

        printf "\nMordell-Weil-Sha information:\n";
        MordellWeilShaInformation(F);

        printf "\nBad primes: %o\n", BadPrimes(F);

        // Verify the congruence
        print "\n--- Verifying Frobenius congruence ---";
        diff_multiset, bad_primes := verify_frobenius_congruence(
            E, F, ELL, FROBENIUS_VERIFICATION_BOUND);
        printf "Difference multiset mod %o: %o\n", ELL, diff_multiset;
        printf "Primes with non-zero difference: %o\n", bad_primes;
    end if;

    print "\n=== Search complete ===";
end procedure;

//=============================================================================
// RUN MAIN
//=============================================================================

main();

/*
    Sample output from a previous run:

    Found partner curve F with:
    a-invariants: [ 0, -1, 0, -2285384981, -42050896792563 ]

    MordellWeilShaInformation:
    [ 2, 2 ]
    [ (-9988484/361 : -5323181/6859 : 1), (-39939547/1444 : 564287421/54872 : 1) ]
    [
    <2, [ 0, 0 ]>
    ]

    BadPrimes: [ 2, 3, 29, 89 ]

    Frobenius trace verification:
    Difference multiset mod 5: {* 0^^167, 1 *}
    Primes with non-zero difference: { 89 }

    (The difference at p=89 is expected since 89 is a prime of bad reduction)
*/
