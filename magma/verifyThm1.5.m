/*
    verifyThm1.5.m

    Verifies Theorem 1.5 for the genus 2 curve
        -2y^2 = x^6 + 2x^4 + 12x^3 + 5x^2 + 6x + 1

    Run with:
        magma -b verifyThm1.5.m
*/

ZZ := Integers();
P<x> := PolynomialRing(Rationals());

// The curve: -2y^2 = x^6 + 2x^4 + 12x^3 + 5x^2 + 6x + 1
// Rewrite as y^2 = (-1/2)(x^6 + 2x^4 + 12x^3 + 5x^2 + 6x + 1)
f := -1/2 * (x^6 + 2*x^4 + 12*x^3 + 5*x^2 + 6*x + 1);

C := HyperellipticCurve(f);
print "Curve C:", C;

//=============================================================================
// 1. Good reduction at 5?
//=============================================================================
print "\n=== 1. Reduction at 5 ===";

// Method A: check discriminant
discC := Discriminant(C);
N := ZZ!(Numerator(discC) * Denominator(discC));
printf "Discriminant (as integer): %o\n", N;
printf "Valuation of discriminant at 5: %o\n", Valuation(N, 5);

// Method B: conductor exponent over Q_5
Q5 := pAdicField(5, 100);
C5 := ChangeRing(C, Q5);
n5 := ConductorExponent(C5);
printf "Conductor exponent at 5: %o\n", n5;

if n5 eq 0 then
    print "RESULT: C has GOOD reduction at 5.";
else
    print "RESULT: C has BAD reduction at 5.";
end if;

//=============================================================================
// 2. Rank of the Jacobian
//=============================================================================
print "\n=== 2. Mordell-Weil rank of Jac(C) ===";

J := Jacobian(C);
SetClassGroupBounds("GRH");
rklow, rkupp := RankBounds(J);
printf "Rank bounds: [%o, %o]\n", rklow, rkupp;

if rklow eq 4 and rkupp eq 4 then
    print "RESULT: Rank is exactly 4.";
elif rklow le 4 and rkupp ge 4 then
    printf "RESULT: Rank is between %o and %o (includes 4 but not proven exact).\n", rklow, rkupp;
else
    printf "RESULT: Rank is between %o and %o (does NOT include 4).\n", rklow, rkupp;
end if;

//=============================================================================
// 3. Torsion subgroup of the Jacobian
//=============================================================================
print "\n=== 3. Torsion subgroup of Jac(C) ===";

C_intmodel := IntegralModel(C);
J_intmodel := Jacobian(C_intmodel);
T, phi := TorsionSubgroup(J_intmodel);
printf "Torsion subgroup: %o\n", T;
printf "Torsion order: %o\n", #T;
printf "Invariants: %o\n", Invariants(T);

//=============================================================================
// 4. Tamagawa numbers
//=============================================================================
print "\n=== 4. Tamagawa numbers of Jac(C) ===";

Attach("../Tamagawa/Tamagawa_pkg2.m");

badprimes := PrimeFactors(N);
printf "Bad primes (from discriminant): %o\n", badprimes;

for p in badprimes do
    try
        R := RegularModel(C_intmodel, p);
        tam := TamagawaNumber(R);
        printf "  Tamagawa number at %o: %o\n", p, tam;
    catch e;
        printf "  Tamagawa number at %o: FAILED (%o)\n", p, e`Object;
    end try;
end for;

print "\n=== Done ===";
quit;
