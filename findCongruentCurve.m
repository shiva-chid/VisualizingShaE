/*
    This scripts helps find a congruent curve to one as output by
    the file `g2RMapproach.m`.
*/

// Define the elliptic curve
E := EllipticCurve([0, -1, 0, -169321, -28379327]);

// The curve needs to be in short Weierstrass form
E := WeierstrassModel(E);
a := Coefficients(E)[4];
b := Coefficients(E)[5];
J := jInvariant(E) / 1728;

// Get the alpha and beta polynomials
alpha, beta := RubinSilverbergPolynomials(5, J);

// Do some conversions of parent rings
R<t> := Parent(alpha);
beta := R!beta;

n := 4;
F := EllipticCurve([a * Evaluate(alpha, n), b * Evaluate(beta, n)]);

partnerCurves := [];
for n in [1..100] do
    // Get a nonsingular curve
    try
        F := EllipticCurve([a * Evaluate(alpha, n), b * Evaluate(beta, n)]);
    catch e;
        continue;
    end try;

    badPrimes := BadPrimes(F);

    if 5 in badPrimes then
        continue;
    end if;

    tamagawaNumbers := TamagawaNumbers(F);

    for x in tamagawaNumbers do
        if x mod 5 eq 0 then
            continue;
        end if;
    end for;

    // We need F to have rank at least 2
    if Rank(F) lt 2 then
        continue;
    end if;

    // If we have got to this point, we have a candidate curve
    Append(~partnerCurves, F);
end for;