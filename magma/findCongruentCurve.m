/*
    This scripts helps find a congruent curve to one as output by
    the file `g2RMapproach.m`.
*/

// Define the elliptic curve
E := EllipticCurve([0, -1, 0, -169321, -28379327]);
printf "Factorisation of conductor of E: %o\n", Factorisation(Conductor(E));

ell := 5;
primesthreshold := 5*10^4;

// The curve needs to be in short Weierstrass form
E := WeierstrassModel(E);
a := Coefficients(E)[4];
b := Coefficients(E)[5];
J := jInvariant(E) / 1728;

// Get Fisher's Hesse polynomials
scal1 := -27; c4 := a/(-27);
scal2 := -54; c6 := b/(-54);
cE := [c4, c6];
DD, cc4, cc6 := HessePolynomials(5,1,cE);


// Get the alpha and beta polynomials
alpha, beta := RubinSilverbergPolynomials(5, J);

// Do some conversions of parent rings
R<t> := Parent(alpha);
beta := R!beta;

partnerCurves := [];
htbd := 5000;
lowhtrationals := Setseq({a/b : a in [-htbd..htbd], b in [1..htbd] | GCD(a,b) eq 1});
lowhtrationals := Sort(lowhtrationals, func<a,b|#Sprint(a) - #Sprint(b)>);
#lowhtrationals;
SetClassGroupBounds("GRH");
count := 0;
// for n in [1..100] do
for n in lowhtrationals do
//    printf "Trying n=%o\n", n;
    count +:= 1;
    if count mod 1000 eq 0 then printf "%o ", count; end if;

    // Get a nonsingular curve
    try
        // Using RubinSilverberg polynomials
        // F := EllipticCurve([a * Evaluate(alpha, n), b * Evaluate(beta, n)]);

        // Using Fisher's Hesse polynomials
        numdenn := [Numerator(n), Denominator(n)];
        F := EllipticCurve([scal1 * Evaluate(cc4, numdenn), scal2 * Evaluate(cc6, numdenn)]);
    catch e;
//        printf "Error in curve creation: %o\n", e;
        continue;
    end try;

    // Check if we have a reasonable curve
    F := MinimalModel(F);
    disc := Discriminant(F);
    disc := ZZ!(Numerator(disc)*Denominator(disc));
    if disc mod ell eq 0 then continue; end if;
    badprimes := [p : p in PrimesUpTo(primesthreshold) | disc mod p eq 0];
    if &*([1] cat [p^Valuation(disc, p) : p in badprimes]) ne Abs(disc) then continue; end if;
    printf "parameter = %o\nbad primes: %o\n", numdenn, badprimes;
//    printf "%o ", #Sprint(aInvariants(F));
//    if #Sprint(aInvariants(F)) gt 200 then continue; end if;


    printf "Computing bad primes...\n";
    time badPrimes := BadPrimes(F);

    if 5 in badPrimes then
        continue;
    end if;

    printf "Computing Tamagawa numbers...\n";
    time tamagawaNumbers := TamagawaNumbers(F);

    for x in tamagawaNumbers do
        if x mod 5 eq 0 then
            continue;
        end if;
    end for;

    // We need F to have rank at least 2
    printf "Computing rank...\n";
//    rkF := Rank(F);
//    time rkF, boo := Rank(F);
//    time anrkF := AnalyticRank(F : Precision := 2);
//    time rkF_low, rkF_upp := RankBounds(F : Effort := 2);
//    rkF := [rkF_low, rkF_upp];
    time rkF, gens, shainfo := MordellWeilShaInformation(F : RankOnly := true);
    printf "Mordell-Weil rank bounds for F = %o\n", rkF;
    if rkF[2] lt 2 then
        continue;
    end if;

    // If we have got to this point, we have a candidate curve
    print F;
    Append(~partnerCurves, F);
end for;
#partnerCurves;


// Found one!
F := partnerCurves[1]; aInvariants(F);
// [ 0, -1, 0, -2285384981, -42050896792563 ]
MordellWeilShaInformation(F);
/*
[ 2, 2 ]
[ (-9988484/361 : -5323181/6859 : 1), (-39939547/1444 : 564287421/54872 : 1) ]
[
<2, [ 0, 0 ]>
]
*/
BadPrimes(F);
// [ 2, 3, 29, 89 ]

apsE := [TraceOfFrobenius(E, p) : p in PrimesUpTo(10^3)];
apsF := [TraceOfFrobenius(F, p) : p in PrimesUpTo(10^3)];
{* GF(5)!(apsE[i] - apsF[i]) : i in [1..#apsE] *};
// {* 0^^167, 1 *}
{ NthPrime(i) : i in [1..#apsE] | GF(5)!(apsE[i] - apsF[i]) ne 0};
// { 89 }

