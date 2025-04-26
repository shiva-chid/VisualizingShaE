/*
   Test comment, update later. Describe what this 
   script is doing.
*/

ZZ := Integers();
P<a,b,c,x> := PolynomialRing(Rationals(),4);
P1<t> := PolynomialRing(Rationals());

htbd := 20;
lowhtrationals := Setseq({a/b : a in [-htbd..htbd], b in [1..htbd] | GCD(a,b) eq 1});
lowhtrationals := Sort(lowhtrationals, func<a,b|#Sprint(a) - #Sprint(b)>);
#lowhtrationals;
lowhtrationals;

primesthreshold := 10^4;
primesbound := 10^2;


// ell := 3; D := 13; E := EllipticCurve([1, 0, 0, -15663, -755809]);
ell := 5; D := 5; E := EllipticCurve([1, -1, 0, -332311, -73733731]);
Es := [E];

filename := Sprintf("ellipticcurvesinfo%o.txt",ell);
fil := Open(filename, "r");
s := Gets(fil);
Es := [];
apEs := [];
while not IsEof(s) do
    s_split := Split(s, ":");
    ainvs := [StringToInteger(x) : x in Split(s_split[1], "[], ")];
    cond := StringToInteger(s_split[2]);
    label := s_split[3];
    aps := [StringToInteger(x) : x in Split(s_split[4], "[], ")];
    apsdict := AssociativeArray();
    primes := PrimesUpTo(primesbound);
    for i := 1 to #aps do
        apsdict[primes[i]] := aps[i];
    end for;
    Append(~Es, EllipticCurve(ainvs));
    Append(~apEs,apsdict);
    s := Gets(fil);
end while;
#Es, #apEs;

/*
apEs := [AssociativeArray() : i in [1..#Es]];
for i in [1..#Es] do
    E := Es[i];
    for p in PrimesUpTo(primesbound) do
        apEs[i][p] := TraceOfFrobenius(E, p);
    end for;
end for;
*/

Dfile := Sprintf("%o.txt", D);
f := eval Read(Dfile);
goodpairs := [**];
count := 0;
for aa, bb, cc in lowhtrationals do
    count +:= 1;
    if count mod 1000 eq 0 then print count; end if;
    g := Evaluate(f, [aa,bb,cc,t]);
    try
        C := HyperellipticCurve(g);
    catch e;
        continue;
    end try;
    discC := Discriminant(C);
    discCmod := ZZ!(Numerator(discC)*Denominator(discC));
    if discCmod mod ell eq 0 then continue; end if;
    badprimes := [p : p in PrimesUpTo(primesthreshold) | discCmod mod p eq 0];
    if &*([1] cat [p^Valuation(discCmod, p) : p in badprimes]) ne Abs(discCmod) then continue; end if;
    indsleft := [1..#apEs];
    for p in PrimesUpTo(primesbound) do
        if discCmod mod p eq 0 then continue; end if;
        try
            Cp := ChangeRing(C, GF(p)); Cp2 := ChangeRing(C, GF(p,2));
        catch e;
            continue;
//            print e, C, p, badprimes;
//            break aa;
        end try;
        n1 := #Cp;
        n2 := #Cp2;
        tp := p+1-n1;
        np := (n1^2+n2)/2-(p+1)*n1-p;
        roo := {r[1] : r in Roots(Polynomial([np,tp,1]), GF(ell))};
        roo := roo join {-r : r in roo};
        indsleft := [i : i in indsleft | apEmodell in roo where apEmodell is GF(ell)!(apEs[i][p])];
        if #indsleft eq 0 then break; end if;
    end for;
    if #indsleft eq 0 then continue; end if;
    print goodpairs;
    Append(~goodpairs, <indsleft, [aa,bb,cc]>);
end for;

///////////////////////////////////////////////

verygoodpairs := [**];
prec := 100;
for x in goodpairs do
    g := Evaluate(f, x[2] cat [t]); 
    C := HyperellipticCurve(g);
/*
    J := Jacobian(C);
    try
        Jell := BaseExtend(J,GF(ell));
    catch e;
        continue;
    end try;
*/
/*
    cond := Conductor(C);
    if cond mod ell eq 0 then continue; end if;
    print Factorisation(cond);
*/
    Qell := pAdicField(ell,prec);
    Cell := ChangeRing(C, Qell);
    nell := ConductorExponent(Cell);
    if nell ne 0 then continue; end if;
    Append(~verygoodpairs, x);
end for;

/*
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

function verify_congruence(E,C,ell : primesbound := 1000);
    F_ell := GF(ell);
//    N := 2*Conductor(E)*Conductor(C);
    discC := Discriminant(C);
    N := ZZ!(Numerator(discC)*Denominator(discC));
    for p in PrimesUpTo(primesbound) do
        if N mod p eq 0 then continue; end if;
        ap := TraceOfFrobenius(E, p);
        try
            Cp := ChangeRing(C, GF(p)); Cp2 := ChangeRing(C, GF(p,2));
        catch e;
            continue;
        end try;
        n1 := #Cp;
        n2 := #Cp2;
        tp := p+1-n1;
        np := (n1^2+n2)/2-(p+1)*n1-p;
        val1 := F_ell!(ap^2+tp*ap+np);
        val2 := F_ell!(ap^2-tp*ap+np);
        if val1 ne 0 and val2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

SetClassGroupBounds("GRH");

verygoodpairs1 := [**];
mwranks := [**];
highprimesbound := 10^3;
for x in verygoodpairs do
    g := Evaluate(f, x[2] cat [t]); 
    C := HyperellipticCurveOfGenus(2,g);
    Es_for_C := [];
    for y in x[1] do
        E := Es[y];
        if verify_congruence(E,C,ell : primesbound := highprimesbound) then
            printf "Verified congruence for primes upto %o\n", highprimesbound;
            Append(~Es_for_C, aInvariants(E));
        end if;
    end for;
    J := Jacobian(C);
    rkbdlow, rkbdupp := RankBounds(J);
    rkbds := [rkbdlow, rkbdupp];
    printf "Rank Bounds for Jac(C): %o\n", rkbds;
    Append(~mwranks, rkbds);
    Append(~verygoodpairs1, <Coefficients(g), rkbds, Es_for_C>);
end for;
/*
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 4, 4 ]
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 0, 2 ]
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 2, 2 ]
Verified congruence for primes upto 1000
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 2, 2 ]
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 0, 2 ]
Verified congruence for primes upto 1000
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 4, 4 ]
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 4, 4 ]
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 2, 4 ]
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 2, 2 ]
Verified congruence for primes upto 1000
Verified congruence for primes upto 1000
Rank Bounds for Jac(C): [ 0, 2 ]
*/

Attach("Tamagawa/Tamagawa_pkg2.m");
for x in verygoodpairs1 do
    C := HyperellipticCurveOfGenus(2,P1!(x[1]));
    discC := Discriminant(C);
    N := ZZ!(Numerator(discC)*Denominator(discC));
    try
        Tamagawa := [<p,TamagawaNumber(RegularModel(C,p))> : p in PrimeFactors(N)];
        printf "Tamagawa numbers for C:\n%o\n", Tamagawa;
    catch e;
        print e;
    end try;
end for;
