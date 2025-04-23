ZZ := Integers();
P<a,b,c,x> := PolynomialRing(Rationals(),4);
P1<t> := PolynomialRing(Rationals());

htbd := 20;
lowhtrationals := [0] cat Setseq({a/b : a in [1..htbd], b in [-htbd..htbd] | b ne 0 and GCD(a,b) eq 1});
lowhtrationals := Sort(lowhtrationals, func<a,b|#Sprint(a) - #Sprint(b)>);
#lowhtrationals;
lowhtrationals;

primesthreshold := 10^4;
primesbound := 10^2;

ell := 3; D := 13; E := EllipticCurve([1, 0, 0, -15663, -755809]);
// ell := 5; D := 5; E := EllipticCurve([1, -1, 0, -332311, -73733731]);
Dfile := Sprintf("%o.txt", D);
f := eval Read(Dfile);
 


Es := [E];
apEs := [AssociativeArray() : i in [1..#Es]];
for i in [1..#Es] do
    E := Es[i];
    for p in PrimesUpTo(primesbound) do
        apEs[i][p] := TraceOfFrobenius(E, p);
    end for;
end for;

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
    J := Jacobian(C);
    try
        Jell := BaseExtend(J,GF(ell));
    catch e;
        continue;
    end try;
    discC := Discriminant(C);
    discCmod := ZZ!(Numerator(discC)*Denominator(discC));
    badprimes := [p : p in PrimesUpTo(primesthreshold) | discCmod mod p eq 0];
    if &*([1] cat [p^Valuation(discCmod, p) : p in badprimes]) ne discCmod then continue; end if;
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

/*
[*
<[ 1 ], [ 0, 1/2, 3 ]>,
<[ 1 ], [ 0, 1/2, -7/6 ]>,
<[ 1 ], [ 0, 1/3, 1/6 ]>,
<[ 1 ], [ 0, 2/3, 2 ]>,
<[ 1 ], [ 0, 2/5, -1/15 ]>,
<[ 1 ], [ 0, 6/5, -2/5 ]>,
<[ 1 ], [ 0, 1/8, 1/6 ]>
*]
*/

///////////////////////////////////////////////

SetClassGroupBounds("GRH");

for x in goodpairs do
    g := Evaluate(f, x[2] cat [t]); 
    C := HyperellipticCurve(g);
    J := Jacobian(C);
    try
        Jell := BaseExtend(J,GF(ell));
    catch e;
        continue;
    end try;
//    RankBounds(J);
    print Factorisation(Conductor(C));
end for;
