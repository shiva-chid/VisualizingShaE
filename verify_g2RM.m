/*
    This is a verification script to test the logic. run with

    magma -bi verify_g2RM.m
*/

ZZ := Integers();
P<a,b,c,x> := PolynomialRing(Rationals(),4);
P1<t> := PolynomialRing(Rationals());

// According to Fisher, we can find an A, P, Q within this set
lowhtrationals := [41/4, 1/5, 0, 4, 6, 1/14];

primesthreshold := 10^4;
primesbound := 10^2;

ell := 7; E := EllipticCurve([1, 1, 0, -38582789, -87778717747]);

Es := [E];
apEs := [AssociativeArray() : i in [1..#Es]];
for i in [1..#Es] do
    E := Es[i];
    for p in PrimesUpTo(primesbound) do
        apEs[i][p] := TraceOfFrobenius(E, p);
    end for;
end for;

function RMby2(A,P,Q)
    R := 4*P;
    B := (Q*(P*A - Q) + 4*P^2 + 1)/P^2;
    C := 4*(P*A - Q) / P;
    cubic := t^3 + A*t^2 + B*t + C;
    Spl<b> := SplittingField(cubic);

    roots := Roots(cubic, Spl);

    g := 1;
    Q1<T> := PolynomialRing(Spl);
    for alpha in roots do
        new_term := T^2 - alpha[1]*T + P*alpha[1]^2 + Q*alpha[1] + R;
        g := g * new_term;
    end for;
    return P1!g;
end function;

function get_goodpairs(apEs,ell)
    goodpairs := [**];
    count := 0;
    for aa, bb, cc in lowhtrationals do
        count +:= 1;
        if count mod 1000 eq 0 then print count; end if;
        //g := Evaluate(f, [aa,bb,cc,t]);
        if bb eq 0 then continue; end if;
        g := RMby2(aa, bb, cc);
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
    return goodpairs;
end function;

ans := get_goodpairs(apEs,ell);

if #ans eq 0 then
    print "FAIL";
else
    print "PASS";
end if;