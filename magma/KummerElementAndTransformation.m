intrinsic KummerElementAndTransformation(f :: RngUPolElt) -> FldNumElt, AlgMatElt
{given an irreducible degree 3 polynomial over Q defining a cubic extension F|Q,
returns an element a in the field K = Q(sqrt(discriminant of f), zeta_3) such that
the normal closure of F over Q(zeta_3) is K(a^(1/3)).
Also, returns a 2x2 matrix gamma over K, such that the corresponding Mobius transformation
sends the 3 cuberoots of a to the 3 roots of f.
The third return value is a cubic polynomial g(x) over K, which is a scalar multiple of
the polynomial x^3-a, and such that f(gamma dot x)*(cx+d)^3 = g(x).}
    P<x> := PolynomialRing(Rationals());
    discf := Discriminant(f);
    D, sqpart := Squarefree(Numerator(discf)*Denominator(discf));
    Qzeta<zeta3> := CyclotomicField(3);
    if D in {1,-3} then
        K := Qzeta;
    else
        K, roots := SplittingField([x^2+x+1,x^2-D]);
        zeta3 := roots[1][1];
    end if;
    fK := ChangeRing(f,K);
    L<alp> := ext<K|fK>;
    A, auts, tau := AutomorphismGroup(L);
//    assert #A eq 3;
//    AK := FixedGroup(L,K); #AK;
    assert exists(sigma){s : s in A | Order(s) ne 1};

    cuberootofa := &+[zeta3^i*(tau(sigma^i))(alp) : i in [0..2]]^-1;
    a := cuberootofa^3;
    assert a in K;

    roots_f := [(tau(sigma^i))(alp) : i in [0..2]];
    roots_g := [cuberootofa*zeta3^i : i in [0..2]];

    fL := ChangeRing(f,L);
    f1, M1 := InverseFractionalLinearTransformation(fL,[[x,1] : x in roots_g] : monic := true);
//    f1;
    f2, M2 := FractionalLinearTransformation(f1,[[x,1] : x in roots_f] : monic := true);
//    f2;

//    return a, cuberootofa, zeta3, M1, M2;

    fK := ChangeRing(f,K);
    M := M2*M1;
    Melts := Eltseq(M);
    assert exists(nonzeroentry){x : x in Melts | x ne 0};
    MK := Matrix(K,2,2,[K!(x/nonzeroentry) : x in Melts]);
    g := ApplyFractionalLinearTransformation(fK,MK : maychangedegree := false);
    assert Degree(g) eq 3 and {Coefficient(g,i) : i in [1,2]} eq {0};
    g0 := Coefficient(g,0);
    g3 := Coefficient(g,3);
    assert g0 eq -a*g3;

//    print D, sqpart;
    return a, MK, g;
end intrinsic;


/*
// initial debugging

AttachSpec("magma/spec");
f := x^3+x+1;
f := x^3+3*x+1;
f := x^3-12*x+8;
a, cbrta, z3, M1, M2 := KummerElementAndTransformation(f);
L := BaseRing(M2);
fL := ChangeRing(f,L);

g := ApplyFractionalLinearTransformation(fL,M1*M2); g;
a in [x[1]^3 : x in Roots(g)];
#{x[1]^3 : x in Roots(g)};
{cbrta*z3^i : i in [0..2]} subset {x[1] : x in Roots(g)};

g := ApplyFractionalLinearTransformation(fL,(M1*M2)^-1); g;
a in [x[1]^3 : x in Roots(g)];
#{x[1]^3 : x in Roots(g)};
{cbrta*z3^i : i in [0..2]} subset {x[1] : x in Roots(g)};
*/

intrinsic TransformationDirect(f :: RngUPolElt) -> FldNumElt, AlgMatElt
{given an irreducible degree 3 polynomial f = x^3 + a1 x + a0 over Q defining a cubic extension F|Q,
returns a 2x2 matrix gamma over K = Q(sqrt(discriminant of f), zeta_3), such that the
corresponding Mobius transformation such that f(gamma dot x)*(cx+d)^3 is a scalar multiple
of x^3 - a for a suitable a in K.}
    P<x> := PolynomialRing(Rationals());
    assert Degree(f) eq 3 and Coefficient(f,3) eq 1 and Coefficient(f,2) eq 0;
    a0 := Coefficient(f,0); a1 := Coefficient(f,1);
    discf := Discriminant(f);
    K, roots := SplittingField([x^2+x+1,x^2+3*discf] : Abs := true, Opt := true);
/*
    pol := DefiningPolynomial(K);
    monicpol := pol/LeadingCoefficient(pol);
    K := NumberField(monicpol);
*/
    zeta3 := K!(Eltseq(roots[1][1]));
    sqrtminus3D := K!(Eltseq(roots[2][1]));
    alpha := (9*a0+sqrtminus3D)/(2*a1^2);
    n := -a1;
    M := Matrix(K,2,2,[[n,n*alpha/3],[n*alpha,1]]);
    return M;
end intrinsic;