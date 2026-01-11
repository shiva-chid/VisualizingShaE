intrinsic CMcurveForIndex3Torsor(L :: SeqEnum : direct := true) -> Crv
{given the coefficients of a smooth plane cubic C (an index 3 torsor of an elliptic curve) 
returns a Cremona-Mazur curve over an at most biquadratic extension whose Jacobian visualizes C}
    F := Universe(L);
    F := FieldOfFractions(F);
    P1<t> := PolynomialRing(F);
    P<x,y,z> := PolynomialRing(F,3);
    mons := [x^3,y^3,z^3,x^2*y,x^2*z,y^2*x,y^2*z,z^2*x,z^2*y,x*y*z];
    f := &+[L[i]*mons[i] : i in [1..#mons]];
/*
    repeat randmat := Matrix(Rationals(),3,3,[Random([-10..10]) : i in [1..9]]); until Determinant(randmat) ne 0;
    f := Evaluate(f,Eltseq(Matrix(P,1,3,[x,y,z])*ChangeRing(randmat,P)));
*/

/*
    L := ChangeUniverse(L,Rationals());
    C1 := GenusOneModel(3,L);
//    RT := RandomTransformation(3 : Unimodular := true);
//    C1 := RT*C1;
    f := DefiningEquations(C1)[1];
*/

    nonappearingmons := [x^2*y,x^2*z,z^2*x,z^2*y];
    othercoeffs := [MonomialCoefficient(f,mon) : mon in nonappearingmons];
    if Set(othercoeffs) eq {0} then
        temp := MonomialCoefficient(f,z^3);
        fnew := -f/temp;
        K<w> := CyclotomicField(3);
        P1<T> := PolynomialRing(K);
        PK<x,y,z> := PolynomialRing(K,3);
        fnew := PK!fnew;
    else
        f := Evaluate(f,z,z-MonomialCoefficient(f,z^2*x)/(3*MonomialCoefficient(f,z^3))*x);
        f := f/MonomialCoefficient(f,z^3);
        cubicpol := Evaluate(f,[1,0,t]);
        if direct then
            M := TransformationDirect(cubicpol);
        else
            a, M, g := KummerElementAndTransformation(cubicpol);
        end if;
        K := BaseRing(M);
        P1<T> := PolynomialRing(K);
        roo := Roots(T^2+T+1); w := roo[1][1];
        PK<x,y,z> := PolynomialRing(K,3);
        newz := M[1,1]*z+M[1,2]*x;
        newx := M[2,1]*z+M[2,2]*x;
        fnew := Evaluate(f,[newx,y,newz]);
        temp := MonomialCoefficient(fnew,z^3);
        fnew := -fnew/temp;
    //    assert MonomialCoefficient(fnew,x^3) eq a;

        fnew := Evaluate(fnew,x,x-MonomialCoefficient(fnew,x^2*y)/(3*MonomialCoefficient(fnew,x^3))*y);
        fnew := Evaluate(fnew,z,z-MonomialCoefficient(fnew,z^2*y)/(3*MonomialCoefficient(fnew,z^3))*y);
        nonappearingmons := [x^2*y,x^2*z,z^2*x,z^2*y];
        othercoeffs := [MonomialCoefficient(fnew,mon) : mon in nonappearingmons];
    //    print othercoeffs;
        assert Set(othercoeffs) eq {0};
        assert MonomialCoefficient(fnew,z^3) eq -1;
    end if;
    a := MonomialCoefficient(fnew,x^3);
    b := MonomialCoefficient(fnew,y^3);
    m := MonomialCoefficient(fnew,x*y*z);
    C := GenusOneModel(fnew);
/*
    E := Jacobian(C); // probably does some minimisation automatically, which we want to avoid.
    R<X,Y,Z> := CoordinateRing(AmbientSpace(E));
    fE := DefiningEquation(E);
*/
    a1,a2,a3,a4,a6 := Explode(aInvariants(C));
    R<X,Y,Z> := PolynomialRing(K,3);
    fE := Y^2 + a1*X*Y + a3*Y - (X^3 + a2*X^2 + a4*X + a6);
    g := Y - w^2*m*X - 3*(1-w)*a*b + (w-w^2)*m^3/9;

    I := ideal<R|fE,Z^3-g>;
    M := EliminationIdeal(I, {X, Z});
    BM := Basis(M);
    assert #BM eq 1;
    A2<x,z> := AffineSpace(K,2);
    pol2 := Evaluate(BM[1],[x,0,z]);
    C2 := ProjectiveClosure(Curve(A2,pol2));
    return C2;
end intrinsic;

intrinsic CMcurveOverBaseFieldForIndex3Torsor(L :: SeqEnum : direct := true) -> Crv
{given the coefficients of a smooth plane cubic C (an index 3 torsor of an elliptic curve) 
returns a Cremona-Mazur curve over an at most biquadratic extension whose Jacobian visualizes C}
    F := Universe(L);
    F := FieldOfFractions(F);
    L := ChangeUniverse(L,F);
    C1 := GenusOneModel(3,L);
//    RT := RandomTransformation(3 : Unimodular := true);
//    C1 := RT*C1;
    f := DefiningEquations(C1)[1];
    P<x,y,z> := Parent(f);
    P1<t> := PolynomialRing(F);

    nonappearingmons := [x^2*y,x^2*z,z^2*x,z^2*y];
    othercoeffs := [MonomialCoefficient(f,mon) : mon in nonappearingmons];
    if Set(othercoeffs) eq {0} then
        temp := MonomialCoefficient(f,z^3);
        fnew := -f/temp;
        K<w> := CyclotomicField(3);
        P1<T> := PolynomialRing(K);
        PK<x,y,z> := PolynomialRing(K,3);
        fnew := PK!fnew;
    else
        f := Evaluate(f,z,z-MonomialCoefficient(f,z^2*x)/(3*MonomialCoefficient(f,z^3))*x);
        f := f/MonomialCoefficient(f,z^3);
        cubicpol := Evaluate(f,[1,0,t]);
        if direct then
            M := TransformationDirect(cubicpol);
        else
            a, M, g := KummerElementAndTransformation(cubicpol);
        end if;
        K := BaseRing(M);
        P1<T> := PolynomialRing(K);
        roo := Roots(T^2+T+1); w := roo[1][1];
        PK<x,y,z> := PolynomialRing(K,3);
        newz := M[1,1]*z+M[1,2]*x;
        newx := M[2,1]*z+M[2,2]*x;
        fnew := Evaluate(f,[newx,y,newz]);
        temp := MonomialCoefficient(fnew,z^3);
        fnew := -fnew/temp;
    //    assert MonomialCoefficient(fnew,x^3) eq a;

        fnew := Evaluate(fnew,x,x-MonomialCoefficient(fnew,x^2*y)/(3*MonomialCoefficient(fnew,x^3))*y);
        fnew := Evaluate(fnew,z,z-MonomialCoefficient(fnew,z^2*y)/(3*MonomialCoefficient(fnew,z^3))*y);
        nonappearingmons := [x^2*y,x^2*z,z^2*x,z^2*y];
        othercoeffs := [MonomialCoefficient(fnew,mon) : mon in nonappearingmons];
    //    print othercoeffs;
        assert Set(othercoeffs) eq {0};
        assert MonomialCoefficient(fnew,z^3) eq -1;
    end if;
    a := MonomialCoefficient(fnew,x^3);
    b := MonomialCoefficient(fnew,y^3);
    m := MonomialCoefficient(fnew,x*y*z);
    C := GenusOneModel(fnew);
/*
    E := Jacobian(C); // probably does some minimisation automatically, which we want to avoid.
    R<X,Y,Z> := CoordinateRing(AmbientSpace(E));
    fE := DefiningEquation(E);
*/
    a1,a2,a3,a4,a6 := Explode(aInvariants(C));
    R<X,Y,Z> := PolynomialRing(K,3);
    fE := Y^2 + a1*X*Y + a3*Y - (X^3 + a2*X^2 + a4*X + a6);
    g := Y - w^2*m*X - 3*(1-w)*a*b + (w-w^2)*m^3/9;
    coes, mons := CoefficientsAndMonomials(g);
    AutK, auts, tau := AutomorphismGroup(K);
    Ng := &*[&+[(tau(a))(coes[i])*mons[i] : i in [1..#coes]]: a in AutK];

    I := ideal<R|fE,Z^3-Ng>;
    M := EliminationIdeal(I, {X, Z});
    BM := Basis(M);
    assert #BM eq 1;
    A2<x,z> := AffineSpace(K,2);
    pol2 := Evaluate(BM[1],[x,0,z]);
    C2 := ProjectiveClosure(Curve(A2,pol2));
    C2 := ChangeRing(C2,F);
    return C2;
end intrinsic;

intrinsic PrimeWitnessForMinimality(C :: Crv : bound := 10^3) -> RngIntElt
{returns a prime witnessing the simplicity of the factor complementary to the natural elliptic curve
inside the Jacobian of C}
    K := BaseRing(C);
    OK := RingOfIntegers(K);
    g := Genus(C);
    for p in PrimesUpTo(bound) do
        if not IsTotallySplit(p,OK) then continue; end if;
        try
            printf "p = %o ", p;
            Ps := PrimeIdealsOverPrime(K,p);
            k, redOK := ResidueClassField(Ps[1]);
/*
            OK := RingOfIntegers(K);
            Ps := Factorisation(p*OK);
            k, redOK := ResidueClassField(Ps[1][1]);
*/
//            Cp := BaseChange(C, redOK);
            z := K.1;
            redK := hom< K -> k | redOK(OK!z) >;
            Cp := BaseChange(C, redK);
            Facs := Factorization(LPolynomial(Cp));
            printf "Factorisation type = %o\n", {*<Degree(fac[1]),fac[2]> : fac in Facs*};
            if #Facs eq 2 and {Degree(fac[1]) : fac in Facs} eq {2,2*g-2} then return p; end if;
        catch e;
//            print e;
            continue;
        end try;
    end for;
    return false;
end intrinsic;

