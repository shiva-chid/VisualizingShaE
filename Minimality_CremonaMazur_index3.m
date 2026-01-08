intrinsic CMcurveForIndex3Torsor(L : SeqEnum) -> Crv
{given the coefficients of a smooth plane cubic C (an index 3 torsor of an elliptic curve) 
returns a Cremona-Mazur curve over an at most biquadratic extension whose Jacobian visualizes C}
    F := Universe(L);
    F := FieldOfFractions(F);
    P1<t> := PolynomialRing(F);
    PP2 := ProjectiveSpace(F,2);
    P<x,y,z> := CoordinateRing(PP2);
    mons := [x^3,y^3,z^3,x^2*y,x^2*z,y^2*x,y^2*z,z^2*x,z^2*y,x*y*z];
    f := &+[L[i]*mons[i] : i in [1..#mons]];

    nonappearingmons := [x^2*y,x^2*z,z^2*x,z^2*y];
    othercoeffs := [MonomialCoefficient(f,mon) : mon in nonappearingmons];
    if Set(othercoeffs) eq {0} then
        temp := MonomialCoefficient(f,z^3);
        fnew := -f/temp;
    else
        cubicpol := Evaluate(f,[1,0,t]);
        a, M, g := KummerElementAndTransformation(cubicpol);
        K := BaseRing(M);
        PK<x,y,z> := PolynomialRing(K,3);
        P1<T> := PolynomialRing(K);
        newz := M[1,1]*z+M[1,2]*x;
        newx := M[2,1]*z+M[2,2]*x;
        fnew := Evaluate(f,[newx,y,newz]);
        temp := MonomialCoefficient(fnew,z^3);
        fnew := -fnew/temp;
        assert MonomialCoefficient(fnew,x^3) eq a;

        fnew := Evaluate(fnew,x,x-MonomialCoefficient(fnew,x^2*y)/(3*MonomialCoefficient(fnew,x^3))*y);
    //    print MonomialCoefficient(fnew,x^2*y);
        fnew := Evaluate(fnew,z,z-MonomialCoefficient(fnew,z^2*y)/(3*MonomialCoefficient(fnew,z^3))*y);
    //    print MonomialCoefficient(fnew,x^2*y), MonomialCoefficient(fnew,z^2*y);
        a := MonomialCoefficient(fnew,x^3);
        b := MonomialCoefficient(fnew,y^3);
        m := MonomialCoefficient(fnew,x*y*z);
        nonappearingmons := [x^2*y,x^2*z,z^2*x,z^2*y];
        othercoeffs := [MonomialCoefficient(fnew,mon) : mon in nonappearingmons];
    //    print othercoeffs;
        assert Set(othercoeffs) eq {0};
        assert MonomialCoefficient(fnew,z^3) eq -1;
    end if;

    C := GenusOneModel(fnew);
/*
    E := Jacobian(C); // probably does some minimisation automatically, which we want to avoid.
    R<X,Y,Z> := CoordinateRing(AmbientSpace(E));
    fE := DefiningEquation(E);
*/
    a1,a2,a3,a4,a6 := Explode(aInvariants(C));
    R<X,Y,Z> := PolynomialRing(K,3);
    fE := Y^2 + a1*X*Y + a3*Y - (X^3 + a2*X^2 + a4*X + a6);
    roo := Roots(T^2+T+1); w := roo[1][1];
    g := Y - w^2*m*X - 3*(1-w)*a*b + (w-w^2)*m^3/9;

    I := ideal<R|fE,Z^3-g>;
    M := EliminationIdeal(I, {X, Z});
    BM := Basis(M);
    assert #BM eq 1;
    A2<X,Z> := AffineSpace(K,2);
    pol2 := Evaluate(BM[1],[X,0,Z]);
    C2 := ProjectiveClosure(Curve(A2,pol2));
    return C2;
end intrinsic;


intrinsic PrimeWitnessForMinimality(C :: Crv) -> RngIntElt
{returns a prime witnessing the simplicity of the factor complementary to the natural elliptic curve
inside the Jacobian of C}
    K := BaseRing(C);
    OK := RingOfIntegers(K);
    g := Genus(C);
    List := [ p : p in PrimesUpTo(1000) | IsTotallySplit(p,OK) ];
    for p in List do
        try
            Ps := PrimeIdealsOverPrime(K,p);
            k, redOK := ResidueClassField(Ps[1]);
//            Cp := BaseChange(C, redOK);
            z := K.1;
            redK := hom< K -> k | redOK(OK!z) >;
            Cp := BaseChange(C, redK);
            Facs := Factorization(Numerator(ZetaFunction(Cp)));
            if #Facs eq 2 and {Degree(fac[1]) : fac in Facs} eq {2,2*g-2} then return p; end if;
        catch e;
//            print e;
            continue;
        end try;
    end for;
end intrinsic;

