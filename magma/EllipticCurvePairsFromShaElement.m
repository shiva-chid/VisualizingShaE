intrinsic EllipticCurvePairsFromShaElement(abcde :: SeqEnum : direct := false)-> EllCrv, EllCrv, CrvHyp
{Given the coefficients [a,b,c,d,e] of a genus one model of degree 2 (that is not itself
isomorphic to an elliptic curve), return two elliptic curves E1, E2 such that E1 is the
Jacobian of the genus one model and E2 is an elliptic curve whose Mordell-Weil group explains
the order-2 element of H^1(Q,E1) corresponding to the given genus one model.}
    P<t> := PolynomialRing(Universe(abcde));
    a,b,c,d,e := Explode(abcde);
    f := P![-4*a*c*e + b^2*e + a*d^2, -(4*a*e-b*d), c, 1];
    try
        E1 := EllipticCurve(f);
    catch e;
        printf "Singular model\n";
//        return false, false;
        return false, false, false;
    end try;

    Delta := b^2-4*a*c;
//    f2 := Evaluate(f,(Delta-t^2)/(4*a));

//    f2 := Evaluate(f,(Delta-(t^2)/(-a))/(4*a));

    f2 := Evaluate(f,(Delta-(t^2)/a)/(4*a));

//    f2 := Evaluate(f,t^2-(c-3*b^2/(8*a)));

//    f2 := Evaluate(f,t^2/(-a)-(c-3*b^2/(8*a)));
//    print f2;
    try
        genus2curve := HyperellipticCurveOfGenus(2,f2);
    catch e;
        printf "Singular genus 2 curve\n";
//        return false, false;
        return false, false, false;
    end try;

    if direct then
        p := -64*a^3*e + 16*a^2*b*d + 16*a^2*c^2 - 16*a*b^2*c + 3*b^4;
        q := (-8*a*c + 3*b^2)*(8*a^2*d - 4*a*b*c + b^3)^2;
        r := (8*a^2*d - 4*a*b*c + b^3)^4;
        try
            E2 := EllipticCurve([0,p,0,q,r]);
            E2 := QuadraticTwist(E2,-a);
        catch e;
            printf "Singular explaining elliptic curve\n";
//            return false, false;
            return false, false, false;
        end try;
//        return E1, E2;
        return E1, E2, genus2curve;
    end if;

    coeffs := Coefficients(f2);
    assert #coeffs eq 7 and coeffs[[2,4,6]] eq [0,0,0];
    coeffsofxsq := [coeffs[i] : i in [1,3,5,7]];
    fE1 := P!coeffsofxsq;
//    print fE1;
    a3 := Coefficient(fE1,3);
    fE1 := Evaluate(fE1,t/a3)*a3^2;
    fE2 := P!(Reverse(coeffsofxsq));
//    print fE2;
    a3 := Coefficient(fE2,3);
    fE2 := Evaluate(fE2,t/a3)*a3^2;
    E1new := EllipticCurve(fE1);
    E2new := EllipticCurve(fE2);
    assert IsIsomorphic(E1,E1new);
/*
    L := RichelotIsogenousSurfaces(genus2curve);
    assert {{x[1],x[2]} : x in L | Type(x) eq SetCart} eq {{E1new,E2new}};
*/
//    return E1, E2new;
    return E1new, E2new, genus2curve;
end intrinsic;

