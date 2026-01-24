intrinsic ApplyFractionalLinearTransformation(f :: RngUPolElt, M :: AlgMatElt : maychangedegree := true) -> RngUPolElt
{apply the fractional linear transformation given by the 2x2 matrix M to transform
the polynomial f of degree 2g+1 or 2g+2}
    deg := Degree(f);
    g := (deg - 1) div 2;
    P<x> := Parent(f);
    a, b, c, d := Explode(Eltseq(M));
    if maychangedegree then
        newdeg := (2*g+2);
    else
        newdeg := deg;
    end if;
    f1 := P!((c*x+d)^(newdeg)*Evaluate(f,(a*x+b)/(c*x+d)));
    return f1;
end intrinsic;

intrinsic FractionalLinearTransformation(f :: RngUPolElt, pts :: SeqEnum : monic := true) -> RngUPolElt, AlgMatElt
{given a polynomial f of degree 2g+1 or 2g+2, and a sequence of 1, 2 or 3 points,
returns a polynomial of degree 2g+1 or 2g+2 obtained by performing a fractional
linear transformation of P^1 that sends the given points in P^1 to infinity, 0 and 1
in order. Also returns a matrix representing the transformation.}
    deg := Degree(f);
    g := (deg - 1) div 2;
    P<x> := Parent(f);
    K := BaseRing(P);
    M := IdentityMatrix(K,2);
    if #pts eq 0 then return f, M; end if;
    if Type(pts[1]) eq Pt then
        pts := [(x[2] ne 0) select [x[1]/x[2],1] else [1,0] where x is Eltseq(pt) : pt in pts];
    end if;
    if Type(pts[1]) eq SeqEnum then
        pts := [(x[2] ne 0) select [x[1]/x[2],1] else [1,0] : x in pts];
    end if;

    if pts[1] ne [1,0] then
        alp := pts[1][1];
        f1 := P ! (x^(2*g+2) * Evaluate(f,alp+1/x));
        pts1 := [(pt ne [1,0]) select [1/(pt[1]-alp),1] else [0,1] : pt in pts[2..#pts]];
        M := M*Matrix(K, 2, 2, [alp,1,1,0]);
    else
        f1 := f; pts1 := pts[2..#pts];
    end if;
//    printf "f1 = %o\npts1 = \n%o\nM = \n%o\n", f1, pts1, M;
    pts1 := [pt[1] : pt in pts1]; // now all remaining points in pts1 are affine.
    if #pts1 ge 1 then
        alp := pts1[1];
//        print alp;
        f1 := Evaluate(f1,x+alp);
        pts1 := [pt-alp : pt in pts1[2..#pts1]];
        M := M*Matrix(K, 2, 2, [1,alp,0,1]);
//        printf "f1 = %o\npts1 = \n%o\nM = \n%o\n", f1, pts1, M;
    end if;
    if #pts1 ge 1 then
        alp := pts1[1];
        f1 := Evaluate(f1,alp*x);
        pts1 := [pt/alp : pt in pts1[2..#pts1]];
        M := M*Matrix(K, 2, 2, [alp,0,0,1]);
//        printf "f1 = %o\npts1 = \n%o\nM = \n%o\n", f1, pts1, M;
    end if;
    if monic then f1 := f1/LeadingCoefficient(f1); end if;
    return f1, M;
end intrinsic;


intrinsic InverseFractionalLinearTransformation(f :: RngUPolElt, pts :: SeqEnum : monic := true) -> RngUPolElt, AlgMatElt
{given a polynomial f of degree 2g+1 or 2g+2, and a sequence of 1, 2 or 3 points,
returns a polynomial of degree 2g+1 or 2g+2 obtained by performing a fractional
linear transformation of P^1 that sends infinity, 0 and 1 to the given points in P^1
in order. Also returns a matrix representing the transformation.}
    deg := Degree(f);
    g := (deg - 1) div 2;
    P<x> := Parent(f);
    K := BaseRing(P);
    M := IdentityMatrix(K,2);
    if #pts eq 0 then return f, M; end if;
    if Type(pts[1]) eq Pt then
        pts := [(x[2] ne 0) select [x[1]/x[2],1] else [1,0] where x is Eltseq(pt) : pt in pts];
    end if;
    if Type(pts[1]) eq SeqEnum then
        pts := [(x[2] ne 0) select [x[1]/x[2],1] else [1,0] : x in pts];
    end if;

    if pts[1] ne [1,0] then
        alp := pts[1][1];
        M := Matrix(K,2,2,[0,1,1,-alp])*M;
        pts1 := [(pt ne [1,0]) select [1/(pt[1]-alp),1] else [0,1]: pt in pts[2..#pts]];
    else
        pts1 := pts[2..#pts];
    end if;
    pts1 := [pt[1] : pt in pts1]; // now all remaining points in pts1 are affine.
//    printf "pts1 = \n%o\nM = \n%o\n", pts1, M;
    if #pts1 ge 1 then
        alp := pts1[1];
        M := Matrix(K,2,2,[1,-alp,0,1])*M;
        pts1 := [pt-alp : pt in pts1[2..#pts1]];
//        printf "pts1 = \n%o\nM = \n%o\n", pts1, M;
    end if;
    if #pts1 ge 1 then
        alp := pts1[1];
        M := Matrix(K,2,2,[1,0,0,alp])*M;
        pts1 := [pt/alp : pt in pts1[2..#pts1]];
//        printf "pts1 = \n%o\nM = \n%o\n", pts1, M;
    end if;
    f1 := ApplyFractionalLinearTransformation(f,M);
    if monic then f1 := f1/LeadingCoefficient(f1); end if;
    return f1, M;
end intrinsic;


/////////////////////////////////////////////////////////////

/*

Example

AttachSpec("magma/spec");
F<a> := FunctionField(Rationals());
P<x> := PolynomialRing(F);
f1 := (x+2)*(x^2-2)*(x^4-4*x^2+a);
f1new, M := FractionalLinearTransformation(f1,[[1,0],[-2,1]]);
Factorisation(f1new);
g1 := x*(x-1)*(x^2-2)*(x^4+a*x^2-a);
g1new, M := FractionalLinearTransformation(g1,[[0,1],[1,1]]);
Factorisation(g1new);

g1old, M2 := InverseFractionalLinearTransformation(g1new,[[0,1],[1,1]]);
g1old; g1; M*M2; M2*M;

*/