E := EllipticCurve("571a1");
F<A,B> := FunctionField(E);
R<X,Y,Z> := PolynomialRing(F, 3);
a := -4;
b := -60;
c := -232;
d := -52;
e := -3;
xi := A;
Q := a*X^2 + b*X*Y + c*Y^2 + d*Y*Z + e*Z^2 + xi*(Y^2 - X*Z);
V := QuadraticSpace(Q);
M := QuadraticFormMatrix(V);
C := CliffordAlgebra(M);
Cplus, emb := EvenSubalgebra(C);
ans, QA, isoCPlusToQA := IsQuaternionAlgebra(Cplus);

P<t> := PolynomialRing(Rationals());
coeffs := [-4*a*c*e + b^2*e + a*d^2, -(4*a*e - b*d), c, 1];
fEE := P!coeffs;
EE := EllipticCurve(fEE);
assert IsIsomorphic(E,EE);
basis_prods := BasisProducts(QA);
basis_prods[4][4] eq -Evaluate(fEE,A);

////////////////////////////////////////////////

// Now we work directly with the model of the elliptic curve
// as given in Section 3.1 of Fisher's "On some algebras" paper.
// The example is from Section 15.1 of Fisher's "Hessian" paper.
P<t> := PolynomialRing(Rationals());
abcdes := [
[-4, -60, -232, -52, -3],
[-11, -68, -52, 164, -64],
[-15, -52, 38, 144, -115],
[-19, 112, -142, -68, -7]
];
for abcde in abcdes do
    a,b,c,d,e := Explode(abcde);
    coeffs := [-4*a*c*e + b^2*e + a*d^2, -(4*a*e - b*d), c, 1];
    fE := P!coeffs;
    E := EllipticCurve(fE);

    F<A,B> := FunctionField(E);
    R<X,Y,Z> := PolynomialRing(F, 3);
    xi := A;
    Q := a*X^2 + b*X*Y + c*Y^2 + d*Y*Z + e*Z^2 + xi*(Y^2 - X*Z);
    V := QuadraticSpace(Q);
    M := QuadraticFormMatrix(V);
    C := CliffordAlgebra(M);
    Cplus, emb := EvenSubalgebra(C);
    ans, QA, isoCPlusToQA := IsQuaternionAlgebra(Cplus);
    basis_prods := BasisProducts(QA);
//    QA;
//    print -basis_prods[4][4]/B^2;
    QAgen, isoQAtoQAgen := ChangeBasis(QA,[1,QA.3/B,QA.1,QA.3*QA.1/B]);
    ans2, QA2, isoQAgentoQA2 := IsQuaternionAlgebra(QAgen);
    QA2;
end for;
/*
Quaternion Algebra with base ring Function Field of E, defined by i^2 = -1, j^2 = 4*A - 28
Quaternion Algebra with base ring Function Field of E, defined by i^2 = -11/4, j^2 = 11*A + 584
Quaternion Algebra with base ring Function Field of E, defined by i^2 = -15/4, j^2 = 15*A + 1246
Quaternion Algebra with base ring Function Field of E, defined by i^2 = -19/4, j^2 = 19*A + 438
*/
