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

// Working directly with the model of the elliptic curve
// as given in Section 3.1 of Fisher
a := -4;
b := -60;
c := -232;
d := -52;
e := -3;
P<t> := PolynomialRing(Rationals());
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
QA;
assert basis_prods[4][4] eq -B^2;
QAgen, isoQAtoQAgen := ChangeBasis(QA,[1,QA.3/B,QA.1,QA.3*QA.1/B]);
ans2, QA2, isoQAgentoQA2 := IsQuaternionAlgebra(QAgen);
QA2;
// Quaternion Algebra with base ring Function Field of E, defined by i^2 = -1, j^2 = 4*A - 28
