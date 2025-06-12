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
allQAs := [];
allEs := [];
for abcde in abcdes do
    a,b,c,d,e := Explode(abcde);
    coeffs := [-4*a*c*e + b^2*e + a*d^2, -(4*a*e - b*d), c, 1];
    fE := P!coeffs;
    E := EllipticCurve(fE); Append(~allEs, E);

    F<A,B> := FunctionField(E);
//    FF<AA,BB> := AlgorithmicFunctionField(F);
    R<X,Y,Z> := PolynomialRing(F, 3);
    xi := A;
    Q := a*X^2 + b*X*Y + c*Y^2 + d*Y*Z + e*Z^2 + xi*(Y^2 - X*Z);
    V := QuadraticSpace(Q);
    M := QuadraticFormMatrix(V);
    C := CliffordAlgebra(M);
    Cplus, emb := EvenSubalgebra(C);
    ans, QA, isoCPlusToQA := IsQuaternionAlgebra(Cplus);
/*
ii := QA.1;
jj := QA.2;
ii^2;
jj^2;
ii @@ isoCPlusToQA;
(ii @@ isoCPlusToQA)@emb;
jj @@ isoCPlusToQA;
*/
    basis_prods := BasisProducts(QA);
//    QA;
//    print -basis_prods[4][4]/B^2;
    QAgen, isoQAtoQAgen := ChangeBasis(QA,[1,QA.3/B,QA.1,QA.3*QA.1/B]);
    ans2, QA2, isoQAgentoQA2 := IsQuaternionAlgebra(QAgen);
    QA2;
    Append(~allQAs, QA2);
end for;
/*
Quaternion Algebra with base ring Function Field of E, defined by i^2 = -1, j^2 = 4*A - 28
Quaternion Algebra with base ring Function Field of E, defined by i^2 = -11/4, j^2 = 11*A + 584
Quaternion Algebra with base ring Function Field of E, defined by i^2 = -15/4, j^2 = 15*A + 1246
Quaternion Algebra with base ring Function Field of E, defined by i^2 = -19/4, j^2 = 19*A + 438
*/

fC := Evaluate(fE,(t^2+28)/4);
C := HyperellipticCurveOfGenus(2,fC);
out := HeuristicDecomposition(C);
E1 := out[3][2][1][1];
E2 := out[3][2][2][1];
E1 := EllipticCurve(E1);
E2 := EllipticCurve(E2);
jInvariant(E1), jInvariant(E2);
Conductor(E1), Conductor(E2);


PQA<W,X,Y,Z> := PolynomialRing(QA,4);
genelt := W+&+[PQA.i*QA.(i-1) : i in [2..4]];
geneltsq_entries := Eltseq(genelt^2);
elt := 19*A+4;
I := ideal<QA | [geneltsq_entries[1] - elt] cat [geneltsq_entries[2..4]]>;
GBI := GroebnerBasis(I);

//////////////////////////////////////

// Eg. A simple example of a Clifford algebra
// whose even subalgebra is the Hamiltonian algebra.

QQ := Rationals();
R<X,Y,Z> := PolynomialRing(QQ, 3);
Q := X^2 + Y^2 + Z^2;
V := QuadraticSpace(Q);
M := QuadraticFormMatrix(V);
C := CliffordAlgebra(M);
Cplus, emb := EvenSubalgebra(C);
Cplus;

//////////////////////////////////////

allFs := [BaseRing(QA) : QA in allQAs];
QA1 := allQAs[1];
F1 := allFs[1];
beta := QA1.2^2;
// beta := F1!beta;
minpolbeta := MinimalPolynomial(beta);
beta := -Coefficient(minpolbeta,0);
beta in F1;
FF1, phi1 := AlgorithmicFunctionField(F1);
PFF1<T> := PolynomialRing(FF1);
KK := ext<FF1|T^2-phi1(beta)>;


F2 := allFs[2];
FF2, phi2 := AlgorithmicFunctionField(F2);
QA2 := allQAs[2];
PFF2<T> := PolynomialRing(FF2);
beta := phi2(F2!(4*(F2.1-60)-28));
KK := ext<FF1|T^2-beta>;


// boo, isoFF1toFF2 := IsIsomorphic(FF1,FF2);
QA2_FF2 := BaseExtend(QA2,phi2);
QA2_KK := BaseExtend(QA2_FF2,KK);
IsSplit(QA2_KK);


PQA1<tt> := PolynomialRing(QA1);
pol := tt^2 - beta;
Factorisation(pol);
