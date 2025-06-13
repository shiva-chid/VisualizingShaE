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
allFs := [];
for abcde in abcdes do
    a,b,c,d,e := Explode(abcde);
    coeffs := [-4*a*c*e + b^2*e + a*d^2, -(4*a*e - b*d), c, 1];
    fE := P!coeffs;
    E := EllipticCurve(fE); Append(~allEs, E);

    F<A,B> := FunctionField(E);
    FF, phi := AlgorithmicFunctionField(F);
    Append(~allFs, FF);
    R<X,Y,Z> := PolynomialRing(FF, 3);
    A := phi(A); B := phi(B);
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
/*
IsMatrixRing(allQAs[3]);
               ^
Runtime error: The base field must be Q, F_q(X) or a number field
*/
fEs := [HyperellipticPolynomials(E) : E in allEs];
fCs := [Evaluate(fEs[1],(t^2+28)/4), Evaluate(fEs[2],(t^2-584)/11), Evaluate(fEs[3],(t^2-1246)/15), Evaluate(fEs[4],(t^2-438)/19)];
Cs := [HyperellipticCurveOfGenus(2,fC) : fC in fCs];
condCs := [];
for C in Cs do
    Append(~condCs, Conductor(C));
    condCs;
end for;

fE := HyperellipticPolynomials(allEs[1]);
fC := Evaluate(fE,(t^2+28)/4);
C := HyperellipticCurveOfGenus(2,fC);
Conductor(C);
out := HeuristicDecomposition(C);
E1 := out[3][2][1][1];
E2 := out[3][2][2][1];
E1 := EllipticCurve(E1);
E2 := EllipticCurve(E2);
jInvariant(E1), jInvariant(E2);
Conductor(E1), Conductor(E2);

//////////////////////////////////////

// trying to check if a given element in the base field of
// a quaternion algebra is a square.

QA := allQAs[1];
PQA<X,Y,Z> := PolynomialRing(QA,3);
genelt := &+[PQA.i*QA.i : i in [1..3]];
gensq := genelt^2; gensq;
assert gensq eq QA.1^2*X^2 + QA.2^2*Y^2 + QA.3^2*Z^2;
FF := allFs[1];
PFF<X,Y,Z> := PolynomialRing(FF,3);
rednrm := (FF!(QA.1^2))*X^2 + (FF!(QA.2^2))*Y^2 + (FF!(QA.3^2))*Z^2;

elt := 19*A+438; // must first choose appropriate isomorphism between elliptic curve models
elt := PFF!elt;
eqntosolve := rednrm-elt;
Factorisation(eqntosolve);
// This is hard to solve, because it is a norm equation over Q(E).
// This question is something like which numbers (functions on E)
// can be written as a specific linear combination of 3 squares (of functions on E)

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

// trying to check if a quaternion algebra over Q(E) is
// split by a specific quadratic extension of Q(E).
// Magma does not seem to support this.
// Can only be done over a global field with char not 2.

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
QA2 := allQAs[2];
PF2<T> := PolynomialRing(F2);
beta := 4*(F2.1-60)-28;
KK := ext<F2|T^2-beta>;
HasEmbedding(KK,QA2);
/*
               ^
Runtime error: The field is not supported.
*/
QA2_KK := BaseExtend(QA2,KK);
IsSplit(QA2_KK);

// boo, isoFF1toFF2 := IsIsomorphic(FF1,FF2);
QA2_FF2 := BaseExtend(QA2,phi2);
QA2_KK := BaseExtend(QA2_FF2,KK);
IsSplit(QA2_KK);

PQA1<tt> := PolynomialRing(QA1);
pol := tt^2 - beta;
Factorisation(pol);
// not supported

/////////////////////////////////////////////////////

E := allEs[1];
jE := jInvariant(E);
assert &and[IsIsomorphic(E,x) : x in allEs[2..#allEs]];
Etors := TorsionSubgroup(E);
if GCD(#Etors,2) eq 1 then printf "E has no rational 2-torsion.\n"; end if;
Emwrk := MordellWeilRank(E);
printf "Mordell-Weil rank of E: %o\n", Emwrk;
K2 := SplittingField(HyperellipticPolynomials(E));
fEs := [HyperellipticPolynomials(E) : E in allEs];
fCs := [Evaluate(fEs[1],(t^2+28)/4), Evaluate(fEs[2],(t^2-584)/11), Evaluate(fEs[3],(t^2-1246)/15), Evaluate(fEs[4],(t^2-438)/19)];
Cs := [HyperellipticCurveOfGenus(2,fC) : fC in fCs];
condCs := [];
for C in Cs do
    Append(~condCs, Conductor(C));
    condCs;
end for;
condCs;
// [ 326041, 15649968, 5216656, 26083280 ]
[cond/571 : cond in condCs];
// [ 571, 27408, 9136, 45680 ]
[Factorisation(x) : x in condCs];

/*
g := (4*t-28)*fEs[1]-1;
Factorisation(g);
*/

for C in Cs do
    time out := HeuristicDecomposition(C);
    E1 := out[3][2][1][1];
    E2 := out[3][2][2][1];
    E1 := EllipticCurve(E1);
    E2 := EllipticCurve(E2);
    Cdecs := [E1,E2];
/*
    jinvs_Cdecs := [jInvariant(E) : E in Cdecs];
    cond_Cdecs := [Conductor(E) : E in Cdecs];
    assert jE in jinvs_Cdecs;
*/
    assert exists(i){i : i in [1..2] | IsIsomorphic(E,Cdecs[i])};
    newind := (i eq 1) select 2 else 1;
    Enew := Cdecs[newind];
    printf "New elliptic curve:\n%o\n", Enew;
    jnew := jInvariant(Enew);
    printf "jInvariant: %o\n", jnew;
    condnew := Conductor(Enew);
    printf "Conductor: %o\n", condnew;
    facscondnew := Factorisation(condnew);
    printf "Conductor factorization: %o\n", facscondnew;
    time lowbd, uppbd := RankBounds(Enew);
    printf "MWRank bounds: %o, %o\n", lowbd, uppbd;
    if {lowbd, uppbd} eq {1} then
        printf "\n\nPartial visualization of H1!\n\n";
    end if;
    // mwgrpnew := MordellWeilGroup(Enew);
    K2new := SplittingField(HyperellipticPolynomials(Enew));
    if IsIsomorphic(K2,K2new) then printf "2-torsion reps are isomorphic.\n"; end if;
    printf "\n==========================================\n\n";
end for;
