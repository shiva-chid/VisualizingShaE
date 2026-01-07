/*
For a suppressed cubic f, the Mobius transformation M outputted by KummerElementAndTransformation
seems to have some sort of pattern:
[1 alpha/3]
[alpha 1/n]
where n is an integer, and alpha is some element in K=Q(zeta_3,sqrt(D)) where D is Disc(f).
We find n and alpha below. Turns out n is in Q, and alpha is in Q(sqrt(-3D))
*/

// Finding n
// m := 3;
m := 2;
F<[a]> := FunctionField(Rationals(),m);
P<alpha,n> := PolynomialRing(F,2);
P1<x> := PolynomialRing(P);
coeffsf := (m eq 3) select a cat [F|1] else a cat [F|0,1];
f := P1 ! coeffsf;
g := P1 ! ((n*alpha*x+1)^3*Evaluate(f,(n*x+n*alpha/3)/(n*alpha*x+1)));
I := ideal<P|[Coefficient(g,i) : i in [1,2]]>;
time GBI := GroebnerBasis(I);
[[Degree(x,P.i) : i in [1,2]] : x in GBI];
Factorisation(GBI[#GBI]);
// There is a linear factor: So n is -a1

// Finding alpha
// m := 3;
// P<alpha,a0,a1,a2> := PolynomialRing(Rationals(),m+1);
m := 2;
P<alpha,a0,a1> := PolynomialRing(Rationals(),m+1);
a := [P.(i+1) : i in [1..m]];
P1<x> := PolynomialRing(P);
coeffsf := (m eq 3) select a cat [P|1] else a cat [P|0,1];
f := P1 ! coeffsf;
if m eq 3 then
    n := -(-3*a[1]*a[3] + a[2]^2)/(a[2] - 1/3*a[3]^2);
elif m eq 2 then
    n := -a[2];
end if;
g := P1 ! ((n*alpha*x+1)^3*Evaluate(f,(n*x+n*alpha/3)/(n*alpha*x+1)));
I := ideal<P|[Coefficient(g,i) : i in [1,2]]>;
time GBI := GroebnerBasis(I);
[[Degree(x,P.i) : i in [1..m+1]] : x in GBI];
Factorisation(GBI[#GBI]);
// It is a1 times a quadratic polynomial in alpha with coefficients:
// [-3*a1,-9*a0,a1^2]
// whose discriminant is -3*D where D is disc(f)=-27*a0^2-4*a1^3
// So alpha is a root of this quadratic polynomial.

//////////////////////////////////////

// Verifying in examples that the Mobius transformation outputted by KummerElementAndTransformation is
// exactly of the expected form, with the choice of n and alpha found by the above computation.
ht := 100;
rationalsuptoht := [a/b : a in [-ht..ht], b in [1..ht] | GCD(a,b) eq 1];
rationalsuptoht := Sort(rationalsuptoht,func<x,y|Maximum(Abs(Numerator(x)),Abs(Denominator(x))) - Maximum(Abs(Numerator(y)),Abs(Denominator(y)))>);

// m := 3;
m := 2;
AttachSpec("spec");
P<x> := PolynomialRing(Rationals());
for count := 1 to 30 do
    if m eq 3 then
        f := P![Random(rationalsuptoht),Random(rationalsuptoht),Random(rationalsuptoht),1];
    elif m eq 2 then
        f := P![Random(rationalsuptoht),Random(rationalsuptoht),0,1];
    end if;
    if not IsIrreducible(f) then continue; end if;
    a, M, g := KummerElementAndTransformation(f);
//    print f, a, M, g;
    if m eq 2 then
        a0 := Coefficient(f,0);
        a1 := Coefficient(f,1);
        if a1 ne 0 then
            assert M[2,2] eq -1/a1; // n is -a1
            minpolalpha := MinimalPolynomial(M[2,1]);
            assert Coefficients(minpolalpha) eq [-3/a1,-9*a0/a1^2,1]; // alpha is a root of this quadratic polynomial
        end if;
    end if;
    if count mod 10 eq 0 then print count; end if;
end for;

///////////////////////////////////////////////////

// Writing the transformed cubic (by applying this Mobius transformation)
// explicitly with coefficients in a0, a1, alpha

F<a0,a1> := FunctionField(Rationals(),2);
P<sqrtminus3D> := PolynomialRing(F);
rels := sqrtminus3D^2-(3*(27*a0^2+4*a1^3));
I := ideal<P|rels>;
R, pi := quo<P|I>;
P1<x> := PolynomialRing(R);
alpha := (9*a0+pi(sqrtminus3D))/(2*a1^2);
coeffsf := [R| a0, a1, 0, 1];
f := P1 ! coeffsf;
n := -a1;
g := &+[coeffsf[i]*(n*x+n*alpha/3)^(i-1)*(n*alpha*x+1)^(4-i) : i in [1..#coeffsf]];
g;

