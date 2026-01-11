AttachSpec("spec");

// The following code produces a Cremona-Mazur cover C ---> E over an at most biquadratic field K
// such that Jac(C) visualizes the torsor.
// To assert minimality of the visualization, we compute zeta function of C at primes p.

// Example 1
L := [1, 1, 3, 5, -7, 15, 2, 19, 1, -2];
C := CMcurveForIndex3Torsor(L);
// C := CMcurveForIndex3Torsor(L : direct := false);
assert Genus(C) eq 4;
p := PrimeWitnessForMinimality(C);
print p;

// Example 2: Diagonal torsor
L := [1, 2, 91, 0, 0, 0, 0, 0, 0, -17]; // A torsor for 182.d1, which is 3-isogenous to 182.d2
C := CMcurveForIndex3Torsor(L);
assert GenusC() eq 1;
PP2<x,y,z> := AmbientSpace(C);
K<w> := BaseRing(PP2);
singptsC := SingularPoints(C); #singptsC;
JacC := EllipticCurve(C,singptsC[1]);
E1 := EllipticCurve([1, 0, 0, -193, -1055]); // 182.d2
E1K := BaseChange(E1,K);
IsIsomorphic(JacC,E1K);


// The following code tries to produce a Cremona-Mazur cover C ---> E over Q such that Jac(C) visualizes the torsor.
// Note that CMcurveOverBaseFieldForIndex3Torsor will not always work, because the model of E
// that Fisher says we need to work with is not necessarily over the base field Q.
// To assert minimality of the visualization, we compute zeta function of C at primes p.
// This takes some time, since the genus of C is 13, which is rather large.
/*
CQ := CMcurveOverBaseFieldForIndex3Torsor(L);
assert Genus(CQ) eq 13;
p := PrimeWitnessForMinimality(CQ);
print p;

PP2<x,y,z> := AmbientSpace(CQ);
CQ;
irrs := IrreducibleComponents(CQ);
irrs := [Curve(x) : x in irrs];
*/

