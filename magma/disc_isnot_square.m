SetLogFile("disc_isnot_square.log");

FF := Rationals();

// ternary cubic with coefficients of x^3, y^3 and z^3 equal to 1
// and discriminant of the cubic polynomial obtained by intersecting with a general line
PF<xi1,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,a2,b1,c1,c2,m> := PolynomialRing(FF,15);
// PF<xi1,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,a,b,a2,b1,c1,c2,m> := PolynomialRing(FF,17);
xi := [xi1, xi2, xi3, xi4, xi5, xi6, xi7, xi8, xi9, xi10];
P<alpha,beta> := PolynomialRing(PF,2);
P3<x,y,z> := PolynomialRing(P,3);
P1<t> := PolynomialRing(P);
f := x^3+y^3+z^3+a2*x^2*y+b1*y^2*x+c1*z^2*x+c2*z^2*y+m*x*y*z;
// f := a*x^3+b*y^3+z^3+a2*x^2*y+b1*y^2*x+c1*z^2*x+c2*z^2*y+m*x*y*z;
g := Evaluate(f,z,alpha*x+beta*y);
g_t := Evaluate(g,[t,1,0]);
disc := Discriminant(g_t);

// conditions for the discriminant to be a square of a polynomial in Qbar[alpha,beta]
mons := &join[MonomialsOfDegree(P,i) : i in [0..3]];
assert #mons eq 10;
general_sqpoly := (&+[xi[i]*mons[i] : i in [1..#mons]])^2;
wantzero := Coefficients(disc - general_sqpoly);
assert #wantzero eq 28;
I := ideal<PF | wantzero>;
X := Scheme(AffineSpace(PF),I);
time dimX := Dimension(X); dimX;
time GBI := GroebnerBasis(I); #GBI;

// Define ideal of relations in a2,b1,c1,c2,m
somepols := [x : x in GBI | { Degree(x,PF.j) : j in [1..10] } eq {0}]; #somepols;
A5<A2,B1,C1,C2,M> := AffineSpace(FF,5);
P5 := CoordinateRing(A5);
h := hom<PF->P5|[0 : i in [1..10]] cat [A2,B1,C1,C2,M]>;
I5 := ideal<P5|[h(pol) : pol in somepols]>;
I5rad := Radical(I5);
time GBI5 := GroebnerBasis(I5rad); #GBI5;
XX := Scheme(A5,I5rad);
dimXX := Dimension(XX); dimXX;
XX;

// Define Genus One Model
Model := GenusOneModel(f);

// This returns the standard degree 12 invariant Delta
Delta_K := h(Discriminant(Model));

// Check against relations
Delta_K in I5;
Delta_K in I5rad;
