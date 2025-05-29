

     //////////////////////////////////////////////////////////
     //                                                      //
     //    On some algebras associated to genus one curves   //
     //                                                      //
     //                   Tom Fisher                         //
     //                                                      //
     //////////////////////////////////////////////////////////

//  This is a MAGMA file, checking the formulae in the paper.
//  Version : June 2017

brac := func<x,y|x*y-y*x>;
Det := Determinant;

function MyDeriv(f,dd)
// We compute D(f) where D is the derivation that sends the
// ith generator to dd[i].
  FF := Parent(f);
  ans := 0;
  for xx in Terms(f) do
    c := Coefficients(xx)[1];
    x := Eltseq(Monomials(xx)[1]);
    word := [FF.x[j]: j in [1..#x]];
    for i in [1..#x] do
      word1 := word;
      word1[i] := dd[x[i]];
      ans +:= c*&*word1;
    end for;
  end for;
  return ans;
end function;

function AlgebraRelations(model)
// INPUT  : A genus one model of degree 2, 3 or 4
// OUTPUT : Relations defining the algebra A_f, together with
//          the central elements xi and eta.
  assert Degree(model) in [2,3,4];
  K := BaseRing(model);
  case Degree(model) :
  when 2 :
    FF<r,s,t> := FreeAlgebra(K,3);
    assert #Eltseq(model) eq 5;
    a,b,c,d,e := Explode(Eltseq(model));
    relns := [
      r^2 - a,
      r*s + s*r - b,
      r*t + t*r + s^2 - c,
      s*t + t*s - d,
      t^2 - e ];
    xi := s^2 - c;
    eta := r*s*t - t*s*r;
  when 3 :
    FF<x,y> := FreeAlgebra(K,2);
    a,b,c,a2,a3,b1,b3,c1,c2,m := Explode(Eltseq(model));
    relns := [
      c*x^3 + c1*x^2 + a3*x + a,
      c*(x^2*y + x*y*x + y*x^2) + c1*(x*y + y*x)
                           + c2*x^2 + m*x + a3*y + a2,
      c*(x*y^2 + y*x*y + y^2*x)+ c2*(x*y + y*x)
                           + c1*y^2 + m*y + b3*x +b1,
      c*y^3 + c2*y^2 + b3*y + b ];
    r := y*(c*x^2 + c1*x + a3);
    s := -(c*y^2 + c2*y + b3);
    t := c*x;
    xi := (c*x*y)^2 - (c*y^2 + c2*y + b3)*(c*x^2 + c1*x + a3)
                                    + (c*m - c1*c2)*x*y + a3*b3;
    eta := (r*s*t + s*t*r + t*r*s)
      + a2*(s*t + t*s) + b3*(t*r + r*t) + c1*(r*s + s*r)
      + (b3*c1 - b1*c)*r + (c1*a2 - c2*a)*s + (a2*b3 - a3*b)*t
      - 6*a*b*c + a2*b3*c1;
  when 4 :
    FF<p,q,r,s> := FreeAlgebra(K,4);
    a11, a12, a13, a14, a22, a23, a24, a33, a34, a44,
    b11, b12, b13, b14, b22, b23, b24, b33, b34, b44
      := Explode(Eltseq(model));
    relns := [
      a11*p^2 + a12*p*q + a22*q^2 + a13*p + a23*q + a33,
      a11*(p*r + r*p) + a12*(p*s + r*q) + a22*(q*s + s*q)
          + a14*p + a24*q + a13*r + a23*s + a34,
      a11*r^2 + a12*r*s + a22*s^2 + a14*r + a24*s + a44,
      b11*p^2 + b12*p*q + b22*q^2 + b13*p + b23*q + b33,
      b11*(p*r + r*p) + b12*(p*s + r*q) + b22*(q*s + s*q)
          + b14*p + b24*q + b13*r + b23*s + b34,
      b11*r^2 + b12*r*s + b22*s^2 + b14*r + b24*s + b44,
      p*q - q*p,
      p*s + r*q - q*r - s*p,
      r*s - s*r ];
    mat := Matrix(2,10,Eltseq(model));
    D := Matrix(10,10,[mat[1,i]*mat[2,j]-mat[1,j]*mat[2,i]: i,j in [1..10]]);
    d12,d13,d14,d15,d16,d17,d18,d19,d110,
        d23,d24,d25,d26,d27,d28,d29,d210,
            d34,d35,d36,d37,d38,d39,d310,
                d45,d46,d47,d48,d49,d410,
                    d56,d57,d58,d59,d510,
                        d67,d68,d69,d610,
                            d78,d79,d710,
                                d89,d810,
                                    d910
       := Explode([D[i,j] : i,j in [1..10] | i lt j]);
    pp := [D[1,i]*p + D[2,i]*q + D[3,i] : i in [1..10]];
    qq := [D[2,i]*p + D[5,i]*q + D[6,i] : i in [1..10]];
    rr := [D[1,i]*r + D[2,i]*s + D[4,i] : i in [1..10]];
    ss := [D[2,i]*r + D[5,i]*s + D[7,i] : i in [1..10]];
    p1,p2,p3,p4,p5,p6,p7,p8,p9,p10 := Explode(pp);
    q1,q2,q3,q4,q5,q6,q7,q8,q9,q10 := Explode(qq);
    r1,r2,r3,r4,r5,r6,r7,r8,r9,r10 := Explode(rr);
    s1,s2,s3,s4,s5,s6,s7,s8,s9,s10 := Explode(ss);
    t := q*r - p*s;
    xi := (p5*s)^2 + (s1*p)^2
      + (d56*p4 + d29*p5 + d37*p5 - d27*p6)*s
      + (d14*s6 + d29*s1 - d37*s1 - d23*s4)*p
      - d56*(d13*r + d23*s - d17*q + d12*t - d19)*s
      - d14*(d27*p + d57*q + d35*r - d25*t - d59)*p;
    eps := d12*(p*r + r*p) + d15*(p*s + q*r + s*p + r*q) + d25*(q*s + s*q);
    dd := Vector([(1/2)*brac(p,eps),(1/2)*brac(q,eps),0,0]);
    eta := (1/2)*MyDeriv(xi,dd);
    c0 := (1/2)*(a12*a34 + a13*a24 + a14*a23)*(b12*b34 + b13*b24 + b14*b23)
      - (a11*b22 - a12*b12 + a22*b11)*(a33*b44 - a34*b34 + a44*b33)
      - (a11*b33 - a13*b13 + a33*b11)*(a44*b22 - a24*b24 + a22*b44)
      - (a11*b44 - a14*b14 + a44*b11)*(a22*b33 - a23*b23 + a33*b22)
      - (a11*a34*b23*b24 + a22*a14*b13*b34
             + a33*a12*b14*b24 + a44*a12*b13*b23)
      - (b11*b34*a23*a24 + b22*b14*a13*a34
             + b33*b12*a14*a24 + b44*b12*a13*a23)
      + 2*(a11*a33*b22*b44 + a22*a44*b11*b33)
      - (2*d15*d810 - d16*d79 + d17*d69 + d17*d78 - d23*d610
         + d24*d78 + d26*d310 - d27*d48 - d28*d210 - d35*d310)
      + (-d16*d79 + d19*d59 - d23*d610 - d27*d48);
  end case;
  if Degree(model) in [2,3] then
    return relns,xi,eta,_;
  else
    return relns,xi,eta,c0;
  end if;
end function;

function GetAlgebra(model)
  relns,xi,eta,c0 := AlgebraRelations(model);
  F := Universe(relns);
  A := quo<F|ideal<F|relns>>;
  if Degree(model) in [2,3] then
    return A,A!xi,A!eta,_;
  else
    return A,A!xi,A!eta,c0;
  end if;
end function;

IsCentral := func<x|forall{i : i in [1..Rank(A)] |
             A.i*x eq x*A.i} where A is Parent(x)>;

// We check in some random numerical examples that xi and eta are
// indeed central elements of Af and that they satisfy a Weierstrass
// equation for the Jacobian of Cf.

for n in [2..4] do
  f := RandomModel(n: RequireNonsingular);
  Af,xi,eta,c0 := GetAlgebra(f);
  assert IsCentral(xi) and IsCentral(eta);
  if n in [2,3] then
    a1,a2,a3,a4,a6 := Explode(aInvariants(f));
    assert eta^2 + a1*xi*eta + a3*eta eq xi^3 + a2*xi^2 + a4*xi + a6;
  else
    MM := Matrices(f);
    P<X> := PolynomialRing(Rationals());
    g := (1/4)*Det(MM[1] + X*MM[2]);
    aa,bb,cc,dd,ee := Explode([Coefficient(g,i): i in [0..4]]);
    F := X^3 + cc*X^2 - (4*aa*ee-bb*dd)*X - 4*aa*cc*ee + bb^2*ee + aa*dd^2;
    assert eta^2 eq Evaluate(F,xi + c0);
  end if;
end for;

////////////////////////////////////////////////////////////////

// The headings (3.1), (3.2), etc, refer to the subsections in the paper.

// (3.1) Binary quartics

K<a,b,c,d,e> := FunctionField(Rationals(),5);
f := GenusOneModel(K,2,[a,b,c,d,e]);
Af,xi,eta := GetAlgebra(f);
assert IsCentral(xi) and IsCentral(eta);
assert eta^2 eq xi^3 + c*xi^2
                   - (4*a*e-b*d)*xi - 4*a*c*e + b^2*e + a*d^2;
F<r,s,t> := PreimageRing(Af);
dd := [brac(s,r),brac(t,r),0];
relns := Basis(DivisorIdeal(Af));
assert forall{rr : rr in relns| Af!MyDeriv(rr,dd) eq 0};
assert MyDeriv(xi,dd) eq 2*eta;
assert MyDeriv(eta,dd) eq 3*xi^2 + 2*c*xi - 4*a*e + b*d;

// (3.2) Ternary cubics

K<a,b,c,a2,a3,b1,b3,c1,c2,m> := FunctionField(Rationals(),10);
f := GenusOneModel(K,3,[a,b,c,a2,a3,b1,b3,c1,c2,m]);
Af,xi,eta := GetAlgebra(f);
assert IsCentral(xi) and IsCentral(eta);
a1,a2,a3,a4,a6 := Explode(aInvariants(f));
assert eta^2 + a1*xi*eta + a3*eta eq xi^3 + a2*xi^2 + a4*xi + a6;
F<x,y> := PreimageRing(Af);
dd := [c*brac(x*y,x),c*brac(y,y*x)];
relns := Basis(DivisorIdeal(Af));
assert forall{rr : rr in relns| Af!MyDeriv(rr,dd) eq 0};
assert MyDeriv(xi,dd) eq 2*eta + a1*xi + a3;
assert MyDeriv(eta,dd) eq 3*xi^2 + 2*a2*xi + a4 - a1*eta;

// (3.3) Quadric intersections

K< a11, a12, a13, a14, a22, a23, a24, a33, a34, a44,
   b11, b12, b13, b14, b22, b23, b24, b33, b34, b44 >
         := FunctionField(Rationals(),20);
f := GenusOneModel(K,4,[
   a11, a12, a13, a14, a22, a23, a24, a33, a34, a44,
   b11, b12, b13, b14, b22, b23, b24, b33, b34, b44]);
relns,xi,eta,c0 := AlgebraRelations(f);
F<p,q,r,s> := Universe(relns);
I := ideal<F|relns>;
print "Case n = 4: Checking xi is central"; // takes a few seconds
time assert forall{ i : i in [1..4] | F.i*xi - xi*F.i in I};
mat := Matrix(2,10,Eltseq(f));
D := Matrix(10,10,[mat[1,i]*mat[2,j]-mat[1,j]*mat[2,i]: i,j in [1..10]]);
d12,d15,d25 := Explode([D[1,2],D[1,5],D[2,5]]);
eps := d12*(p*r + r*p) + d15*(p*s + q*r + s*p + r*q) + d25*(q*s + s*q);
dd := [(1/2)*brac(p,eps),(1/2)*brac(q,eps),0,0];
print "Case n = 4: Checking D is well defined";
time assert forall{rr : rr in relns| MyDeriv(rr,dd) in I};

// (4.1) Binary quartics

K<a,b,c,d,e,m11,m12,m21,m22,la> := FunctionField(Rationals(),10);
f := GenusOneModel(2,[a,b,c,d,e]);
_,gamma := IsTransformation(2,<la,Matrix(2,2,[m11,m12,m21,m22])>);
f1 := gamma*f;
R2<x,z> := PolynomialRing(f);
subst := [m11*x + m21*z,m12*x + m22*z];
assert Equation(f1) eq la^2*Evaluate(Equation(f),subst);
relns,xi := AlgebraRelations(f);
relns1,xi1 := AlgebraRelations(f1);
FF<r,s,t> := Universe(relns);
I := ideal<FF|relns>;
psi := func<x|Evaluate(x,[
  la*(m11^2*r + m11*m12*s + m12^2*t),
  la*(2*m11*m21*r + (m11*m22 + m12*m21)*s + 2*m12*m22*t),
  la*(m21^2*r + m21*m22*s + m22^2*t)])>;
assert forall{ eqn : eqn in relns1 | psi(eqn) in I};
rho := -la^2*(2*m11^2*m21^2*a + m11*m21*(m11*m22 + m12*m21)*b
   + 2*m11*m12*m21*m22*c + m12*m22*(m11*m22 + m12*m21)*d + 2*m12^2*m22^2*e);
assert psi(xi1) - (Det(gamma)^2*xi + rho) in I;
kappa := la*(m11*m21*r + m12*m21*s + m12*m22*t);
D := func<x|MyDeriv(x,[brac(s,r),brac(t,r),0])>;
for x in [r,s,t] do
  assert psi(D(x)) eq Det(gamma)*D(psi(x)) + brac(kappa,psi(x));
end for;

// (4.2) Ternary cubics

K<a,b,c,a2,a3,b1,b3,c1,c2,m,
   m11,m12,m13,m21,m22,m23,m33> := FunctionField(Rationals(),17);
f := GenusOneModel(3,[a,b,c,a2,a3,b1,b3,c1,c2,m]);
R3<x,y,z> := PolynomialRing(f);
M := [Matrix(3,3,[m11,m12,m13,m21,m22,m23,0,0,m33]),
      Matrix(3,3,[0,1,0,0,0,1,1,0,0])];
cobs := [[m11*x + m21*y,m12*x + m22*y,m13*x + m23*y + m33*z],[z,x,y]];
for ctr in [1,2] do
  _,gamma := IsTransformation(3,<1,M[ctr]>);
  f1 := gamma*f;
  assert Equation(f1) eq Evaluate(Equation(f),cobs[ctr]);
  relns,xi := AlgebraRelations(f);
  relns1,xi1 := AlgebraRelations(f1);
  FF<x,y> := Universe(relns);
  I := ideal<FF|relns>;
  if ctr eq 1 then
    psi := func<a|Evaluate(a,
               [(1/m33)*(m11*x + m12*y - m13),
                (1/m33)*(m21*x + m22*y - m23)])>;
    kappa := c*m33*(m23*(m11*x + m12*y) - m13*(m21*x + m22*y));
  else
    xinv := -(1/a)*(c*x^2 + c1*x + a3);
    psi := func<a|Evaluate(a,[-y*xinv,xinv])>;
    kappa := c*y*x + c1*y;
  end if;
  assert forall{ eqn : eqn in relns1 | psi(eqn) in I};
  rho := K!NormalForm(psi(xi1) - Det(gamma)^2*xi,I);
  assert Denominator(rho) eq 1;
  cc := Eltseq(f1)[3];
  D := func<a|MyDeriv(a,[c*brac(x*y,x),c*brac(y,y*x)])>;
  D1 := func<a|MyDeriv(a,[cc*brac(x*y,x),cc*brac(y,y*x)])>;
  for a in [x,y] do
    assert psi(D1(a)) - Det(gamma)*D(psi(a)) - brac(kappa,psi(a)) in I;
  end for;
end for;

// (4.3) Quadric intersections

K < a11, a12, a13, a14, a22, a23, a24, a33, a34, a44,
    b11, b12, b13, b14, b22, b23, b24, b33, b34, b44,
    A11,A12,A21,A22 > := FunctionField(Rationals(),24);
f := GenusOneModel(4,[
    a11, a12, a13, a14, a22, a23, a24, a33, a34, a44,
    b11, b12, b13, b14, b22, b23, b24, b33, b34, b44]);
R4<x1,x2,x3,x4> := PolynomialRing(f);
A := Matrix(2,2,[A11,A12,A21,A22]);
I2 := IdentityMatrix(K,2);
MM := [BlockMatrix(2,2,[A^(-1),0,0,I2]),
       BlockMatrix(2,2,[I2,0,A,I2]),
       Matrix(K,4,4,[<1,3,1>,<2,4,1>,<3,1,1>,<4,2,1>])];

for ctr in [1..3] do
  M := MM[ctr];

  if ctr eq 3 then
    K := Rationals(); // Switching to a numerical example
    f := RandomModel(4:CoeffSet := [-30..30],RequireNonsingular);
    f := ChangeRing(f,K);
    M := ChangeRing(M,K);
  end if;

  _,gamma := IsTransformation(4,<IdentityMatrix(K,2),M>);
  f1 := gamma*f;
  relns,xi := AlgebraRelations(f);
  relns1,xi1 := AlgebraRelations(f1);
  FF<p,q,r,s> := Universe(relns);
  I := ideal<FF|relns>;

  mat := Matrix(2,10,Eltseq(f));
  D := Matrix(10,10,[mat[1,i]*mat[2,j]-mat[1,j]*mat[2,i]: i,j in [1..10]]);
  eps := D[1,2]*(p*r + r*p) + D[1,5]*(p*s + q*r + s*p + r*q)
                    + D[2,5]*(q*s + s*q);
  mat := Matrix(2,10,Eltseq(f1));
  D1 := Matrix(10,10,[mat[1,i]*mat[2,j]-mat[1,j]*mat[2,i]: i,j in [1..10]]);
  eps1 := D1[1,2]*(p*r + r*p) + D1[1,5]*(p*s + q*r + s*p + r*q)
                    + D1[2,5]*(q*s + s*q);
  if ctr eq 1 then
    psi := func<x|Evaluate(x,
        [A11*p + A21*q,A12*p + A22*q,A11*r + A21*s,A12*r + A22*s])>;
    assert psi(eps1) eq Det(gamma)*eps;
    kappa := 0;
  elif ctr eq 2 then
    psi := func<x|Evaluate(x,[p - A11,q - A12,r - A21, s - A22])>;
    kappa := A11*(D[1,2]*r + D[1,5]*s) + A12*(D[1,5]*r + D[2,5]*s);
  elif ctr eq 3 then
    pp := [D[1,i]*p + D[2,i]*q + D[3,i] : i in [1..10]];
    qq := [D[2,i]*p + D[5,i]*q + D[6,i] : i in [1..10]];
    rr := [D[1,i]*r + D[2,i]*s + D[4,i] : i in [1..10]];
    ss := [D[2,i]*r + D[5,i]*s + D[7,i] : i in [1..10]];
    t := q*r - p*s;
    ttt := - D[8,9]*(ss[1]*r + ss[4])
           - D[8,10]*(rr[5]*q + rr[6] + pp[5]*s + pp[7] + D[2,9])
           - D[9,10]*(qq[1]*p + qq[3]);
    tinv := 1/(D[8,10]^2 - D[8,9]*D[9,10])*ttt;
    subst := [-s*tinv,q*tinv,r*tinv,-p*tinv];
    subst := [NormalForm(v,I): v in subst];
    psi := func<x|Evaluate(x,subst)>;
    Delta := D[8,10]^2 - D[8,9]*D[9,10];
    la := (2*D[4,8]*D[8,10] - D[3,8]*D[9,10] - D[8,9]*D[4,9]
                 + D[8,9]*D[3,10])/(2*Delta);
    mu := (2*D[3,10]*D[8,10] - D[4,10]*D[8,9] - D[9,10]*D[3,9]
                 + D[9,10]*D[4,8])/(2*Delta);
    kappa := la*(p*(ss[1]*r + ss[4]) + qq[10])
           + mu*(r*(qq[1]*p + qq[3]) + ss[8])
       + r*(D[1,2]*p + D[1,5]*q) - (1/2)*(D[2,3]*r + D[2,6]*s);
  end if;

  if ctr lt 3 then
    assert forall{ eqn : eqn in relns1 | psi(eqn) in I};
    rho := K!NormalForm(psi(xi1) - Det(gamma)^2*xi,I);
    assert Denominator(rho) in [1,Det(A)^4];
  else
    xit := (D[1,5]^2 - D[1,2]*D[2,5])*t^2 + (D[1,5]*D[3,7]
      - D[1,2]*D[6,7] - D[1,5]*D[4,6] - D[2,5]*D[3,4])*t
      + (D[3,7]*D[8,10] - D[3,6]*D[9,10] + D[4,6]*D[8,10]
      - D[4,7]*D[8,9])*tinv + (D[8,10]^2 - D[8,9]*D[9,10])*tinv^2;
    c1 := K!NormalForm(xi - xit,I);
    assert Denominator(c1) eq 1;
    assert psi(t) - tinv in I;
  end if;

  D := func<x|MyDeriv(x,[(1/2)*brac(p,eps),(1/2)*brac(q,eps),0,0])>;
  D1 := func<x|MyDeriv(x,[(1/2)*brac(p,eps1),(1/2)*brac(q,eps1),0,0])>;
  for x in [p,q,r,s] do
    assert psi(D1(x)) - Det(gamma)*D(psi(x)) - brac(kappa,psi(x)) in I;
  end for;

end for;

// (5.1) Binary quartics

K<a2,a4,a6> := FunctionField(Rationals(),3);
R2<x,y> := PolynomialRing(K,2);
F := y^2 - (x^3 + a2*x^2 + a4*x + a6);
f := GenusOneModel(a6*x^4 + a4*x^3*y + a2*x^2*y^2 + x*y^3);
relns,xi,eta := AlgebraRelations(f);
FF<r,s,t> := Universe(relns);

ZZ<x0,y0> := quo<R2|F>;
E := [[Matrix(ZZ,2,2,[<i,j,1>]): j in [1..2]]: i in [1..2]];
r1 := Matrix(ZZ,2,2,[-y0,x0^2 + a2*x0 + a4,-x0,y0]);
s1 := Matrix(ZZ,2,2,[0,x0 + a2,1,0]);
t1 := Matrix(ZZ,2,2,[0,1,0,0]);
assert forall{reln : reln in relns | Evaluate(reln,[r1,s1,t1]) eq 0};
assert Evaluate(xi,[r1,s1,t1]) eq x0;
assert Evaluate(eta,[r1,s1,t1]) eq y0;
assert Evaluate(t,[r1,s1,t1]) eq E[1,2];
assert Evaluate(s - s^2*t,[r1,s1,t1]) eq E[2,1];

A<r,s,t> := quo<FF|relns>;
xi := A!xi;
eta := A!eta;
E12 := t;
E21 := s - s^2*t;
assert IsCentral(xi);
assert IsCentral(eta);
assert Evaluate(F,[xi,eta]) eq 0;
assert E12^2 eq 0;
assert E21^2 eq 0;
assert E12*E21 + E21*E12 eq 1;
assert -eta*E12*E21 + (xi^2 + a2*xi + a4)*E12 - xi*E21 + eta*E21*E12 eq r;
assert (xi + a2)*E12 + E21 eq s;
assert E12 eq t;

// (5.2) Ternary cubics

K<a1,a2,a3,a4,a6> := FunctionField(Rationals(),5);
R2<x,y> := PolynomialRing(K,2);
F := y^2 + a1*x*y + a3*y - (x^3 + a2*x^2 + a4*x + a6);
R3<x,y,z> := PolynomialRing(K,3);
cubic := R3!(x^3*Evaluate(F,[z/x,y/x]));
f := GenusOneModel(cubic);
relns,xi,eta := AlgebraRelations(f);
FF<x,y> := Universe(relns);

ZZ<x0,y0> := quo<R2|F>;
E := [[Matrix(ZZ,3,3,[<i,j,1>]): j in [1..3]]: i in [1..3]];
x1 := Matrix(ZZ,3,3,[-x0-a2,-1,0,x0^2+a2*x0+a4,0,y0,y0+a1*x0+a3,0,x0]);
y1 := Matrix(ZZ,3,3,[0,0,0,-a1,0,-1,1,0,0]);
assert forall{reln : reln in relns | Evaluate(reln,[x1,y1]) eq 0};
assert Evaluate(xi,[x1,y1]) eq x0;
assert Evaluate(eta,[x1,y1]) eq y0;
assert Evaluate(-x*y^2*(x + xi + a2),[x1,y1]) eq E[1,2];
assert Evaluate(-y^2*(x*y - a1),[x1,y1]) eq E[2,3];
assert Evaluate((y*x - a1)*y^2,[x1,y1]) eq E[3,1];

A<x,y> := quo<FF|relns>;
xi := A!xi;
eta := A!eta;
E12 := -x*y^2*(x + xi + a2);
E23 := -y^2*(x*y - a1);
E31 := (y*x - a1)*y^2;
assert IsCentral(xi);
assert IsCentral(eta);
assert Evaluate(F,[xi,eta]) eq 0;
assert [E12^2,E23^2,E31^2] eq [0,0,0];
assert [E12*E31,E23*E12,E31*E23] eq [0,0,0];
assert E12*E23*E31 + E23*E31*E12 + E31*E12*E23 eq 1;
assert -(xi + a2)*E12*E23*E31 - E12 + (xi^2 + a2*xi + a4)*E23*E31
       + eta*E23 + (eta + a1*xi + a3)*E31 + xi*E31*E12*E23 eq x;
assert -a1*E23*E31 - E23 + E31 eq y;

// (5.3) Quadric intersections

K<a,b,c,d,e> := FunctionField(Rationals(),5);
R2<x,z> := PolynomialRing(K,2);
quartic := a*x^4 + b*x^3*z + c*x^2*z^2 + d*x*z^3 + e*z^4;
f1 := GenusOneModel(quartic);
relns1,xi1,eta1 := AlgebraRelations(f1);
FF<r1,s1,t1> := Universe(relns1);
R4<x1,x2,x3,x4> := PolynomialRing(K,4);
f := GenusOneModel([x1^2 - x3*x4,
        x2^2 - (a*x3^2 + b*x1*x3 + c*x3*x4 + d*x1*x4 + e*x4^2)]);
relns,xi,eta,c0 := AlgebraRelations(f);
FF<p,q,r,s> := Universe(relns);

Af1 := quo<Universe(relns1)|relns1>;
pp := Matrix(Af1,2,2,[0,1,0,0]);
qq := Matrix(Af1,2,2,[r1,s1,0,r1]);
rr := Matrix(Af1,2,2,[0,0,1,0]);
ss := Matrix(Af1,2,2,[t1,0,s1,t1]);
assert forall{reln : reln in relns | Evaluate(reln,[pp,qq,rr,ss]) eq 0};
assert Evaluate(xi,[pp,qq,rr,ss]) eq Af1!(xi1 + c);
assert Evaluate(eta,[pp,qq,rr,ss]) eq -Af1!eta1;
assert pp*qq*rr + rr*qq*pp eq Af1!r1;
assert pp*rr*qq*rr + rr*qq*rr*pp eq Af1!s1;
assert pp*ss*rr + rr*ss*pp eq Af1!t1;

A<p,q,r,s> := quo<FF|relns>;
E12 := p;
E21 := r;
r1 := p*q*r + r*q*p;
s1 := p*r*q*r + r*q*r*p;
t1 := p*s*r + r*s*p;
assert forall{reln : reln in relns1 | Evaluate(reln,[r1,s1,t1]) eq 0};
assert forall{[x,y] : x in [r1,s1,t1],y in [E12,E21] | x*y eq y*x};
assert E12^2 eq 0;
assert E21^2 eq 0;
assert r1 + s1*E12 eq q;
assert t1 + s1*E21 eq s;

AA,BB := Explode(Matrices(f));
f1 := GenusOneModel((1/4)*Det(AA*x + BB*z));
assert Evaluate(Equation(Jacobian(f1)),[A!xi,A!eta,1]) eq 0;

// (6.3) Binary quartics

K<a,c,d,e> := FunctionField(Rationals(),4);
R2<x,z> := PolynomialRing(K,2);
f := GenusOneModel(a*x^4 + c*x^2*z^2 + d*x*z^3 + e*z^4);
P<X> := PolynomialRing(K);
L<aa> := quo<P|X^2 - a>;
E := Jacobian(f);
C,E,pi := nCovering(f:E := E);
Q := C(L)![1,0,aa];
P := E(L)![-c,d*aa];
assert pi(Q) in [P,-P];

// (6.4) Ternary cubics

K0<z3> := CyclotomicField(3);
K<a,b,b1,b3,m> := FunctionField(K0,5);
R2<x,y> := PolynomialRing(K,2);
g := y - z3^2*m*x - 3*(1 - z3)*a*b + (1/9)*(z3 - z3^2)*m^3;
R3<x,y,z> := PolynomialRing(K,3);
f := GenusOneModel(a*x^3 + b*y^3 - z^3 + b1*x*y^2 + b3*y^2*z + m*x*y*z);
relns,xi,eta := AlgebraRelations(f);
A<x,y> := quo<Universe(relns)|relns>;
v := y*x - z3*x*y - (1/3)*(1-z3)*m;
assert x^3 eq a;
assert x*v eq z3*v*x;
assert v^3 eq Evaluate(g,[A!xi,A!eta]);
P<X> := PolynomialRing(K);
L<aa> := quo<P|X^3 - a>;
E := Jacobian(f);
C,E,pi := nCovering(f:E := E);
xR := -(1/3)*m^2 + b1*aa - b3*aa^2;
yR := -z3^2*b3*m*aa^2 + z3^2*b1*m*aa + 3*(1 - z3)*a*b + 1/9*(2 + z3)*m^3;
R := E(L)![xR,yR];
sigma := hom<L->L|z3*aa>;
R1 := E(L)![sigma(xR),sigma(yR)];
R2 := E(L)![sigma(sigma(xR)),sigma(sigma(yR))];
Q := C(L)![1,0,aa];
assert pi(Q) in [R1 - R2,R2 - R1];
assert forall{P : P in [R,R1,R2] | Evaluate(g,[P[1],P[2]]) eq 0};

// (6.6) Quadric intersections

K<a,b,c,d1,d2,d3,d4,e> := FunctionField(Rationals(),8);
delta := a*c*(b^2 - 4*a*c);
P<X> := PolynomialRing(K);
K1<sqrtdelta> := quo<P|X^2 - delta>;
R2<x,y> := PolynomialRing(K1,2);
m := c*d1^2 - b*d1*d3 + a*d3^2 + (b^2 - 4*a*c)*d4;
n := b*c*d1^2 + a*c*(d2^2 - 4*d1*d3) + a*b*d3^2;
g := x^2 - ((4*a*c*d2)/sqrtdelta)*y
       + (2*(b*m + 2*a*c*d2^2)/(b^2 - 4*a*c))*x
       + 8*a*c*e*d2^2 + (m^2 + d2^2*n)/(b^2 - 4*a*c);
P<X> := PolynomialRing(K1);
L<theta> := quo<P|a*X^4 + b*X^2 + c>;
theta1 := (a*b*theta^3 + (b^2 - 2*a*c)*theta)/sqrtdelta;
assert a*theta1^4 + b*theta1^2 + c eq 0;
P<X> := PolynomialRing(K);
FF<phi> := quo<P|(X^2 + a*b)^2 - 4*a^3*c>;
xP := (2*a*phi^2*(d2*d3*phi + m) - d2*(d1*phi + a*d2)*(b*phi^2
                       + a*(b^2 - 4*a*c)))/(2*a^2*(b^2 - 4*a*c));
iota := hom<FF->L|a*(theta + theta1)>;

R4<x1,x2,x3,x4> := PolynomialRing(K,4);
f := GenusOneModel([a*x1^2 + b*x1*x3 + c*x3^2
    + (d1*x1 + d2*x2 + d3*x3 + d4*x4)*x4,x2^2 - x1*x3 - e*x4^2]);
fL := ChangeRing(f,L);
relns,xi,eta,c0 := AlgebraRelations(fL);
A<p,q,r,s> := quo<Universe(relns)|relns>;
q1 := (a*b*q^3 + (b^2 - 2*a*c)*q)/sqrtdelta;
qinv := -(a*q^3 + b*q)/c;
assert q*qinv eq 1;
v := a*(r*q1 - q*r)
     + (a/sqrtdelta)*(a*q^3*q1 - c)*(d1*q + d2 + d3*qinv);
vbar := a*(-r*q1 - q*r)
     - (a/sqrtdelta)*(-a*q^3*q1 - c)*(d1*q + d2 + d3*qinv);
assert a*q^4 + b*q^2 + c eq 0;
assert v*q eq q1*v;
assert v^4 eq Evaluate(g,[A!(xi + c0),A!eta]);

phi0 := a*(q + q1);
dd := a^2*(b^2 - 4*a*c);
xP0 := &+[dd*Eltseq(xP)[i]*phi0^(i-1): i in [1..4]]/dd;
assert v*vbar eq A!(xi + c0) - xP0;

C := Curve(f);
Q := C(L)![theta^2,theta,1,0];
Q1 := C(L)![theta1^2,theta1,1,0];
M := Matrix(4,4,[a,phi,(phi^2 + a*b)/(2*a),0,
   a,-phi,(phi^2 + a*b)/(2*a),0,
   -a*d1,-a*d2,-a*d3,-a*d4 - e*phi^2,
   0,0,0,1]) where phi is iota(phi);
Minv := M^(-1);
S<z1,z2> := PolynomialRing(L,2);
P<T> := PolynomialRing(S);
subst := [&+[Minv[i,j]*[T*z1,z2,z1,T*z2][j]
                          : j in [1..4]]: i in [1..4]];

f1,f2 := Explode(Equations(fL));
f0 := (b^2 - 4*a*c)*f1 + (4*a*sqrtdelta*theta^2
               + 2*b*sqrtdelta - 4*a*b*c + b^3)*f2;
assert Evaluate(f0,subst) eq 0;
qq := Det(M)*Evaluate(f2,subst);
aa,bb,cc := Explode([Coefficient(qq,i): i in [2,1,0]]);
quartic := GenusOneModel(bb^2 - 4*aa*cc);
Q := [&+[M[i,j]*Q[j] : j in [1..4]]: i in [1..4]];
Q1 := [&+[M[i,j]*Q1[j] : j in [1..4]]: i in [1..4]];
assert Q[2] eq 0 and Q[4] eq 0;
assert Q1[2] eq 0 and Q1[4] eq 0;
assert Eltseq(quartic)[1] eq a^2*(theta - theta1)^2;

BQ := GenusOneModel(FF,2,[ -phi^2 - 2*a*b, -2*d1*phi^2 - 2*a*d2*phi -
4*a^2*d3 - 2*a*b*d1, -6*a*d4*phi^2 + 8*a^3*c*e - 2*a^2*b^2*e -
4*a^2*b*d4 - 4*a^2*d1*d3 + a^2*d2^2, -6*a*d2*e*phi^3 + (4*a^2*d3*e -
2*a*b*d1*e - 2*a*d1*d4)*phi^2 + (-8*a^2*b*d2*e + 2*a^2*d2*d4)*phi +
8*a^3*c*d1*e - 4*a^3*d3*d4 - 2*a^2*b^2*d1*e - 2*a^2*b*d1*d4,
-2*a*d1*d2*e*phi^3 + (-4*a^3*c*e^2 + a^2*b^2*e^2 + 2*a^2*d2^2*e -
a^2*d4^2)*phi^2 + (-4*a^3*d2*d3*e - 2*a^2*b*d1*d2*e)*phi -
8*a^4*c*d4*e + 4*a^4*d3^2*e + 2*a^3*b^2*d4*e - 4*a^3*b*d1*d3*e +
2*a^3*b*d2^2*e - 2*a^3*b*d4^2 + 4*a^3*c*d1^2*e ]);
assert [iota(x): x in Eltseq(BQ)] eq Eltseq(quartic);

A,B := Explode(Matrices(f));
R2<x,z> := PolynomialRing(K,2);
f1 := GenusOneModel((1/4)*Det(A*x + B*z));
E := BaseChange(Jacobian(f1),FF);
C,E,pi := nCovering(BQ : E := E);

P<X> := PolynomialRing(FF);
LL<beta> := quo<P|X^2 + phi^2 + 2*a*b>;
assert Evaluate(Equation(C),[1,0,beta]) eq 0;
eqns := DefiningEquations(pi);
P := E(LL)![Evaluate(e,[1,0,beta]): e in eqns];
assert P[1] eq xP;
xx := iota(P[1]);
yy := &+[iota(Eltseq(P[2])[j])*[1,-a*(theta-theta1)][j] : j in [1..2]];
assert Evaluate(g,[xx,yy]) eq 0 or Evaluate(g,[xx,-yy]) eq 0;
