/*
  To run this code, do

      magma -bi GluingGenus2.m
*/

R<x> := PolynomialRing(Integers());
print "This script computes the gluing data for a genus 2 curve over Q.";
E := EllipticCurve([1, 1, 1, -352, -2689]);  // taken from the LMFDB, this has label 66.b1
A := -456219;
B := -118606410;

f := x^3 + A*x + B;

a := 1;
b := 2;
c := 3;

g := (x - a)*(x - b)*(x - c);

h := f*g;
C := HyperellipticCurve(h);
J := Jacobian(C);
desc := HeuristicEndomorphismAlgebra(C);  // this shows that J is simple

// we now glue E and J
// TODO

print "Finished.";
//exit;