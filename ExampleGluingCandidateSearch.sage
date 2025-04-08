# cd into gluing_utilities/src/

load("compatible_curves.py")
R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, -1, 2, -2, 1]), R([1]));
C
N = 997
ell = 3

L = compatible_curves(C,N,ell)
L
[x.conductor() for x in L[0] if x.conductor() % 3 != 0]
[x for x in L[0] if x.conductor() % 3 != 0]
[x.rank() for x in L[0] if x.conductor() % 3 != 0]
