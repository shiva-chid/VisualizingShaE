
AttachSpec("spec");

n := 2;

Es := [];
// A random example from lmfdb
E := EllipticCurve([1, 0, 0, -3507, 6507]); // 462.g2
Append(~Es,E);
E := EllipticCurve([1, 0, 0, -130000, -18051943]); // 195.a1
Append(~Es,E);
// Examples with no 2-torsion
E := EllipticCurve([0, -1, 1, -929, -10595]); // 571.b1 // good
Append(~Es,E);
E := EllipticCurve([1, 0, 1, -60582, -5745352]); // 1058.b1 // problem
Append(~Es,E);
E := EllipticCurve([1, 0, 1, 253, -26862]); // 1058.b2 // problem
Append(~Es,E);
E := EllipticCurve([0, -1, 1, -208, -307]); // 1325.f1 // problem
Append(~Es,E);
E := EllipticCurve([0, -1, 1, -92, -311]); // 1613.b1 // good
Append(~Es,E);
E := EllipticCurve([0, 0, 1, -243, -1519]); // 1701.j1 // good
Append(~Es,E);

for E in Es do
A, phi := TwoSelmerGroup(E);
Sel2elts := [a @@ phi : a in A];
alltwocovers := [TwoCover(x) : x in Sel2elts];
allgenusonemodels := [gx : x in alltwocovers | not IsEquivalent(gx,GenusOneModel(2,E)) where gx := GenusOneModel(x)];
#allgenusonemodels;

allsha2 := [GenusOneModel(2,E)] cat allgenusonemodels;
assert forall{x : x in allsha2 | IsLocallySoluble(x)};
// allsha2 := [RandomTransformation(2 : Unimodular := true, Size := 5)*x : x in allsha2];
allsha2 := [Minimise(Reduce(x)) : x in allsha2];
for psi in allsha2 do
    P<x,z> := Parent(DefiningEquation(psi));
end for;
allsha2;
allsha2abcdes := [Eltseq(psi) : psi in allsha2];
allsha2abcdes;

// newabcdes := [[-y/x[1] : y in x] : x in allsha2abcdes | x[1] ne 0];
// newsha2 := [GenusOneModel(2,x) : x in newabcdes];

///////////////////////////////////////////////////


for psi in allsha2[2..#allsha2] do
// for psi in newsha2 do
    abcde := Eltseq(psi);
//    E1, E2, C := EllipticCurvePairsFromShaElement(abcde);
    E1, E2, C := EllipticCurvePairsFromShaElement(abcde : direct := true);
    if Type(E1) ne CrvEll then
        printf "Could not construct elliptic curves from Sha element\n";
        continue;
    end if;
/*
    _, E2old, _ := EllipticCurvePairsFromShaElement(abcde);
    assert IsIsomorphic(E2,E2old);
*/
    assert IsIsomorphic(E1,E);
//    assert IsIsomorphic(E1,E) or IsQuadraticTwist(E1,E);
    print RankBounds(Jacobian(C));
//    E2 := QuadraticTwist(E2,-abcde[1]);
    rk, boo := Rank(E2);
    tors := AbelianInvariants(TorsionSubgroup(E2));
    printf "j(E1): %o\nCond(E2): %o, Rank(E2) = %o, TorsionSubgroup(E2) = %o\n\n", jInvariant(E1), Conductor(E2), <rk,boo>, tors;
//    printf "RichelotIsogenous surfaces of Jac(C): \n%o\n", RichelotIsogenousSurfaces(C);
end for;

end for;
