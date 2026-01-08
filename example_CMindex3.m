AttachSpec("spec");
L := [1, 1, 3, 5, -7, 15, 2, 19, 1, -2];
C := CMcurveForIndex3Torsor(L);
assert Genus(C) eq 4;
p := PrimeWitnessForMinimality(C);
print p;
