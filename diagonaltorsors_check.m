SetLogFile("diagonaltorsors_check.log");
AttachSpec("spec");

fil := Open("diagonaltorsors.txt","r");
s := Gets(fil);
count := 0;
while not IsEof(s) do
    count +:= 1; if count mod 100 eq 0 then print count; end if;
    try
        C1 := GenusOneModel(3,s);
        L1 := Eltseq(C1);
    catch e;
        print s;
        print e;
        s := Gets(fil);
        continue;
    end try;
    C2 := CMcurveForIndex3Torsor(L1);
    assert Genus(C2) eq 1;
    PP2<x,y,z> := AmbientSpace(C2);
    K<w> := BaseRing(PP2);
    KK<ww> := NumberFieldExtra(DefiningPolynomial(K));
    phi := hom<K->KK|ww>;
    C2 := BaseChange(C2,phi);
    singptsC2 := SingularPoints(C2);
    JacC2 := EllipticCurve(C2,singptsC2[1]);
    tau2 := PeriodMatrix(HyperellipticCurve(JacC2));

    JacC1 := Jacobian(C1);
    tau1 := PeriodMatrix(HyperellipticCurve(JacC1));

    boo, M := IsIsogenousPeriodMatrices(tau1,tau2);
    assert boo and Determinant(M) in {3,-3};

    jJacC2 := jInvariant(JacC2);
    assert jJacC2 in Rationals();
    E2j := EllipticCurveFromjInvariant(jJacC2);
    boo, twist_fac := IsQuadraticTwist(JacC2,E2j);
    assert boo and twist_fac in Rationals();

    s := Gets(fil);
end while;
