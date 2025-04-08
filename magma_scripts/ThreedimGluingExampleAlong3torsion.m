E := EllipticCurve([0, -1, 1, -13934, -614414]); E := WeierstrassModel(E);
E;
a := -18058896;
b := -28882792944;
J := jInvariant(E)/1728;
alpha, beta := RubinSilverbergPolynomials(3,J);


for d in [1..10] do
    for n in [-10..10] do
        if GCD(n,d) gt 1 then continue; end if;
        Enew := EllipticCurve([a*Evaluate(alpha,n/d),b*Evaluate(beta,n/d)]);
        print Enew;
        if AnalyticRank(Enew) eq 0 then
            print n, d;
            // break;
        end if;
    end for;
end for;
