intrinsic getcharpols(C :: CrvHyp : primesend := 500) -> SeqEnum
{returns a sequence of tuples <p,charpol_p> where charpol_p is the characteristic polynomial 
of Frobenius at p on the Tate module of the Jacobian of the given genus 2 hyperelliptic curve C,
and p ranges over the good primes NthPrime(N) for all N within the given bounds.}
    P<x> := PolynomialRing(Rationals());
    C := SimplifiedModel(C);
    f := HyperellipticPolynomials(C);
    primebounds := " [1.." cat IntegerToString(NthPrime(primesend)) cat "]";
    fstr := " \"" cat Sprint(f) cat "\"";
    filename := Sprintf("temp_%o", Getpid());
    cmd := "lpdata2 " cat filename cat fstr cat primebounds cat " 1";
    _ := System(cmd);
    newfilename := filename cat "_lpdata.txt";
    fil := Open(newfilename, "r");
    s := Gets(fil);
    s := Gets(fil);
    charpols := [];
    P<x> := PolynomialRing(Integers());
    while not IsEof(s) do
        s_split := [StringToInteger(x) : x in Split(s,",")];
        p := s_split[1]; ap := s_split[2]; bp := s_split[3];
        Append(~charpols,<p,P![p^2,p*ap,bp,ap,1]>);
        s := Gets(fil);
    end while;
    delete(fil);
    System("rm " cat newfilename);
    return charpols;
end intrinsic;

