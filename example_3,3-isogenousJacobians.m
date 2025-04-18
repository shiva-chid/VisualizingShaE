// searched lmfdb g2c_curves for analytic rank 2, good prime 3, non-maximal prime 3
// and eyeballed a mod-3 image allowing a (3,3) isogeny
// This is our F' in Fisher's notation.
lbl := "28900.a.57800.1";
cond := 28900;
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 1, 4, 0, -3], R![1, 1, 0, 1]);
// {* IsSquare(ChangeRing(EulerFactor(C,p),GF(3))) : p in PrimesUpTo(10000) | cond mod p ne 0 and p ne 3 *};
Attach("getcharpols.m");
L := getcharpols(C);
{* IsSquare(ChangeRing(x[2],GF(3))) : x in L | cond mod x[1] ne 0 and x[1] ne 3 *};
delete(C);


// finding the isogenous Jacobian from isogeny-complete database
// This is F in Fisher's notation.
isogclasslbl := "28900.a";
cond := 28900;
fil := Open("big_g2database_isogcomplete.txt","r"); // Remember that this database only has typical curves, i.e., End = Z.
s := Gets(fil);
print s;
all_data := [];
while not IsEof(s) do
    s_split := Split(s,":");
    N := StringToInteger(s_split[1]);
    label := s_split[2];
    if N eq cond and label eq isogclasslbl then
        Append(~all_data,s);
    end if;
    s := Gets(fil);
end while;
#all_data;
// 1
all_data;
// [ 28900:28900.a:1781550892102302750:[[58225,17710345,57800,-18091400]]:[[[0,-2,4,-1,-3],[1,1,0,1]]]:[3]:[[-2654975,6651955585,5167377800,-5598689318600],[58225,17710345,57800,-18091400]]:[[[85140,-319781,84841,306026,-903429,-185679,-8148],[1,1,0,1]],[[0,-2,4,-1,-3],[1,1,0,1]]]:[[0,3],[3,0]] ]


R<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(R![85140,-319781,84841,306026,-903429,-185679,-8148],R![1,1,0,1]);
cond := Conductor(C);
// {* IsSquare(ChangeRing(EulerFactor(C,p),GF(3))) : p in PrimesUpTo(10000) | cond mod p ne 0 and p ne 3 *};
L := getcharpols(C);
{* IsSquare(ChangeRing(x[2],GF(3))) : x in L | cond mod x[1] ne 0 and x[1] ne 3 *};


// Need to still ensure that Tamagawa numbers of F are coprime to 3.
// The Tamagawa number of F'=28900.a.57800.1 at 2 is 3.
Attach("Tamagawa/Tamagawa_pkg2.m");
[<p,TamagawaNumber(RegularModel(C,p))> : p in PrimeFactors(cond)];
// [ <2, 3>, <5, 1>, <17, 1> ]
// Sadly the Tamagawa number at 2 is not coprime to 3.

// In general, the roles of F and F' can be swapped, as long as we ensure that Fisher's conditions hold for F.
