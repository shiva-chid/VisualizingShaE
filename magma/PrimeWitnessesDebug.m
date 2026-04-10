/*PrimeWitnesses.m

  Use this script to generate prime witness data:

  ls data/sha_order3_processed | parallel -j30 "magma -b InputFileName:={} PrimeWitnesses.m"

  Run from magma/ directory.

*/

// Timeout in seconds for PrimeWitnessFromCoefficients. Set to 0 for no timeout.
PRIME_WITNESS_COMPUTATION_S := 0;

RealInputFileName := "../data/sha_order3_processed/" cat InputFileName;

AttachSpec("spec");

LinesOfInputFile := Split(Read(RealInputFileName), "\n");

function PrimeWitnessFromCoefficients(L)
    // Given a list L of 10 integers (ternary cubic coefficients),
    // compute the CM curve and return the prime witness for minimality.
    C := CMcurveForIndex3Torsor(L);
    p := PrimeWitnessForMinimality(C);
    return p;
end function;

function MagmaListFromProcessedLine(MyLine)
    /*
        Takes e.g. "5 6 8 -9 -16 27 -17 40 -7 18 "
        to the magma list of ints [5, 6, 8, -9, -16, 27, -17, 40, -7, 18]
    */
    parts := Split(MyLine, " ");
    L := [];
    for s in parts do
        if s ne "" then
            Append(~L, StringToInteger(s));
        end if;
    end for;
    return L;
end function;

for i -> MyLine in LinesOfInputFile do
    L := MagmaListFromProcessedLine(MyLine);
    printf "Processing line %o: %o\n", i, L;
    timed_out := false;
    p := 0;
    try
        p := PrimeWitnessFromCoefficients(L);
    catch e
        printf "FAILURE for input %o: %o\n", L, e;
    end try;
end for;

quit;