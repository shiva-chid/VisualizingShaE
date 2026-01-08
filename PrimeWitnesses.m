/*PrimeWitnesses.m

  Use this script to generate prime witness data:

  ls data/sha_order3_processed | parallel -j30 "magma -b InputFileName:={} PrimeWitnesses.m"

  This will take the data files in `./data/curve_data_new_may_23`, and for each curve in each file,
  will generate a new file (appended with `with_torsion`) in the directory `torsion_data_out` containing
  the torsion primes. 

  In the pipeline for running the main algorithm on lots of curves, one would then use the `produce_final_output.py`
  script in the `helper_scripts` folder to combine the output of the `nonmaximal.py` algorithm with these torsion
  primes.

*/

// Timeout in seconds for PrimeWitnessFromCoefficients. Set to 0 for no timeout.
PRIME_WITNESS_COMPUTATION_S := 2;

RealInputFileName := "data/sha_order3_processed/" cat InputFileName;
OutputFileName := "data/sha_order3_processed_witnessed/" cat "witnessed." cat InputFileName;

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
    try
        if PRIME_WITNESS_COMPUTATION_S gt 0 then
            SetAlarm(PRIME_WITNESS_COMPUTATION_S);
        end if;
        p := PrimeWitnessFromCoefficients(L);
        if PRIME_WITNESS_COMPUTATION_S gt 0 then
            SetAlarm(0);  // Clear the alarm
        end if;
    catch e
        if PRIME_WITNESS_COMPUTATION_S gt 0 then
            SetAlarm(0);  // Clear the alarm
        end if;
        printf "FAILURE for input %o: %o\n", L, e;
        p := 0;
    end try;
    to_print := MyLine cat " PrimeWitness: " cat IntegerToString(p) cat "\n";
    fprintf OutputFileName, to_print;
end for;

quit;