#
# Char0Gauss: Linear algebra stuff
#
# Implementations
#

LoadPackage("meataxe64");

#
# Exists as "RatNumberFromModular" in the edim package as well.
#
RationalReconstruction := function(N, r)
    local a0, b0, a1, b1, a2, b2, q;

    a0 := N;
    b0 := 0;
    a1 := r;
    b1 := 1;

    while 2 * a1^2 > N - 1 do
        q := QuoInt(a0, a1);

        a2 := a0 - q * a1;
        b2 := b0 - q * b1;

        a0 := a1; b0 := b1;
        a1 := a2; b1 := b2;
    od;
    if (2 * b1^2 <= N - 1) and GcdInt(a1, b1) = 1 then
        return a1 / b1;
    fi;

    return fail;
end;

## testing the thing
#Display(RationalReconstruction(101, 77));
#for n in [10..1000] do
#    count := 0;
#    for i in [0..n-1] do
#        q := RationalReconstruction(n, i);
#        if q = fail then
#        else
#    #        Print(count, " : ", i, ": ", q, "\n");
#            count := count + 1;
#        fi;
#    od;
#    ratio := 1. * (count / n);
#    Print("[0..", n, ") size: ", count, " ratio: ", ratio, "\n");
#od;

ErrorTolerantLifting := function(N, r)
    local a0, b0, a1, b1, q;
    a0 := N;
    b0 := 0;
    a1 := r;
    b1 := 1;
    repeat
        q := [
            a0 - QuoInt(LCM_INT(a0, b0), a1^2) * a1,
            b0 - QuoInt(LCM_INT(a1, b1), b1^2) * b1
        ];
        a0 := a1;
        b0 := b1;
        a1 := q[1];
        b1 := q[2];
    until a1^2 + b1^2 >= a0^2 + b0^2;
    if a0^2 + b0^2 < N then
        return a0/b0;
    fi;
    return fail;
end;

Char0Gauss := function(A, b)
    local s, nom, den, solns, primes, v, p, i, j, t, prod;
    # scan the whole matrix to choose the primes and find the total denominator
    nom := 1;
    den := 1;
    primes := Set([]);
    # scanning A for primes in numerators and denominators separately
    for i in [1..Length(A)] do
        for j in [1..Length(A[i])] do
            if A[i][j] <> 0 then
                nom := LCM_INT(nom, NumeratorRat(A[i][j]));
                den := LCM_INT(den, DenominatorRat(A[i][j]));
            fi;
        od;
    od;
    # scanning b for primes in numerators and denominators separately
    for i in [1..Length(b)] do
        nom := LCM_INT(nom, NumeratorRat(b[i]));
        den := LCM_INT(den, DenominatorRat(b[i]));
    od;
    Print("lcm numerator: ", nom, "\n");
    Print("lcm denominator: ", den, "\n");
    # add primes until prod becomes 0
    prod := (nom * den)^2;
    primes := [2];
    while prod > 0 do
        p := NextPrimeInt(NextPrimeInt(primes[Length(primes)]));
        prod := QuoInt(prod, p);
        primes[Length(primes) + 1] := p;
        Print("increasing boundaries with ", p, "\n");
    od;
#    primes[Length(primes) + 1] := 756839;
    # solution vector
    solns := [];
    A := den * A;
    b := den * b;
    i := 1;
    Display(primes);
    # debug info
    Display("SOLVING");
    Display(A);
    Display("WITH");
    Display(b);
    # iterate the array of primes
    repeat
        p := primes[i];
        # solve the equation modulo ith prime
        v := SolutionMat(One(Z(p)) * A, One(Z(p)) * b);
        # if the failure is detected a posteriori, choose a different prime
        if v = fail or p in FactorsInt(den) then
            Print("skipped: prime ", p);
            # find new prime that will replace p and hope that it gets solved for it.
            while true do
                p := NextPrimeInt(p);
                if i < Length(primes) and not (p in primes) or (i = Length(primes)) then
                    break;
                fi;
            od;
            primes[i] := p;
            Sort(primes);
            Print("; added ", primes[i], " instead\n");
        else
            # otherwise, proceed
            Print("approved: prime ", p, "; solution ", List(v, x -> Int(x)), "\n");
            i := i + 1;
            Append(solns, [v]);
        fi;
    until i > Length(primes);
    if Length(solns) = 0 then
        return fail;
    fi;
    v := 0 * [1..Length(solns[1])];
    Display(primes);
    Display(solns);
    for i in [1..Length(solns[1])] do
        v[i] := 0 * [1..Length(solns)];
        for j in [1..Length(primes)] do
          v[i][j] := Int(solns[j][i]);
        od;
        Print("remainders: ", v[i], "\n");
        v[i] := RationalReconstruction(Product(primes), ChineseRem(primes, v[i]));
    od;
    return v;
end;

#primes := [ 3, 5, 11, 19, 23, 31, 41, 47, 59, 71, 73, 83, 97, 103 ];
#Print("***********************************************\n");
#Display(primes);
#for r in [
#    List([1..Length(primes)], i -> Int(One(Z(primes[i])) * 2)),
#    List([1..Length(primes)], i -> Int(One(Z(primes[i])) * 87/182)),
#    List([1..Length(primes)], i -> Int(One(Z(primes[i])) * -4784805/416143))
#] do
#    Print("r = ", r, "\n");
#    r := ChineseRem(primes, r);
#    Print("modulus ", Product(primes), ": ", r, "\n");
#    Display(RationalReconstruction(Product(primes), r));
#od;
#Print("***********************************************\n");

A := [
  [1/2, 2/91, 391/2],
  [0, 2, 5],
  [0, 0, 4573/134]
];
b := [1, 1, 1];
Display(A);
Display(b);
Display(SolutionMat(A, b));
Display(Char0Gauss(A, b));
