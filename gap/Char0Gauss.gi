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

Char0Gauss := function(A, b, primes)
    local s, nom, den, solns, step, no_append, v, p, soln_i, i, j, t, prod, tries, Ap, reason;
    # scan the whole matrix to choose the primes and find the total denominator
    nom := 1;
    den := 1;
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
        if b[i] <> 0 then
            nom := LCM_INT(nom, NumeratorRat(b[i]));
            den := LCM_INT(den, DenominatorRat(b[i]));
        fi;
    od;
    Print("lcm numerator: ", nom, "\n");
    Print("lcm denominator: ", den, "\n");

    # steps := number_of_primes_between
    tries := 1;
    step := 1;
    if primes = [] then
        tries := 5;
        primes := [2];
        # add primes until prod becomes 0
        prod := ((nom * den + 1) * 2)^2;

        while prod > 0 do
            step := step + 1;
            p := primes[Length(primes)];
            for i in [1..step] do
                p := NextPrimeInt(p);
            od;
            prod := QuoInt(prod, p);
            primes[Length(primes) + 1] := p;
            Print("increasing boundaries with ", p, "\n");
        od;
    fi;
#    primes[Length(primes) + 1] := 756839;
    # solution vector
    Print("primes ", primes, "\n");
    # debug info
    A := den * A;
    b := den * b;
    Display("SOLVING");
    Display(A);
    Display("WITH");
    Display(b);

    soln_i := 1;
    solns := [];
    no_append := 1;
    repeat
        Print("primes: ", primes, "\n");
        # iterate the array of primes
        v := [];
        i := soln_i;
        repeat
            p := primes[i];
            # solve the equation modulo ith prime
            # need to check the dimension
            Ap := Zero(Z(p)) * A;
            v := fail;
#            if Rank(Ap) = Minimum(Length(A), Length(A[1])) then
            reason := "no reason";
            if true then
                v := SolutionMat(One(Z(p)) * A, One(Z(p)) * b);
                if v = fail then
                    reason := "unsolvable";
                fi;
            else
                reason := "free variable";
#                reason := StringPrint("free variables (rank = ", Rank(A), ")");
            fi;
            # if the failure is detected a posteriori, choose a different prime
            #  or v = (Zero(Z(p)) * [1..Length(v)])
            if v = fail then
#                if not v = fail then
#                    reason := "zero solution vector";
#                fi;
                Print("skipped: prime ", p, " because ", reason);
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
        soln_i := Length(primes) + 1;
        if Length(solns) = 0 then
            return fail;
        fi;
        v := 0 * [1..Length(solns[1])];
        prod := Product(primes);
        Display(primes);
        Display(solns);
        for i in [1..Length(solns[1])] do
            v[i] := 0 * [1..Length(solns)];
            for j in [1..Length(primes)] do
                v[i][j] := Int(solns[j][i]);
            od;
            Print("remainders: ", v[i], "\n");
            v[i] := RationalReconstruction(prod, ChineseRem(primes, v[i]));
        od;

        # add more primes before retrying
        step := step + 1;
        no_append := no_append + 1;
#        step := QuoInt(9 * step, 2);
        for i in [1..no_append] do
            p := primes[Length(primes)];
    #        prod := QuoInt(prod, p);
            for j in [1..step] do
                p := NextPrimeInt(p);
            od;
            primes[Length(primes) + 1] := p;
        od;
        Print("reached a solution ", v, "\n");
        tries := tries - 1;
    until not fail in v or tries = 0;
    return v;
end;

#primes := [ 2, 5, 11, 17, 23, 31, 41, 47, 59, 67, 73, 83, 97, 103, 109, 127,
#137, 149, 157, 167, 179, 191, 197, 211, 227 ];
#Print("***********************************************\n");
#Display(primes);
#for r in [
#    List([1..Length(primes)], i -> Int(One(Z(primes[i])) * 65536)),
##    List([1..Length(primes)], i -> Int(One(Z(primes[i])) * -130981/182)),
#    List([1..Length(primes)], i -> Int(One(Z(primes[i])) * -1622888073626375/4323946081))
#] do
#    Print("r = ", r, "\n");
#    r := ChineseRem(primes, r);
#    Print("modulus ", Product(primes), ": ", r, "\n");
#    Display(RationalReconstruction(Product(primes), r));
#od;
#Print("***********************************************\n");

#A := [
#  [1/65536, 2/91, 391/2],
#  [0, 2, 5],
#  [0, 0, 4573/134]
#];
#b := [1, 1, 54/1236473];
#Display(A);
#Display(b);
#Display(SolutionMat(A, b));
#Display(Char0Gauss(A, b, Primes2{[2..10]}));

run_tests_char0gauss := function(m, n, times)
    local A, b, v, w, i, t1, t2;
    for i in [1..times] do
        repeat
            A := RandomMat(m, n, Rationals);
        until Rank(A) = Minimum(m, n);
#        Display(A);
#        Display(Rank(A));
        b := RandomMat(1, n, Rationals)[1];
        t1 := NanosecondsSinceEpoch();
        v := SolutionMat(A, b);
        t1 := NanosecondsSinceEpoch() - t1;
        Display(v);
        t2 := NanosecondsSinceEpoch();
        w := Char0Gauss(A, b, []);
        t2 := NanosecondsSinceEpoch() - t2;
        Print("Standard: ", t1 / 1000000000., "\n");
        Print("ChRemThm: ", t2 / 1000000000., "\n");
        if not v = w then
            Display("FUCKUP");
            Display(A);
            Display(b);
            Display(v);
            Display(w);
        fi;
    od;
end;

run_tests_char0gauss(50, 50, 10);
# TODO1: a chosen prime is in the denominator of the solution.
# TODO2: free variables and dimension of the kernel
