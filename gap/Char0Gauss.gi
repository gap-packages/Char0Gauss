#
# Char0Gauss: Linear algebra stuff
#
# Implementations
#

LoadPackage("meataxe64");
LoadPackage("gauss");

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

DeclareInfoClass("InfoChar0GaussChineseRem");
Char0Gauss := function(A, b, primes, die)
    local
        Ap, bp,
        nom, den,
        solns, soln_i,
        prod, step,
        no_used,
        indices, indices_lookahead,
        v, mtx_soln, r, result, p,
        i, j, k, t,
        finish, failures,
        get_random_inds,
        get_next_prime;
    Info(InfoChar0GaussChineseRem, 10, "primes ", primes, "\n");
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
#    Info(InfoChar0GaussChineseRem, 10, 10, "Hello, world");
    Info(InfoChar0GaussChineseRem, 10, "lcm numerator: ", nom, "\n");
    Info(InfoChar0GaussChineseRem, 10, "lcm denominator: ", den, "\n");

    # some functions
    get_random_inds := function(k, len)
        Info(InfoChar0GaussChineseRem, 10, "get_random_inds ", k, " ", len, "\n");
        return SortedList(Shuffle([1..len]){[1..k]});
    end;
    get_next_prime := function(n, step)
        local i, p;
        p := n;
        for i in [1..step] do
            p := NextPrimeInt(p);
        od;
        return p;
    end;
    # steps := number_of_primes_between
    step := 1;
    Info(InfoChar0GaussChineseRem, 10, "primes ", primes, "\n");
    if primes = [] then
        primes := [2];
        # add primes until prod becomes 0
#        prod := ((nom * den + 1) * 2)^2;
#        while prod > 0 do
##            step := step + 1;
#            p := get_next_prime(primes[Length(primes)], step);
#            prod := QuoInt(prod, p);
#            primes[Length(primes) + 1] := p;
#            Info(InfoChar0GaussChineseRem, 10, "increasing boundaries with ", p, "\n");
#        od;
    fi;
#    primes[Length(primes) + 1] := 756839;
    A := den * A;
    b := den * b;
    Info(InfoChar0GaussChineseRem, 10, "primes ", primes, "\n");
    # debug info
    Info(InfoChar0GaussChineseRem, 10, "SOLVING\n", A, "\nWITH\n", b, "\n");

    soln_i := 1; # has to be 1
    solns := []; # must be []
    r := List([1..Length(A)], x -> []); # remainders
    no_used := Maximum(1, Int(Length(primes) * .5)); # number of variables used to reconstruct a solution
    indices := []; # indices of the solution primes
    failures := 0; # how many failures left before we reshuffle
    result := List([1..Length(A)], x -> fail);
    # iterate the array of primes
    # since solving the matrix is the bottleneck, checks will occur
    # after each prime, rather than incremental chunks
    while true do
        i := soln_i;
        p := primes[i];
        # solve the equation modulo ith prime
        # need to check the dimension
        v := fail;
#        if Rank(Ap) = Minimum(Length(A), Length(A[1])) then
        # solve modulo p
        Ap := One(GF(p)) * A;
        bp := One(GF(p)) * b;
        if p < 2^60 then
#            Display(p);
#            finish := SolutionMat(Ap, bp);
            mtx_soln := MTX64_SolutionsMat(MTX64_Matrix(Ap), MTX64_Matrix([bp]));
            Info(InfoChar0GaussChineseRem, 10, MTX64_LengthOfBitString(mtx_soln[1]), ", " , MTX64_Matrix_NumRows(mtx_soln[2]), 'x', MTX64_Matrix_NumCols(mtx_soln[2]), "\n");
            if(MTX64_Matrix_NumRows(mtx_soln[2])) = 0 then
                v := fail;
            else
                v := List([1..Length(A)], ind -> -One(Z(p)) * MTX64_ExtractFieldElement(MTX64_GetEntry(mtx_soln[2], 0, ind - 1)));
#                Display(List(finish, x -> Int(x)));
#                Display(List(v, x -> Int(x)));
#                if not v = finish then
#                    return fail;
#                fi;
            fi;
        else
            v := SolutionMat(Ap, bp);
        fi;
        # if the failure is detected a posteriori, choose a different prime
        if v = fail or Gcd(2 * 3 * 5 * 7 * 11 * den, p) <> 1 then
            Info(InfoChar0GaussChineseRem, 10, "skipped: prime ", p);
            # find new prime that will replace p and hope that it gets solved for it.
            while true do
                p := NextPrimeInt(p);
                if i < Length(primes) and not (p in primes) or (i = Length(primes)) then
                    break;
                fi;
            od;
            primes[i] := p;
            Sort(primes);
            Info(InfoChar0GaussChineseRem, 10, "; added ", primes[i], " instead\n");
            continue;
        else # prime approved
            Info(InfoChar0GaussChineseRem, 10, "approved: prime ", p, "; solution ", List(v, x -> Int(x)), "\n");
            Append(solns, [v]);
            # see if using previous indices + this prime we succeed
            for k in [1..Length(v)] do
                r[k][i] := Int(v[k]);
            od;
            indices_lookahead := indices;
            Append(indices_lookahead, [i]);
            finish := true;
#            if 2*(nom * den)^2 <= Product(primes{indices_lookahead}) then
            if true then
                if failures >= 3 then
                    failures := 0;
                    # can't happen on 1st iteration since failures > 0
                    indices := get_random_inds(Minimum(no_used, Length(indices)), i - 1);
                    Info(InfoChar0GaussChineseRem, 10, "INLOOP RESHUFFLE ", no_used, "/", i,"\n");
                fi;
                for k in [1..Length(v)] do
                    result[k] := RationalReconstruction(Product(primes{indices_lookahead}), ChineseRem(primes{indices_lookahead}, r[k]{indices_lookahead}));
#                    Info(InfoChar0GaussChineseRem, 10, k, "th - ", ChineseRem(primes{indices_lookahead}, r[k]{indices_lookahead}), " -> ", result[k], "\n");
                    if result[k] = fail then
                        finish := false;
                        failures := failures + 1;
                        Info(InfoChar0GaussChineseRem, 10, "FAILURES ", failures, "\n");
                        break;
                    fi;
                od;
                if finish and result * A = b then
                    Info(InfoChar0GaussChineseRem, 10, r, "\n");
                    return result;
                fi;
            fi;
            i := i + 1;
        fi;
        soln_i := i;
        no_used := Int(.95 * soln_i);
        # since not successful, make a different subset
        if i > Length(primes) then
            if die then
                break;
            elif RandomList([1..2]) < 2 then
                Info(InfoChar0GaussChineseRem, 10, "ROUTINE SHUFFLE ", no_used, "/", Length(primes), "\n");
                indices := get_random_inds(no_used, Length(primes));
            fi;
            primes[Length(primes) + 1] := get_next_prime(Int(primes[Length(primes)] * 1.0), step);
            Info(InfoChar0GaussChineseRem, 10, "added a prime[", i, "] ", p, "\n");
        fi;
    od;
    return result;
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
#Display(Char0Gauss(A, b, Primes2{[2..10]}), false);

first_n_primes := function(n)
    local primes, p, i;
    primes := [];
    for i in [1..n] do
        if i = 1 then
            p := 2;
        else
            p := NextPrimeInt(p);
        fi;
        Append(primes, [p]);
    od;
    return primes;
end;

#Read("gap/mat15");
#mat := TransposedMat(ConvertSparseMatrixToMatrix(mat));
#Read("gap/vec15");
#vec := ConvertSparseMatrixToMatrix(vecs);
#vec := List(vec, ir -> ir[15]);


run_tests_char0gauss := function(m, n, times)
    local A, b, v, w, i, times1, times2, t1, t2, small_primes;
    times1 := [];
    times2 := [];
    for i in [1..times] do
        repeat
            A := RandomMat(m, n, Rationals);
        until Rank(A) = Minimum(m, n);
#        A := mat;
#        b := vec;
#        Display(A);
#        Display(Rank(A));
        b := RandomMat(1, n, Rationals)[1];
#        Display(b);

#        small_primes := first_n_primes(2000);
        t2 := NanosecondsSinceEpoch();
        w := Char0Gauss(A, b, [], false);
#        w := Char0Gauss(A, b, Filtered(Primes2, x -> x > 1000000){[2..50]}, false);
#        w := Char0Gauss(A, b, Filtered(small_primes, x -> x > 100), false);
#        w := Char0Gauss(A, b, small_primes, false);
        t2 := NanosecondsSinceEpoch() - t2;

        t1 := NanosecondsSinceEpoch();
        v := SolutionMat(A, b);
        t1 := NanosecondsSinceEpoch() - t1;
        if not v = w then
            Display("FUCKUP");
            Display(A);
            Display(b);
            Print("correct: ", v, "\n");
            Print("wrong:   ", w, "\n");
            Display("FAIL");
            return fail;
        fi;
        Append(times1, [t1 / 1000000000.]);
        Append(times2, [t2 / 1000000000.]);
    od;
    Print("Standard\t", Int((m*n)^.5), "\t", Average(times1), "\n");
    Print("ChRemThm\t", Int((m*n)^.5), "\t", Average(times2), "\n");
    return true;
end;

SetInfoLevel(InfoChar0GaussChineseRem, 1);
#t := NanosecondsSinceEpoch();
#Char0Gauss(RandomMat(512, 512, Rationals), RandomMat(1, 512, Rationals)[1], first_n_primes(2000), false);
#t := (NanosecondsSinceEpoch() - t) / 1000000000.;
#Display(t);

#A := [
#  [75843, 237842, 2311],
#  [764, 762634, 1923904],
#  [3, 4, 5]
#];
#b := [1, 1, 1];
#Display(Char0Gauss(A, b, [], false));

for i in [1, 2, 4, 8, 16, 32, 64, 128, 256] do
    if run_tests_char0gauss(i, i, 1) = fail then
        break;
    fi;
od;

# TODO1: a chosen prime is in the denominator of the solution.
# TODO2: free variables and dimension of the kernel
# TODO3: quick finish fails because full reconstruction does not imply correct solution
# solved: wait for 6 verifications before approving
