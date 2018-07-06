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

#ErrorTolerantLifting := function(N, r)
#    local a0, b0, a1, b1, q;
#    a0 := N;
#    b0 := 0;
#    a1 := r;
#    b1 := 1;
#    repeat
#        q := [
#            a0 - QuoInt(LCM_INT(a0, b0), a1^2) * a1,
#            b0 - QuoInt(LCM_INT(a1, b1), b1^2) * b1
#        ];
#        a0 := a1;
#        b0 := b1;
#        a1 := q[1];
#        b1 := q[2];
#    until a1^2 + b1^2 >= a0^2 + b0^2;
#    if a0^2 + b0^2 < N then
#        return a0/b0;
#    fi;
#    return fail;
#end;
#
DeclareInfoClass("InfoChar0GaussChineseRem");
Char0Gauss := function(A, b, moduli, die)
    local
        Ap, bp,
        nom, den,
        solns, soln_i,
        prod, step,
        no_used,
        indices, indices_lookahead,
        v, mtx_soln, r, alphas, result, p,
        i, j, k, t,
        get_random_inds,
        get_next_prime,
        postpone_this_prime;
    Info(InfoChar0GaussChineseRem, 10, "moduli ", moduli, "\n");
    # scan the whole matrix to choose the moduli and find the total denominator
    nom := 1;
    den := 1;
    # scanning A for moduli in numerators and denominators separately
    for i in [1..Length(A)] do
        for j in [1..Length(A[i])] do
            if A[i][j] <> 0 then
                nom := LCM_INT(nom, NumeratorRat(A[i][j]));
                den := LCM_INT(den, DenominatorRat(A[i][j]));
            fi;
        od;
    od;
    # scanning b for moduli in numerators and denominators separately
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
    # steps := number_of_moduli_between
    step := 1;
    A := den * A;
    b := den * b;
    # debug info
    Info(InfoChar0GaussChineseRem, 10, "SOLVING\n", A, "\nWITH\n", b, "\n");

    soln_i := 1; # has to be 1
    solns := []; # must be []
    r := List([1..Length(A)], x -> []); # remainders
    no_used := Maximum(1, Int(Length(moduli) * .5)); # number of variables used to reconstruct a solution
    indices := []; # indices of the solution moduli
    result := List([1..Length(A)], x -> fail);
    # iterate the array of moduli
    # since solving the matrix is the bottleneck, checks will occur
    # after each prime, rather than incremental chunks
    while true do
        i := soln_i;
        # make a different subset
        if i > Length(moduli) then
            if die then
                break;
            elif Length(moduli) > 0 and RandomList([1..2]) < 2 then
                Info(InfoChar0GaussChineseRem, 10, "ROUTINE SHUFFLE ", no_used, "/", Length(moduli), "\n");
                indices := get_random_inds(no_used, Length(moduli));
            fi;
            # add another prime
            if Length(moduli) = 0 then
                p := 1;
            else
                p := Int(moduli[Length(moduli)] * 1.0);
            fi;
            moduli[Length(moduli) + 1] := get_next_prime(p, step);
            Info(InfoChar0GaussChineseRem, 10, "added a prime[", i, "] ", moduli[Length(moduli)], "\n");
        fi;
        p := moduli[i];
        # solve the equation modulo ith prime
        # need to check the dimension
        v := fail;
        postpone_this_prime := function(p)
            # find new prime that will replace p and hope that it gets solved for it.
            Info(InfoChar0GaussChineseRem, 10, "skipped: prime ", p, "\n");
            while true do
                p := NextPrimeInt(p);
                if i < Length(moduli) and not (p in moduli) or (i = Length(moduli)) then
                    break;
                fi;
            od;
            moduli[i] := p;
            Sort(moduli);
            Info(InfoChar0GaussChineseRem, 10, "added ", moduli[i], " instead\n");
        end;
#        if Rank(Ap) = Minimum(Length(A), Length(A[1])) then
        # initialize modulo p
        if IsPrime(p) then
            Ap := One(GF(p)) * A;
            bp := One(GF(p)) * b;
        else
            Ap := One(ZmodnZ(p)) * A;
            bp := One(ZmodnZ(p)) * b;
        fi;
        # detect a problem a priori
        if Gcd(2 * 3 * 5 * 7 * 11 * den, p) <> 1 then
            postpone_this_prime(p);
            continue;
        fi;
        # solve the system modulo p
        if IsPrime(p) and p < 2^60 then
            # meataxe64 version if available
            mtx_soln := MTX64_SolutionsMat(MTX64_Matrix(Ap), MTX64_Matrix([bp]));
            Info(InfoChar0GaussChineseRem, 10, MTX64_LengthOfBitString(mtx_soln[1]), ", " , MTX64_Matrix_NumRows(mtx_soln[2]), 'x', MTX64_Matrix_NumCols(mtx_soln[2]), "\n");
            if(MTX64_Matrix_NumRows(mtx_soln[2])) = 0 then
                v := fail;
            else
                v := -MTX64_ExtractMatrix(mtx_soln[2])[1];
            fi;
        else
            # standard version if available
            v := SolutionMat(Ap, bp);
        fi;
        # if the failure is detected a posteriori, choose a different prime
        if v = fail then
            postpone_this_prime(p);
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

            prod := Product(moduli{indices_lookahead});
            alphas := List(indices_lookahead, ind -> Gcdex(prod / moduli[ind], moduli[ind]).coeff1 mod moduli[ind]);
#            Info(InfoChar0GaussChineseRem, 10, "ALPHAS: ", List([1..Length(indices_lookahead)], x -> [moduli[x], alphas[x]]), "\n");
            for k in [1..Length(v)] do
#                result[k] := RationalReconstruction(prod, ChineseRem(moduli{indices_lookahead}, r[k]{indices_lookahead}));
                result[k] := Sum(List([1..Length(indices_lookahead)], ind -> alphas[ind] * r[k][indices_lookahead[ind]] * prod / moduli[indices_lookahead[ind]])) mod prod;
                result[k] := RationalReconstruction(prod, result[k]);
#                Print("CHREM: ", ChineseRem(moduli{indices_lookahead}, r[k]{indices_lookahead}), " VS ", t, "\n");
#                Info(InfoChar0GaussChineseRem, 10, k, "th - ", ChineseRem(moduli{indices_lookahead}, r[k]{indices_lookahead}), " -> ", result[k], "\n");
                if result[k] = fail then
                    break;
                fi;
                if k = Length(v) and result * A = b then
                    Info(InfoChar0GaussChineseRem, 10, r, "\n");
                    return result;
                fi;
            od;
            i := i + 1;
        fi;
        soln_i := i;
        no_used := Int(.95 * soln_i);
    od;
    return result;
end;

#moduli := [ 2, 5, 11, 17, 23, 31, 41, 47, 59, 67, 73, 83, 97, 103, 109, 127,
#137, 149, 157, 167, 179, 191, 197, 211, 227 ];
#Print("***********************************************\n");
#Display(moduli);
#for r in [
#    List([1..Length(moduli)], i -> Int(One(Z(moduli[i])) * 65536)),
##    List([1..Length(moduli)], i -> Int(One(Z(moduli[i])) * -130981/182)),
#    List([1..Length(moduli)], i -> Int(One(Z(moduli[i])) * -1622888073626375/4323946081))
#] do
#    Print("r = ", r, "\n");
#    r := ChineseRem(moduli, r);
#    Print("modulus ", Product(moduli), ": ", r, "\n");
#    Display(RationalReconstruction(Product(moduli), r));
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
#Display(Char0Gauss(A, b, moduli2{[2..10]}), false);

first_n_moduli := function(n)
    local moduli, p, i;
    moduli := [];
    for i in [1..n] do
        if i = 1 then
            p := 2;
        else
            p := NextPrimeInt(p);
        fi;
        Append(moduli, [p]);
    od;
    return moduli;
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

SetInfoLevel(InfoChar0GaussChineseRem, 10);

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

for i in [1..50] do
    if run_tests_char0gauss(i, i, 1) = fail then
        break;
    fi;
od;

#p := 17;
#p := p^2;
#v := List(SolutionMat(One(ZmodnZ(p)) * [[1]], One(ZmodnZ(p)) * [-1]), x -> Int(x));
#Display(v);
#Display(RationalReconstruction(p, ChineseRem([p], v)));

# TODO1: a chosen prime is in the denominator of the solution.
# TODO2: free variables and dimension of the kernel
# TODO3: quick finish fails because full reconstruction does not imply correct solution
# solved: wait for 6 verifications before approving
