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

DeclareInfoClass("InfoChar0GaussChineseRem");
Char0Gauss := function(A, b, moduli, die)
    local
        Ap, bp,
#        nom,
#        den, # factors,
        solns, soln_i,
        prod, step, segment, segment_size, segment_i, segment_lcm_den,
        no_used,
        indices, indices_lookahead,
        v, mtx_soln, r, power, alphas, result, p,
        i, j, k, t,
        get_random_inds,
        get_next_prime,
        get_next_segment,
        get_segment_den_lcm,
        postpone_this_prime;
    Info(InfoChar0GaussChineseRem, 10, "moduli ", moduli, "\n");
    # scan the whole matrix to choose the moduli and find the total denominator
#    factors := Set([]);
#    nom := 1;
#    den := 1;
    # scanning A for moduli in numerators and denominators separately
#    for i in [1..Length(A)] do
#        for j in [1..Length(A[i])] do
#            if A[i][j] <> 0 then
#                factors := Union(factors, Set(FactorsInt(NumeratorRat(A[i][j]))));
##                nom := LCM_INT(nom, NumeratorRat(A[i][j]));
#                den := LCM_INT(den, DenominatorRat(A[i][j]));
#            fi;
#        od;
#        Print("A... ", i, "\n");
#    od;
    # scanning b for moduli in numerators and denominators separately
#    for i in [1..Length(b)] do
#        if b[i] <> 0 then
#            factors := Union(factors, Set(FactorsInt(NumeratorRat(b[i]))));
##            nom := LCM_INT(nom, NumeratorRat(b[i]));
#            den := LCM_INT(den, DenominatorRat(b[i]));
#        fi;
#        Print("b.. ", i, "\n");
#    od;
#    factors := Set(Union(FactorsInt(nom), FactorsInt(den)));
#    factors := Set(Union(factors, FactorsInt(den)));
#    Info(InfoChar0GaussChineseRem, 10, 10, "Hello, world");
#    Info(InfoChar0GaussChineseRem, 10, "lcm numerator: ", nom, "\n");
#    Info(InfoChar0GaussChineseRem, 10, "lcm denominator: ", den, "\n");

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
    get_next_segment := function(n, size, step)
        local i, cur, segment;
        segment := 0 * [1..size];
        cur := n;
        for i in [1..size] do
            cur := get_next_prime(cur, step);
            segment[i] := cur;
        od;
        return segment;
    end;
    get_segment_den_lcm := function(segment)
        local i, j, n, g, prod, lcm_den, break_n;
        prod := Product(segment);
        lcm_den := 1;
        for i in [1..Length(A)] do
            for j in [1..Length(A[i])] do
                g := GcdInt(DenominatorRat(A[i][j]), prod);
                if g <> 1 then
                    prod := prod / g;
                    lcm_den := lcm_den * g;
                fi;
            od;
        od;
        for i in [1..Length(b)] do
            g := GcdInt(DenominatorRat(b[i]), prod);
            if g <> 1 then
                prod := prod / g;
                lcm_den := lcm_den * g;
            fi;
        od;
        return lcm_den;
    end;
    # steps := number_of_moduli_between
    step := 1;
    segment_i := 1;
    segment_size := 512;
    segment := moduli;
#    A := den * A;
#    b := den * b;
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
        power := 1;
        i := soln_i;
        # make a different subset
        if i > Length(moduli) then
            if die then
                break;
            fi;
            # reshuffle
            if Length(moduli) > 0 and RandomList([1..2]) < 2 then
                Info(InfoChar0GaussChineseRem, 10, "ROUTINE SHUFFLE ", no_used, "/", Length(moduli), "\n");
                indices := get_random_inds(no_used, Length(moduli));
            fi;
            # add another prime
            if Length(moduli) = 0 then
                p := 1;
            else
                p := Int(RootInt(moduli[Length(moduli)], power) * 1.0);
            fi;
#            Print("segment_i: ", segment_i, ", segment: ", segment, "\n");
            if segment_i > Length(segment) then
                segment := get_next_segment(p, segment_size, step);
                segment_i := 1;
                segment_lcm_den := get_segment_den_lcm(segment);
#                Print("new segment: ", segment, "\n");
#                Print("segment lcm: ", segment_lcm_den, "\n");
            fi;
            moduli[Length(moduli) + 1] := segment[segment_i]^power;
            Info(InfoChar0GaussChineseRem, 10, "added a prime[", i, "] ", moduli[Length(moduli)], "\n");
        fi;
        p := moduli[i];
        # solve the equation modulo ith prime
        # need to check the dimension
        v := fail;
        postpone_this_prime := function(p)
            # find new prime that will replace p and hope that it gets solved for it.
            Info(InfoChar0GaussChineseRem, 10, "skipped: prime ", p, "\n");
            p := RootInt(p, power);
            while true do
                p := NextPrimeInt(p);
                if i < Length(moduli) and not (p^power in moduli) or (i = Length(moduli)) then
                    break;
                fi;
            od;
            if p in segment then
                # move segment index to match p
                while segment[segment_i] <> p do
                  segment_i := segment_i + 1;
                od;
            else
                # segment doesn't contain p
                segment[segment_i] := p;
                Sort(segment);
                segment_lcm_den := get_segment_den_lcm(segment);
            fi;
            moduli[i] := segment[segment_i] ^ power;
            Sort(moduli);
            Info(InfoChar0GaussChineseRem, 10, "added ", moduli[i], " instead\n");
        end;
        # if Rank(Ap) = Minimum(Length(A), Length(A[1])) then
        # detect a problem a priori
        if GcdInt(2 * 3 * 5 * 7 * 11 * segment_lcm_den, p) <> 1 then
            postpone_this_prime(p);
            continue;
        fi;
        # initialize modulo p
        if IsPrime(p) then
            Ap := One(GF(p)) * A;
            bp := One(GF(p)) * b;
        else
            Ap := One(ZmodnZ(p)) * A;
            bp := One(ZmodnZ(p)) * b;
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
            alphas := List(moduli{indices_lookahead}, m -> Gcdex(prod / m, m).coeff1 mod m);
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
                    Info(InfoChar0GaussChineseRem, 10, "Primes: ", Length(indices_lookahead), "/", Length(moduli), "\n");
                    Info(InfoChar0GaussChineseRem, 10, "Result: ", result, "\n");
                    return result;
                fi;
            od;
            i := i + 1;
        fi;
        soln_i := i;
        segment_i := segment_i + 1;
        no_used := Int(.90 * soln_i);
    od;
    return result;
end;


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


genmat_invertible := function(m, n)
    local A, b;
    if m = n then
        A := RandomInvertibleMat(n, Rationals);
    else
        repeat
            A := RandomMat(m, n, Rationals);
        until Rank(A) = Minimum(m, n);
    fi;
    b := RandomMat(1, n, Rationals)[1];
    return [A, b];
end;

genmat_invertible_small_sln := function(m, n)
    local A, v, bound;
    if m * n < 2^12 then
        if m = n then
            Display("generating invertible matrix...");
            A := RandomInvertibleMat(n, Rationals);
        elif m * n < 2^12 then
            repeat
                Display("attempting to generate a matrix...");
                A := RandomMat(m, n, Rationals);
            until Rank(A) = Minimum(m, n);
            Display("done");
        fi;
    else
        A := RandomMat(m, n, Rationals);
    fi;
    bound := (m*n)^2;
    v := List([1..n], x -> Random([-bound..bound]) / Random([1..bound]));
    return [A, v * A];
end;


run_tests_char0gauss := function(m, n, gen_mat, times)
    local obj, A, b, v, w, i, times1, times2, t1, t2, small_primes;
    times1 := [];
    times2 := [];
    for i in [1..times] do
        obj := gen_mat(m, n);
        Display("generated matrices");
        A := obj[1];
        b := obj[2];

        t2 := NanosecondsSinceEpoch();
        w := Char0Gauss(A, b, [], false);
        t2 := NanosecondsSinceEpoch() - t2;

#        t1 := NanosecondsSinceEpoch();
#        v := SolutionMat(A, b);
#        t1 := NanosecondsSinceEpoch() - t1;
#        Append(times1, [t1 / 1000000000.]);
        Append(times2, [t2 / 1000000000.]);
    od;
#    Print("Standard\t", Int((m*n)^.5), "\t", Average(times1), "\n");
    Print("ChRemThm\t", Int((m*n)^.5), "\t", Average(times2), "\n");
    return true;
end;

SetInfoLevel(InfoChar0GaussChineseRem, 1);

for i in [128,256,512,1024] do
    if run_tests_char0gauss(i, i, genmat_invertible_small_sln, 3) = fail then
        break;
    fi;
od;
