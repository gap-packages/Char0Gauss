#
# Char0Gauss_SolutionMat: Linear algebra stuff
#
# Implementations
#
#
DeclareInfoClass("InfoChar0GaussChineseRem");
Char0Gauss_SolutionMat := function(A, b, moduli, die)
    local
        Ap, bp,
        solns, soln_i,
        prod, step, segment, segment_size, segment_i, segment_lcm_den,
        no_used,
        indices, indices_lookahead,
        v, v_solvables,
        mtx_echel, mtx_soln, mtx_solvables, mtx_solvables_i, solvables, c_i,
        r, r_start, power, alphas, result, p,
        i, j, k, t,
        get_random_inds,
        get_next_prime,
        get_next_segment,
        get_segment_den_lcm,
        postpone_this_prime,
        bitstring_to_list,
        fill_with_junk,
        tmp;
    Info(InfoChar0GaussChineseRem, 10, "moduli ", moduli, "\n");
    # scan the whole matrix to choose the moduli and find the total denominator

    # some functions
    get_random_inds := function(l, k, len)
        Info(InfoChar0GaussChineseRem, 10, "get_random_inds ", l, "..", k, " ", len, "\n");
        return SortedList(Shuffle([l..k]){[1..len]});
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
    fill_with_junk := function(vec, fixed)
        return List([1..Length(vec)], function(ind)
            local bound;
            if ind in fixed then
                return vec[ind];
            else
                # bound := 9000;
                return 0; # Random([-bound..bound]) / Random([1..bound]);
            fi;
        end);
    end;
    bitstring_to_list := function(bitstring)
        local len, lst;
        lst := List([1..MTX64_LengthOfBitString(bitstring)], i -> i * MTX64_GetEntryOfBitString(bitstring, i - 1));
        lst := Set(lst);
        RemoveSet(lst, 0);
        return SortedList(lst);
    end;
    # steps := number_of_moduli_between
    step := 1;
    segment_i := 1;
    segment_size := 512;
    segment := moduli;
    # debug info
    Info(InfoChar0GaussChineseRem, 10, "SOLVING\n", A, "\nWITH\n", b, "\n");

    soln_i := 1; # has to be 1
    i := soln_i;
    r_start := soln_i;
    solns := []; # must be []
    r := List([1..Length(A)], x -> []); # remainders
    no_used := Maximum(1, Int(Length(moduli) * .5)); # number of variables used to reconstruct a solution
    indices := []; # indices of the solution moduli
    result := List([1..Length(A)], x -> fail);
    solvables := [];
    # iterate the array of moduli
    # since solving the matrix is the bottleneck, checks will occur
    # after each prime, rather than incremental chunks
    while true do
        power := 1;
#        # make a different subset
        if i > Length(moduli) then
            if die then
                break;
            fi;
            # reshuffle
            if Length(moduli) > 0 and RandomList([1..2]) < 2 then
#                Info(InfoChar0GaussChineseRem, 10, "ROUTINE SHUFFLE ", no_used, "/", Length(moduli), "\n");
                indices := get_random_inds(r_start, soln_i, no_used);
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
        # detect a problem a priori
        if GcdInt(2 * 3 * 5 * 7 * 11 * segment_lcm_den, p) <> 1 then
            postpone_this_prime(p);
            continue;
        fi;
        # initialize modulo p
        if IsPrime(p) then
            Ap := A * One(GF(p));
            bp := b * One(GF(p));
        else
            Ap := A * One(ZmodnZ(p));
            bp := b * One(ZmodnZ(p));
        fi;
        # solve the system modulo p
        if IsPrime(p) and p < 2^60 then
            # meataxe64 version if available
            Ap := MTX64_Matrix(Ap);
            bp := MTX64_Matrix([bp]);

            mtx_echel := MTX64_Echelize(Ap);
            v_solvables := MTX64_EmptyBitString(Length(A));
#            Display(-MTX64_ExtractMatrix(mtx_echel.remnant)[1]);
            for c_i in [0..MTX64_NumRows(mtx_echel.remnant) - 1] do
                if MTX64_DNzl(mtx_echel.remnant, c_i) = fail then
                    MTX64_SetEntryOfBitString(v_solvables, c_i);
                fi;
            od;
            v_solvables := bitstring_to_list(v_solvables);
            Info(InfoChar0GaussChineseRem, 10, "v_solvables: ", v_solvables, "\n");

            if not IsSubset(v_solvables, solvables) then
                Info(InfoChar0GaussChineseRem, 10, "not a subset: ", solvables, "\n");
                v := fail;
            else
#                Print("mult:", -MTX64_ExtractMatrix(mtx_echel.multiplier)[1], "\n");
#                Print("remn:", -MTX64_ExtractMatrix(mtx_echel.remnant)[1], "\n");
                mtx_soln := MTX64_SolutionsMat(Ap, bp)[2];
                if MTX64_NumRows(mtx_soln) = 0 then
                    Info(InfoChar0GaussChineseRem, 10, "can't solve\n");
                    v := fail;
                else
                    tmp := -MTX64_ExtractMatrix(mtx_soln)[1];
                    v := List([1..Length(A)], function(ind)
                        local i;
                        for i in [1..Length(v_solvables)] do
                            if v_solvables[i] = ind then
                                return tmp[i];
                            fi;
                        od;
                        return fail;
                    end);
                fi;
#                Print("solution: ", v, "\n");
                if not v = fail and not v_solvables = solvables then
                    Info(InfoChar0GaussChineseRem, 10, "proper subset: ", solvables, "\n");
                    r := List([1..Length(A)], x -> []); # remainders
                    r_start := soln_i;
                    no_used := 0;
                    indices := [];
                fi;
            fi;
            solvables := SortedList(Union(Set(solvables), Set(v_solvables)));
#            Info(InfoChar0GaussChineseRem, 10, MTX64_LengthOfBitString(mtx_soln[1]), ", " , MTX64_NumRows(mtx_soln[2]), 'x', MTX64_NumCols(mtx_soln[2]), "\n");
        # else
            # standard version if available
            # v := SolutionMat(Ap, bp);
        fi;
#        # if the failure is detected a posteriori, choose a different prime
        if v = fail then
            postpone_this_prime(p);
            continue;
        else # prime approved
            Info(InfoChar0GaussChineseRem, 10, "approved: prime ", p, "; solution ",
                 List(v, function(x)
                    if x = fail then
                        return fail;
                    else
                        return Int(x);
                    fi;
                end),
            "\n");
            Append(solns, [v]);
            # see if using previous indices + this prime we succeed
            for k in [1..Length(v)] do
                if k in v_solvables then
                    r[k][soln_i - r_start + 1] := Int(v[k]);
                else
                    r[k][soln_i - r_start + 1] := fail;
                fi;
            od;
            indices_lookahead := indices;
            Append(indices_lookahead, [i]);

            prod := Product(moduli{indices_lookahead});
            alphas := List(moduli{indices_lookahead}, m -> Gcdex(prod / m, m).coeff1 mod m);
#            Info(InfoChar0GaussChineseRem, 10, "ALPHAS: ", List([1..Length(indices_lookahead)], x -> [moduli[x], alphas[x]]), "\n");
            for k in [1..Length(v)] do
#                result[k] := RationalReconstruction(prod, ChineseRem(moduli{indices_lookahead}, r[k]{indices_lookahead}));
                if k in v_solvables then
                    result[k] := Sum(
                        List([1..Length(indices_lookahead)], function(ind)
                            # Print("r_ind: ", indices_lookahead, ":", r_start, ":", soln_i, " -> ", List(indices_lookahead, x -> x - r_start), "\n");
                            return alphas[ind] * r[k][indices_lookahead[ind] - r_start + 1] * prod
                                  / moduli[indices_lookahead[ind]];
                        end)
                    ) mod prod;
                    result[k] := RationalReconstruction(prod, result[k]);
                    if result[k] = fail then
                        break;
                    fi;
                else
                    result[k] := fail;
                fi;
#                Print("CHREM: ", ChineseRem(moduli{indices_lookahead}, r[k]{indices_lookahead}), " VS ", t, "\n");
#                Info(InfoChar0GaussChineseRem, 10, k, "th - ", ChineseRem(moduli{indices_lookahead}, r[k]{indices_lookahead}), " -> ", result[k], "\n");
                if k = Length(v) then
#                    Print("result:", result{solvables}, "\n");
                    if fill_with_junk(result, solvables) * List(A, vec -> fill_with_junk(vec, solvables)) = b then
                        Info(InfoChar0GaussChineseRem, 10, r, "\n");
                        Info(InfoChar0GaussChineseRem, 10, "Primes: ", Length(indices_lookahead), "/", Length(moduli), "\n");
                        Info(InfoChar0GaussChineseRem, 10, "Result: ", result, "\n");
                        return result;
                    fi;
                fi;
            od;
            i := i + 1;
        fi;
        soln_i := i;
        segment_i := segment_i + 1;
        no_used := Int(.90 * (soln_i - r_start + 1));
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
#
#genmat_random := function(m, n)
#    return [RandomMat(m, n, Rationals), RandomMat(1, n, Rationals)[1]];
#end;
#
#genmat_invertible := function(m, n)
#    local A, b;
#    if m * n < 2^12 then
#        if m = n then
#            Display("generating invertible matrix...");
#            A := RandomInvertibleMat(n, Rationals);
#        else
#            repeat
#                Display("attempting to generate a matrix...");
#                A := RandomMat(m, n, Rationals);
#            until Rank(A) = Minimum(m, n);
#        fi;
#    else
#        A := RandomMat(m, n, Rationals);
#    fi;
#    Display("done");
#    b := RandomMat(1, n, Rationals)[1];
#    return [A, b];
#end;
#
#genmat_invertible_small_sln := function(m, n)
#    local A, v, bound;
#    if m * n < 2^12 then
#        if m = n then
#            Display("generating invertible matrix...");
#            A := RandomInvertibleMat(n, Rationals);
#        else
#            repeat
#                Display("attempting to generate a matrix...");
#                A := RandomMat(m, n, Rationals);
#            until Rank(A) = Minimum(m, n);
#            Display("done");
#        fi;
#    else
#        A := RandomMat(m, n, Rationals);
#    fi;
#    bound := (m*n)^2;
#    v := List([1..n], x -> Random([-bound..bound]) / Random([1..bound]));
#    return [A, v * A];
#end;
#
#
#run_tests_char0gauss := function(m, n, gen_mat, times)
#    local obj, A, b, v, w, u, i, times1, times2, times3, t1, t2, t3, small_primes;
#    times1 := [];
#    times2 := [];
#    times3 := [];
#    for i in [1..times] do
#        obj := gen_mat(m, n);
##        Display("generated matrices");
#        A := obj[1];
#        b := obj[2];
#
#        t2 := NanosecondsSinceEpoch();
#        w := Char0Gauss_SolutionMat(A, b, [], false);
#        t2 := NanosecondsSinceEpoch() - t2;
#
##        t3 := NanosecondsSinceEpoch();
##        u := C0GAUSS_SolutionMatVec_Padic(A, [b]);
##        t3 := NanosecondsSinceEpoch() - t3;
#
#        t1 := NanosecondsSinceEpoch();
#        v := SolutionMat(A, b);
#        t1 := NanosecondsSinceEpoch() - t1;
#
#        Append(times1, [t1 / 1000000000.]);
#        Append(times2, [t2 / 1000000000.]);
##        Append(times3, [t3 / 1000000000.]);
#        if not v = w then
#            Display("fail");
#            Display(v);
#            Display(w);
#            return false;
#        fi;
#    od;
#    Print("Standard\t", Int((m*n)^.5), "\t", Average(times1), "\n");
#    Print("ChineseR\t", Int((m*n)^.5), "\t", Average(times2), "\n");
##    Print("PadicSln\t", Int((m*n)^.5), "\t", Average(times3), "\n");
#    return true;
#end;
#
#SetInfoLevel(InfoChar0GaussChineseRem, 10);

#run_tests_char0gauss(2, 3, genmat_invertible, 100);
#for i in [1..20] do
#    if run_tests_char0gauss(i, i, genmat_invertible, 1000) = fail then
#        break;
#    fi;
#od;
