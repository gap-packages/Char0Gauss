
#
# Exists as "RatNumberFromModular" in the edim package as well.
#
InstallGlobalFunction( RationalReconstruction,
function(N, r)
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
end );

InstallGlobalFunction( MatRationalReconstruction,
                       {N, mat} -> List(mat, r -> List(r, e -> RationalReconstruction(N, e) ) ) );

# FIXME:
PadicLess := function(a,b)
    local fam, p, q_a, q_b, r_a, r_b, div;

    fam := FamilyObj(a);
    p := fam!.prime;

    r_a := p^(a![1]) * a![2];
    r_b := p^(b![1]) * b![2];
    div := p ^ (fam!.precision - 1);

    repeat
        q_a := QuoInt(r_a, div);
        q_b := QuoInt(r_b, div);
        if r_a < r_b then
            return true;
        elif r_a > r_b then
            return false;
        fi;
        r_a := QuoInt(r_a, div);
        r_b := QuoInt(r_b, div);
        div := div / p;
    until div = 1;
    return false;
end;


# FIXME: better way of detecting insufficient progress and abort?
PadicDenominator := function(number, max_iter)
    local n, thresh, tmp, big, little, bigf, littlef, biggest, fam,
          is_int;

    # Threshold where we consider something an integer
    # This should probably not be computed every time
    fam := FamilyObj(number);
    thresh := fam!.prime ^ QuoInt(fam!.precision, 2);

    Info(InfoChar0GaussPadics, 10, " n: ", number, "\n");

    is_int := function(n)
        local negn;
        negn := -n;
        return (n![2] < thresh) or
               (negn![2] < thresh);
    end;

    if is_int(number) then
        return 1;
    fi;

    little := number;
    littlef := 1;
    big := number;
    bigf := 1;

    n := 0;
    while n < max_iter do
        n := n + 1;

        tmp := little + big;
        Info(InfoChar0GaussPadics, 10
             , " lf: ", String(littlef, 16)
             , " bf: ", String(bigf, 16));
        Info(InfoChar0GaussPadics, 10
             , " little: ", little
             , " big:    ", big);

        if is_int(tmp) then
            Info(InfoChar0GaussPadics, 1, "PadicDenominator iterations: ", n);
            return bigf + littlef;
        fi;

        if PadicLess(tmp, little) then
            little := tmp;
            littlef := bigf + littlef;
        elif PadicLess(big, tmp) then
            big := tmp;
            bigf := bigf + littlef;
        else
            Info(InfoChar0GaussLinearEq, 1, "little <= tmp <= big: "
                 , little, " "
                 , tmp, " "
                 , big);
            Error("This shouldn't happen");
        fi;
    od;

    Info(InfoChar0GaussPadics, 1
         , " failed to compute denominator after ", n, " iterations, giving up");
    return fail;
end;

# Compute LCM of denominators of a list of p-adics
# TODO: do we have to do them all?
PadicDenominatorList := function(list, max_iter)
    local denom, old_denom, k, iter, found;

    found := false;
    old_denom := 1;
    denom := 1;
    k := 1;

    repeat
        # denom := PadicDenominator(old_denom * list[k], max_iter);
        denom := PadicDenominator(list[k], max_iter);

        if (denom <> fail) and (denom > 1) then
            found := true;
            old_denom := Lcm(old_denom, denom);
        fi;

        k := k + 1;
        Info(InfoChar0GaussLinearEq, 10, "current denominator: ", old_denom);
    until k > Length(list);

    if found then
        Info(InfoChar0GaussLinearEq, 10, "found denominator: ", old_denom);
        return old_denom;
    else
        Info(InfoChar0GaussLinearEq, 10, "failed to find");
        return fail;
    fi;
end;

